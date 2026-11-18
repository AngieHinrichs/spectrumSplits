#!/usr/bin/env python3
"""Rank spectrumSplits.py splits by null margin and bootstrap support.

This is step 10 of parameter tuning: deciding which of the accepted splits
are worth interpreting.  It reads only files the search and the null already
produce, so it adds no new computation on the tree.

Two independent questions, because they fail independently:

  1. MARGIN.  Is this split's chi-squared larger than permutation noise for a
     subtree of THIS size?  A single --min_chi cannot answer that, because
     chi-squared scales with the counts in the table: a threshold calibrated
     on the whole tree is far too strict for a deep split, and one calibrated
     on a deep split is too permissive at the root.  So each split is scored
     against null replicates drawn from subtrees of comparable size, and what
     is reported is the ratio.  This is the check that justifies running the
     search at a low --min_chi and filtering afterwards.

     Getting that size right is the whole difficulty, and it is not the size
     of the subtree named in the log.  The search carves accepted splits out
     as it goes, so the table a split was tested against is its search root
     MINUS every split already accepted inside it.  That shrinks across
     iterations, so using the search root's full size systematically inflates
     the size -- and therefore the null ceiling -- for every late-found
     split, which penalises exactly the deep splits the low threshold was
     meant to reach.  On the 129k MTBC tree this understated the margin of
     node_44869 by a factor of 3.7 and moved 8 splits between tiers.
     searched_size() models the carving, and is validated against the exact
     table totals the cascade null records.

  2. SUPPORT.  Does the split survive resampling the alignment positions?
     --bootstrap_splits writes one spectra file per replicate; a split that
     depends on a handful of recurrent sites will vanish from most of them
     even when its chi-squared is enormous.

A split can pass one and fail the other, and the two failures mean different
things.  A big margin with low support says the signal is real but rests on
few positions -- the artifact pattern.  A small margin with high support says
a consistent but weak difference.  Neither is reported as a single number.

THE PLACEMENT WOBBLE.  Bootstrap replicates often place a real split one node
up or down, so counting exact node IDs understates support.  Given
--bootstrap_replicates (per-replicate membership) this is handled by asking,
for each replicate, whether ANY node within --neighbourhood of the observed
node was split, which is a fraction of replicates and so cannot exceed 100%.
Summing the per-node counts instead would double-count replicates and can
exceed the number of replicates; that is an upper bound, not a support value,
and this script will not print it.  Note also that numeric ID adjacency is a
proxy for topological adjacency, not the same thing -- usher assigns IDs at
read time in traversal order, so neighbours are usually but not always close.

Usage:
    split_confidence.py --splits_log splits.log --spectra spectra.tsv \
        --null null_sites_cascade.tsv \
        --bootstrap_replicates bootstrap_by_replicate.tsv --output conf.tsv

    ### the per-replicate file is made from the replicate outputs with:
    ### for f in bootstrap_*_splits_output.tsv; do n=${f#bootstrap_};
    ###   n=${n%%_*}; cut -f1 "$f" | grep ^node_ | sed "s/^/$n\t/"; done
"""

import argparse
import collections
import re
import sys

### The search prints its own progress to stderr and then each bootstrap
### replicate prints to the SAME stream, so a redirected log holds both.  Only
### the lines before this marker describe the real tree.
BOOTSTRAP_MARKER = "Bootstrapping splits with"

RE_SUBTREE = re.compile(
    r"Computing distances between splits in \(sub\)tree (node_\S+)")
RE_SPLIT = re.compile(r"New split found at (node_\S+) with x2 ([\d.eE+-]+)")


def parse_splits_log(lines):
    """Return (chi, parent) for the real run's splits, ignoring bootstraps.

    The log interleaves, per iteration, a "(sub)tree node_P" line followed by
    the splits found inside node_P.  That is how a split's parent -- the
    subtree whose table it was tested against, and so the size its null must
    match -- is recovered; the TSV does not record it.
    """
    chi, parent, order = {}, {}, []
    current = None
    for line in lines:
        if line.startswith(BOOTSTRAP_MARKER):
            break
        m = RE_SUBTREE.match(line)
        if m:
            current = m.group(1)
            continue
        m = RE_SPLIT.match(line)
        if m:
            ### A split seen twice would mean the log was concatenated from
            ### two runs; that silently corrupts every total below it.
            if m.group(1) in chi:
                raise ValueError(f"{m.group(1)} appears twice before the "
                                 f"bootstrap marker; is this log two runs?")
            chi[m.group(1)] = float(m.group(2))
            parent[m.group(1)] = current
            ### Log order is the order the search accepted them, which is what
            ### says how much had already been carved out of each search root.
            order.append((m.group(1), current))
    return chi, parent, order


def parse_spectra(lines):
    """Node -> RESIDUAL mutation count, as spectrumSplits.py reports it.

    Each row covers only what is left of a subtree after its children were
    carved out, so these do not nest; cumulative_mutations() adds them back up.
    """
    res = {}
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if i == 0 or len(f) < 4 or not f[0]:
            continue
        res[f[0]] = int(f[1])
    return res


def read_clade_samples(lines):
    """'node<TAB>sample' per line -> node: set of samples."""
    sets = collections.defaultdict(set)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] and not f[0].startswith("#"):
            sets[f[0]].add(f[1])
    return dict(sets)


def nesting_from_samples(sets):
    """node -> set of splits strictly inside it, by sample membership.

    Membership is the only durable statement of which clade contains which,
    and it does not agree with the log: a split accepted late can turn out to
    be an ANCESTOR of one accepted earlier, since both were found while
    searching the same root.  On the MTBC tree that happened for 11 of 62.
    """
    nodes = list(sets)
    inside = {n: {m for m in nodes if m != n and sets[m] < sets[n]}
              for n in nodes}
    return inside


def nesting_from_parents(parent):
    """Fallback nesting from the log's search roots alone.

    Only used without --clade_samples.  It cannot see a late split that
    contains an earlier one, so the sizes it yields are upper bounds.
    """
    kids = collections.defaultdict(list)
    for child, par in parent.items():
        kids[par].append(child)
    inside = {}

    def walk(node, seen):
        if node in inside:
            return inside[node]
        if node in seen:
            raise ValueError(f"cycle in split parentage at {node}")
        seen = seen | {node}
        out = set()
        for child in kids[node]:
            out.add(child)
            out |= walk(child, seen)
        inside[node] = out
        return out

    for node in set(parent) | set(parent.values()):
        if node is not None:
            walk(node, frozenset())
    return inside


def subtree_mutations(node, inside, residual):
    """Mutations in a node's whole subtree.

    spectrumSplits.py reports RESIDUAL counts -- each row covers only what is
    left of a subtree after its children were carved out -- so a subtree is
    its own residual plus the residuals of everything inside it.
    """
    return residual.get(node, 0) + sum(residual.get(m, 0)
                                       for m in inside.get(node, ()))


def searched_size(order, k, inside, residual, root="node_1"):
    """Mutations in the table the k'th split was tested against.

    Its search root, less every split already accepted strictly inside that
    root.  Only the MAXIMAL prior splits are subtracted: a prior split nested
    inside another prior split was already removed along with its ancestor,
    and subtracting it again would undercount the table.
    """
    node, search_root = order[k]
    if search_root is None:
        return None
    contained = inside.get(search_root)
    if contained is None:
        contained = set(inside) if search_root == root else set()
    prior = [p for p, _ in order[:k] if p in contained]
    maximal = [p for p in prior
               if not any(p in inside.get(q, ()) for q in prior if q != p)]
    return (subtree_mutations(search_root, inside, residual)
            - sum(subtree_mutations(m, inside, residual) for m in maximal))


def parse_null(lines):
    """(subtree_mutations, chi2) for NULL replicates only.

    The 'observed' rows are the real splits of whichever run produced the
    file; including them would compare observations against themselves.
    """
    out = []
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if i == 0 or len(f) < 7 or f[0] == "observed":
            continue
        out.append((int(f[5]) + int(f[6]), float(f[4])))
    return out


def ceiling(nulls, size, window):
    """Highest null chi-squared from subtrees within `window`-fold of `size`.

    A maximum, not a quantile: with a few dozen replicates in the window the
    upper tail is not resolved well enough for a quantile to mean anything,
    and the maximum is the conservative reading.
    """
    inside = [c for s, c in nulls if size / window <= s <= size * window]
    return (max(inside), len(inside)) if inside else (None, 0)


def parse_support(lines):
    """Aggregate `uniq -c` support counts: node -> replicates containing it."""
    sup = {}
    for line in lines:
        f = line.split()
        if len(f) == 2 and f[1].startswith("node_"):
            sup[f[1]] = int(f[0])
    return sup


def parse_replicates(lines):
    """replicate -> set of node IDs split in it."""
    reps = collections.defaultdict(set)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].startswith("node_"):
            reps[f[0]].add(f[1])
    return reps


def node_number(node):
    """Trailing integer of a node ID, or None when it has none."""
    m = re.search(r"(\d+)$", node)
    return int(m.group(1)) if m else None


def neighbourhood_support(reps, node, window):
    """Replicates splitting `node` or any node within `window` of its number.

    A fraction of replicates: each replicate counts once however many of its
    splits fall in the window, which is what keeps this from exceeding 100%.
    """
    n = node_number(node)
    if n is None or not reps:
        return None
    wanted = range(n - window, n + window + 1)
    hits = 0
    for nodes in reps.values():
        for other in nodes:
            k = node_number(other)
            if k is not None and k in wanted:
                hits += 1
                break
    return hits


def tier(margin, support, strong_margin, strong_support):
    """A: both criteria.  B: one.  C: neither.  Unknown inputs never pass."""
    ok_m = margin is not None and margin >= strong_margin
    ok_s = support is not None and support >= strong_support
    if ok_m and ok_s:
        return "A"
    return "B" if (ok_m or ok_s) else "C"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Rank splits by size-matched null margin and bootstrap "
                    "support.")
    p.add_argument("--splits_log", help="stderr log from the spectrumSplits.py "
                                        "run (bootstrap lines are ignored)")
    p.add_argument("--spectra", help="--output_spectrum TSV from the same run")
    p.add_argument("--null", help="cascade null TSV from permutation_null.py")
    p.add_argument("--bootstrap_support",
                   help="aggregate counts, as from "
                        "'cut -f1 bootstrap_*_splits_output.tsv | grep ^node_ "
                        "| sort | uniq -c'")
    p.add_argument("--clade_samples",
                   help="'node<TAB>sample' per line, from 'matUtils extract "
                        "-I <node> -u' over the accepted splits. Establishes "
                        "the true nesting; without it the tested-against "
                        "sizes are upper bounds and margins are understated "
                        "for late-found splits.")
    p.add_argument("--bootstrap_replicates",
                   help="per-replicate membership, 'replicate<TAB>node' per "
                        "line; needed for neighbourhood support")
    p.add_argument("--neighbourhood", type=int, default=10,
                   help="node-number window for wobble-tolerant support. "
                        "Default 10; 0 disables.")
    p.add_argument("--null_window", type=float, default=1.5,
                   help="fold-range of subtree size accepted as size-matched. "
                        "Default 1.5.")
    p.add_argument("--strong_margin", type=float, default=5.0,
                   help="margin for tier A. Default 5.")
    p.add_argument("--strong_support", type=float, default=70.0,
                   help="support percent for tier A. Default 70.")
    p.add_argument("--output", help="output TSV (default stdout)")
    p.add_argument("--selftest", action="store_true",
                   help="run self-tests and exit")
    return p.parse_args(argv)


def build_rows(args, log_lines, spectra_lines, null_lines,
               support_lines=None, replicate_lines=None, sample_lines=None):
    chi, parent, order = parse_splits_log(log_lines)
    residual = parse_spectra(spectra_lines)
    nulls = parse_null(null_lines)
    sup = parse_support(support_lines) if support_lines else {}
    reps = parse_replicates(replicate_lines) if replicate_lines else {}
    nreps = len(reps)

    sets = read_clade_samples(sample_lines) if sample_lines else None
    if sets:
        inside = nesting_from_samples(sets)
        ### node_1 is the whole tree; it is a search root but never a split.
        inside.setdefault("node_1", set(sets))
        true_parent = {}
        for n in sets:
            sup_n = [m for m in sets if m != n and sets[n] < sets[m]]
            true_parent[n] = (min(sup_n, key=lambda m: len(sets[m]))
                              if sup_n else "node_1")
    else:
        inside = nesting_from_parents(parent)
        inside.setdefault("node_1", set(parent))
        true_parent = {}

    sizes = {order[k][0]: searched_size(order, k, inside, residual)
             for k in range(len(order))}

    rows = []
    for node in sorted(chi, key=lambda n: -chi[n]):
        size = sizes.get(node)
        ceil, nnull = (None, 0) if size is None else ceiling(
            nulls, size, args.null_window)
        margin = None if (ceil is None or ceil == 0) else chi[node] / ceil
        exact = sup.get(node)
        if reps:
            exact = 100.0 * sum(1 for s in reps.values() if node in s) / nreps
        elif exact is not None:
            exact = float(exact)
        nb = None
        if reps and args.neighbourhood > 0:
            nb = 100.0 * neighbourhood_support(
                reps, node, args.neighbourhood) / nreps
        best = nb if nb is not None else exact
        rows.append({
            "node": node, "search_root": parent[node],
            "clade_parent": true_parent.get(node, "NA"),
            "chi2": chi[node], "tested_in_mutations": size,
            "null_ceiling": ceil, "null_replicates": nnull, "margin": margin,
            "support_exact": exact, "support_neighbourhood": nb,
            "tier": tier(margin, best, args.strong_margin,
                         args.strong_support),
        })
    return rows, nreps


COLUMNS = ["node", "search_root", "clade_parent", "chi2",
           "tested_in_mutations", "null_ceiling", "null_replicates",
           "margin", "support_exact", "support_neighbourhood", "tier"]


def write_rows(rows, handle):
    print("\t".join(COLUMNS), file=handle)
    for r in rows:
        out = []
        for c in COLUMNS:
            v = r[c]
            if v is None:
                out.append("NA")
            elif isinstance(v, float):
                out.append(f"{v:.4g}")
            else:
                out.append(str(v))
        print("\t".join(out), file=handle)


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    missing = [n for n in ("splits_log", "spectra", "null")
               if not getattr(args, n)]
    if missing:
        print("error: " + ", ".join("--" + m for m in missing) +
              " required (or pass --selftest)", file=sys.stderr)
        return 2
    if args.bootstrap_support and args.bootstrap_replicates:
        print("note: both bootstrap inputs given; the per-replicate file is "
              "used, since the aggregate cannot give neighbourhood support",
              file=sys.stderr)

    def read(path):
        if not path:
            return None
        with open(path) as fh:
            return fh.readlines()

    rows, nreps = build_rows(
        args, read(args.splits_log), read(args.spectra), read(args.null),
        read(args.bootstrap_support), read(args.bootstrap_replicates),
        read(args.clade_samples))
    if not args.clade_samples:
        print("  WARNING: no --clade_samples, so nesting comes from the log's "
              "search roots alone", file=sys.stderr)
        print("  A split accepted late can contain one accepted earlier, "
              "which that cannot see; sizes are then upper bounds and the "
              "margins of late splits are understated.", file=sys.stderr)

    if not rows:
        print("error: no splits found before the bootstrap marker in "
              f"{args.splits_log}", file=sys.stderr)
        return 1
    nonull = [r["node"] for r in rows if r["margin"] is None]
    if nonull:
        print(f"  WARNING: no size-matched null for {len(nonull)} split(s): "
              f"{', '.join(nonull[:8])}{' ...' if len(nonull) > 8 else ''}",
              file=sys.stderr)
        print("  Widen --null_window, or run permutation_null.py deeper so "
              "the cascade reaches these subtree sizes.", file=sys.stderr)
    if nreps:
        print(f"  {nreps} bootstrap replicates; support is a percent of them",
              file=sys.stderr)
        if args.neighbourhood > 0:
            print(f"  neighbourhood support allows the split to move "
                  f"+-{args.neighbourhood} in node number; ID adjacency is a "
                  f"proxy for topology, so read it as an upper estimate",
                  file=sys.stderr)
    disagree = [r for r in rows if r["clade_parent"] not in ("NA", None)
                and r["clade_parent"] != r["search_root"]]
    if disagree:
        print(f"  note: {len(disagree)} split(s) sit inside a clade other "
              f"than the subtree they were found in, which is expected and is "
              f"why membership decides the nesting", file=sys.stderr)
    counts = collections.Counter(r["tier"] for r in rows)
    print(f"  {len(rows)} splits: {counts['A']} tier A, {counts['B']} tier B, "
          f"{counts['C']} tier C", file=sys.stderr)
    weakest = min((r for r in rows if r["margin"] is not None),
                  key=lambda r: r["margin"], default=None)
    if weakest:
        print(f"  weakest margin {weakest['margin']:.2f}x at "
              f"{weakest['node']}", file=sys.stderr)

    if args.output:
        with open(args.output, "w") as fh:
            write_rows(rows, fh)
        print(f"Wrote {args.output}", file=sys.stderr)
    else:
        write_rows(rows, sys.stdout)
    return 0


def selftest():
    failures = []

    def check(ok, label):
        if ok:
            print(f"ok: {label}")
        else:
            print(f"FAIL: {label}")
            failures.append(label)

    ### 1. Parent assignment, and the bootstrap boundary.  Replicate lines
    ### look identical to real ones, so a parser that reads the whole log
    ### invents splits that were never accepted on the real tree.
    log = [
        "Computing distances between splits in (sub)tree node_1\n",
        "New split found at node_5 with x2 1000.0\n",
        "Computing distances between splits in (sub)tree node_5\n",
        "New split found at node_9 with x2 50.0\n",
        "Bootstrapping splits with 100 replicates using 8 threads.\n",
        "Computing distances between splits in (sub)tree node_1\n",
        "New split found at node_77 with x2 999.0\n",
    ]
    chi, parent, order = parse_splits_log(log)
    check(sorted(chi) == ["node_5", "node_9"]
          and parent == {"node_5": "node_1", "node_9": "node_5"}
          and chi["node_5"] == 1000.0
          and order == [("node_5", "node_1"), ("node_9", "node_5")],
          "splits log gives search roots, order, and stops at the marker")

    ### 2. A concatenated log would double every cumulative total.
    try:
        parse_splits_log(log[:2] + log[:2])
        check(False, "a repeated split raises instead of double-counting")
    except ValueError:
        check(True, "a repeated split raises instead of double-counting")

    ### 3. Residuals do not nest; the subtree total is their sum.  Here
    ### node_1 keeps 10 of its own, node_5 keeps 20 and node_9 keeps 30, so
    ### node_5's subtree is 50 and the whole tree is 60.
    spectra = ["Node_ID\tTotal_Mutations\tNumber_Tips\tMutations:Tips\n",
               "node_1\t10\t5\t2.0\n", "node_5\t20\t5\t4.0\n",
               "node_9\t30\t5\t6.0\n", "\n"]
    residual = parse_spectra(spectra)
    inside = nesting_from_parents(parent)
    inside.setdefault("node_1", set(parent))
    check(subtree_mutations("node_9", inside, residual) == 30
          and subtree_mutations("node_5", inside, residual) == 50
          and subtree_mutations("node_1", inside, residual) == 60,
          "subtree totals add residuals back up")

    ### Membership overrides the log, because a split accepted late can be an
    ### ancestor of one accepted earlier.  Here node_7 was found second yet
    ### contains node_9, and both name node_5 as their search root.
    sets = {"node_5": {"a", "b", "c"}, "node_7": {"b", "c"}, "node_9": {"c"}}
    ins = nesting_from_samples(sets)
    check(ins["node_5"] == {"node_7", "node_9"} and ins["node_7"] == {"node_9"}
          and ins["node_9"] == set(),
          "membership finds a late split that contains an earlier one")

    ### The carving.  Three splits found in node_1 in this order; by the time
    ### the third is tested the first two are gone.  node_9 is nested inside
    ### node_5, so subtracting it again would undercount the table -- the bug
    ### the maximal-prior rule exists to stop.
    ins2 = {"node_1": {"node_5", "node_9", "node_7"},
            "node_5": {"node_9"}, "node_9": set(), "node_7": set()}
    res2 = {"node_1": 100, "node_5": 20, "node_9": 30, "node_7": 40}
    ord2 = [("node_5", "node_1"), ("node_9", "node_5"), ("node_7", "node_1")]
    sizes = [searched_size(ord2, k, ins2, res2) for k in range(3)]
    ### whole tree 190; node_5's subtree is 50, so the third sees 140
    check(sizes == [190, 50, 140],
          "tested-against size shrinks as earlier splits are carved out")

    ### Without the maximal rule this would read 190-50-30 = 110.
    check(searched_size([("node_5", "node_1"), ("node_9", "node_5"),
                         ("node_7", "node_1")], 2, ins2, res2) == 140,
          "a prior split nested in another prior split is not subtracted twice")

    ### 4. Comparing an observation against itself would make every margin 1.
    nulls = parse_null([
        "replicate\titeration\tsplit_root_id\tnode_id\tchi2\tmutations_below"
        "\tmutations_above\n",
        "observed\t1\tnode_1\tnode_5\t9999.0\t10\t50\n",
        "1\t1\tnode_1\tnode_5\t12.0\t10\t50\n",
        "2\t1\tnode_1\tnode_5\t15.0\t20\t40\n",
    ])
    check(nulls == [(60, 12.0), (60, 15.0)],
          "null parsing drops the observed rows")

    ### 5. The window is multiplicative in both directions, and an unmatched
    ### size reports nothing rather than borrowing a distant ceiling.
    c1, n1 = ceiling([(100, 5.0), (140, 9.0), (1000, 99.0)], 100, 1.5)
    c2, n2 = ceiling([(1000, 99.0)], 100, 1.5)
    check(c1 == 9.0 and n1 == 2 and c2 is None and n2 == 0,
          "size-matched ceiling takes the window maximum, or nothing")

    ### 6. The failure this script exists to avoid.  Replicate 1 splits three
    ### nodes inside the window; summing per-node counts would score it 3, and
    ### over 100 replicates such sums exceed 100% and are not support at all.
    reps = parse_replicates([
        "1\tnode_50\n", "1\tnode_51\n", "1\tnode_52\n", "2\tnode_51\n",
        "3\tnode_900\n",
    ])
    nb = neighbourhood_support(reps, "node_50", 5)
    exact = sum(1 for s in reps.values() if "node_50" in s)
    check(len(reps) == 3 and nb == 2 and nb <= len(reps) and exact == 1,
          "neighbourhood support counts replicates once, never summing nodes")

    ### 7. An unknown margin or support must not be read as passing.
    check(tier(10.0, 90.0, 5, 70) == "A" and tier(10.0, 10.0, 5, 70) == "B"
          and tier(1.0, 90.0, 5, 70) == "B" and tier(1.0, 10.0, 5, 70) == "C"
          and tier(None, 90.0, 5, 70) == "B" and tier(None, None, 5, 70) == "C"
          and tier(10.0, None, 5, 70) == "B",
          "tiers need both criteria for A and treat unknowns as failing")

    ### 8. End to end, including a split whose subtree size the null misses.
    args = parse_args(["--selftest"])
    nulltsv = [
        "replicate\titeration\tsplit_root_id\tnode_id\tchi2\tmutations_below"
        "\tmutations_above\n",
        "1\t1\tnode_1\tnode_5\t20.0\t30\t30\n",
        "1\t2\tnode_5\tnode_9\t10.0\t25\t25\n",
    ]
    reptsv = ["1\tnode_5\n", "2\tnode_5\n", "2\tnode_9\n", "3\tnode_6\n"]
    rows, nreps = build_rows(args, log, spectra, nulltsv, None, reptsv, None)
    by = {r["node"]: r for r in rows}
    check(nreps == 3 and len(rows) == 2
          and abs(by["node_5"]["margin"] - 1000.0 / 20.0) < 1e-9
          and abs(by["node_5"]["support_exact"] - 200.0 / 3) < 1e-6
          and abs(by["node_5"]["support_neighbourhood"] - 100.0) < 1e-6
          and by["node_5"]["tier"] == "A",
          "end to end: margin, support and tier for a strong split")
    ### node_9 was tested inside node_5, so its size is node_5's subtree (50)
    ### and not the root's 60.  At the default 1.5-fold window both null rows
    ### are in range, so the ceiling is their maximum.
    check(by["node_9"]["tested_in_mutations"] == 50
          and by["node_9"]["null_ceiling"] == 20.0
          and abs(by["node_9"]["margin"] - 2.5) < 1e-9,
          "a deep split is sized by its own subtree, not by the root")

    ### Narrowing the window drops the 60-mutation null and leaves only the
    ### one drawn at node_9's own scale, which is the whole point of matching.
    narrow = parse_args(["--null_window", "1.05"])
    rows_n, _ = build_rows(narrow, log, spectra, nulltsv, None, reptsv, None)
    deep = {r["node"]: r for r in rows_n}["node_9"]
    check(deep["null_ceiling"] == 10.0 and deep["null_replicates"] == 1
          and abs(deep["margin"] - 5.0) < 1e-9,
          "a narrower window keeps only nulls at the split's own scale")

    print("FAILED" if failures else "all self-tests passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
