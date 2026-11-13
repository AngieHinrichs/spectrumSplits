#!/usr/bin/env python3
"""Calibrate spectrumSplits.py's --min_chi against a permutation null.

This is step 4 of parameter tuning.  spectrumSplits.py scores every node in a
subtree and keeps the single best chi-squared value, so the largest chi-squared
it finds is a maximum over hundreds of thousands of correlated tests.  A
nominal significance cutoff is therefore meaningless: --min_chi has to clear
whatever the search turns up when there is no clade-specific spectrum at all.

This script measures that.  It destroys any association between mutation type
and position in the tree, then runs the same search and the same statistic
spectrumSplits.py uses, and reports the chi-squared values it still finds.
Set --min_chi comfortably above the maximum observed here.

Null models (--mode):

  types  Shuffle mutation-type labels across all mutations in the tree, keeping
         the topology and the number of mutations per branch fixed.  This is
         the strict independence null: it removes clade-specific spectra and
         also removes the type clumping caused by recurrent sites, so it tends
         to give the more permissive (lower) threshold.

  sites  Give every site a single mutation type (its most common one), then
         shuffle that site-to-type assignment across sites.  Which branches
         carry mutations at which site is untouched, so a site that mutates
         500 times inside one clade still dumps 500 counts of one type into
         that clade - the type is just randomized.  This preserves the
         independence violation that recurrent sites cause, and so gives the
         more conservative (higher) threshold.  This is usually the one to
         calibrate against, especially for a tree with convergent selection.

The statistic matches spectrumSplits.py: a 2 x 12 contingency table of the
candidate subtree's spectrum against the rest of the current split, scored
with Pearson's chi-squared (no continuity correction, as scipy applies none
for tables larger than 2 x 2).  Columns whose total is zero are skipped, where
scipy would raise; this only arises when a mutation type is entirely absent.
Pass --fold to score the strand-folded 6-type spectrum instead.

Usage:
    # headline calibration: one pass over the whole tree, per replicate
    permutation_null.py --input_tree tree.pb.gz --min_mutations 500 \
        --replicates 20 --mode sites --include_observed

    # explore the full greedy cascade rather than just the first split
    permutation_null.py --input_tree tree.pb.gz --min_mutations 500 \
        --replicates 5 --max_splits 25 --min_chi 1
"""

import argparse
import math
import random
import sys
from array import array
from collections import Counter, defaultdict

# Must match the ordering used in spectrumSplits.py find_splits().
MUT_TYPES = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
MUT_INDEX = {m: i for i, m in enumerate(MUT_TYPES)}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def folded_label(mut_type):
    """Label for the strand-complementary pair a mutation class belongs to."""
    ref, alt = mut_type[0], mut_type[1]
    reverse = COMPLEMENT[ref] + COMPLEMENT[alt]
    if ref in ("C", "T"):
        return mut_type + "|" + reverse
    return reverse + "|" + mut_type


def folded_index():
    """Map each of the 12 type indices onto one of 6 strand-folded indices.

    Labels are sorted so the category order matches profile_tree.py and
    spectrumSplits.py; the statistic is order-independent but the reported
    order should not differ between the tools.
    """
    labels = sorted({folded_label(m) for m in MUT_TYPES})
    mapping = [labels.index(folded_label(m)) for m in MUT_TYPES]
    return mapping, labels


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calibrate spectrumSplits.py --min_chi against a permutation null."
    )
    parser.add_argument("--input_tree", type=str, required=True,
                        help="Input tree file (protobuf format)")
    parser.add_argument("--output", type=str, default="permutation_null.tsv",
                        help="Output TSV of every split scored, per replicate")
    parser.add_argument("--replicates", type=int, default=20,
                        help="Number of permutation replicates")
    parser.add_argument("--mode", type=str, default="sites", choices=["sites", "types"],
                        help="Null model (see module docstring)")
    parser.add_argument("--min_mutations", type=int, default=500,
                        help="Match the --min_mutations you intend to use for the real run")
    parser.add_argument("--max_branch_length", type=int, default=100000,
                        help="Match the --max_branch_length you intend to use")
    parser.add_argument("--min_chi", type=float, default=1.0,
                        help="Permissive threshold used during the null search; only "
                             "matters when --max_splits > 1")
    parser.add_argument("--max_splits", type=int, default=1,
                        help="Splits to accept per replicate. 1 (default) scores a single "
                             "pass over the whole tree, which is the quantity the first "
                             "real split must beat. Higher values follow the greedy "
                             "cascade into subtrees, where chi-squared is necessarily smaller.")
    parser.add_argument("--fold", action="store_true",
                        help="Score the strand-folded 6-type spectrum instead of 12 types")
    parser.add_argument("--include_observed", action="store_true",
                        help="Also score the unpermuted tree through the same code path")
    parser.add_argument("--seed", type=int, default=1,
                        help="Base random seed; replicate r uses seed + r")
    parser.add_argument("--nthreads", type=int, default=1,
                        help="Run replicates concurrently (requires fork; falls back to serial)")
    parser.add_argument("--selftest", action="store_true",
                        help="Check the chi-squared implementation and exit")
    return parser.parse_args()


# --------------------------------------------------------------------------
# chi-squared
# --------------------------------------------------------------------------

def chi2_2xk(below, above):
    """Pearson chi-squared for a 2 x k contingency table.

    Matches scipy.stats.chi2_contingency([below, above]) for k > 2, which
    applies no continuity correction.  Columns with a zero total are skipped.
    Returns 0.0 if either row total is zero.
    """
    rb = 0
    ra = 0
    for v in below:
        rb += v
    for v in above:
        ra += v
    n = rb + ra
    if rb <= 0 or ra <= 0 or n <= 0:
        return 0.0
    chi2 = 0.0
    for j in range(len(below)):
        col = below[j] + above[j]
        if col == 0:
            continue
        e_b = rb * col / n
        e_a = ra * col / n
        d_b = below[j] - e_b
        d_a = above[j] - e_a
        chi2 += d_b * d_b / e_b + d_a * d_a / e_a
    return chi2


def selftest():
    """Verify chi2_2xk against a naive implementation and against scipy."""
    rng = random.Random(12345)
    worst = 0.0
    for trial in range(200):
        k = rng.choice([6, 12])
        below = [rng.randint(0, 500) for _ in range(k)]
        above = [rng.randint(0, 5000) for _ in range(k)]
        mine = chi2_2xk(below, above)

        # naive reference
        rows = [below, above]
        row_tot = [sum(r) for r in rows]
        col_tot = [below[j] + above[j] for j in range(k)]
        n = sum(row_tot)
        ref = 0.0
        for i in range(2):
            for j in range(k):
                if col_tot[j] == 0:
                    continue
                e = row_tot[i] * col_tot[j] / n
                ref += (rows[i][j] - e) ** 2 / e
        if abs(mine - ref) > 1e-9 * max(1.0, ref):
            print(f"FAIL naive: trial {trial} mine={mine} ref={ref}", file=sys.stderr)
            return 1
        worst = max(worst, abs(mine - ref))
    print(f"ok: matches naive reference on 200 random tables (max abs diff {worst:.3e})")

    try:
        from scipy.stats import chi2_contingency
    except ImportError:
        print("note: scipy not importable here; skipped the scipy comparison. "
              "Re-run --selftest in your spectrumSplits environment to confirm.")
        return 0
    worst = 0.0
    for trial in range(200):
        k = rng.choice([6, 12])
        below = [rng.randint(1, 500) for _ in range(k)]
        above = [rng.randint(1, 5000) for _ in range(k)]
        mine = chi2_2xk(below, above)
        ref = chi2_contingency([below, above])[0]
        if abs(mine - ref) > 1e-8 * max(1.0, ref):
            print(f"FAIL scipy: trial {trial} mine={mine} ref={ref}", file=sys.stderr)
            return 1
        worst = max(worst, abs(mine - ref))
    print(f"ok: matches scipy.stats.chi2_contingency on 200 random tables "
          f"(max abs diff {worst:.3e})")
    return 0


# --------------------------------------------------------------------------
# tree extraction
# --------------------------------------------------------------------------

class FlatTree(object):
    """Topology plus per-branch mutations, independent of bte.

    Parent indices are strictly less than child indices, so iterating indices
    in reverse visits children before parents.
    """

    def __init__(self):
        self.node_ids = []
        self.parent = array("i")
        self.children = []
        self.raw_len = array("i")      # len(node.mutations), as spectrumSplits sees it
        self.positions = []            # array('i') per node, parsed SNVs only
        self.types = []                # array('b') per node, aligned with positions


def parse_mutation(mutation):
    if len(mutation) < 3:
        return None
    ref, alt = mutation[0].upper(), mutation[-1].upper()
    if ref not in COMPLEMENT or alt not in COMPLEMENT or ref == alt:
        return None
    try:
        pos = int(mutation[1:-1])
    except ValueError:
        return None
    return pos, MUT_INDEX[ref + alt]


def extract(root):
    tree = FlatTree()
    stack = [(root, -1)]
    while stack:
        node, p = stack.pop()
        i = len(tree.node_ids)
        tree.node_ids.append(node.id)
        tree.parent.append(p)
        tree.children.append([])
        if p >= 0:
            tree.children[p].append(i)
        muts = node.mutations
        tree.raw_len.append(len(muts))
        pos_arr = array("i")
        type_arr = array("b")
        for mutation in muts:
            parsed = parse_mutation(mutation)
            if parsed is None:
                continue
            pos_arr.append(parsed[0])
            type_arr.append(parsed[1])
        tree.positions.append(pos_arr)
        tree.types.append(type_arr)
        for child in node.children:
            stack.append((child, i))
    return tree


# --------------------------------------------------------------------------
# permutation
# --------------------------------------------------------------------------

def permute_types(tree, rng):
    """Shuffle type labels across all mutations in the tree."""
    flat = array("b")
    for t in tree.types:
        flat.extend(t)
    shuffled = list(flat)
    rng.shuffle(shuffled)
    out = []
    k = 0
    for t in tree.types:
        n = len(t)
        out.append(array("b", shuffled[k:k + n]))
        k += n
    return out


def site_modal_types(tree):
    """Most common mutation type at each site."""
    per_site = defaultdict(Counter)
    for pos_arr, type_arr in zip(tree.positions, tree.types):
        for pos, t in zip(pos_arr, type_arr):
            per_site[pos][t] += 1
    return {pos: counter.most_common(1)[0][0] for pos, counter in per_site.items()}


def permute_sites(tree, rng, modal):
    """Shuffle the site-to-type assignment, keeping site positions in place."""
    sites = list(modal.keys())
    assigned = [modal[p] for p in sites]
    rng.shuffle(assigned)
    new_map = dict(zip(sites, assigned))
    out = []
    for pos_arr in tree.positions:
        out.append(array("b", [new_map[p] for p in pos_arr]))
    return out


# --------------------------------------------------------------------------
# the search, mirroring spectrumSplits.py find_splits()
# --------------------------------------------------------------------------

def own_spectra(tree, types, ncat, fold_map, max_branch_length):
    """Per-branch spectrum vectors, honouring --max_branch_length."""
    n = len(tree.node_ids)
    out = [None] * n
    zero = [0] * ncat
    for i in range(n):
        t = types[i]
        if not t or tree.raw_len[i] > max_branch_length:
            out[i] = zero
            continue
        vec = [0] * ncat
        if fold_map is None:
            for x in t:
                vec[x] += 1
        else:
            for x in t:
                vec[fold_map[x]] += 1
        out[i] = vec
    return out


def scan_subtree(tree, own, ncat, split_root, stop_set, min_mutations):
    """Score every candidate node in one subtree; return (max_chi2, node, below, above).

    Mirrors the inner loop of find_splits(): descend from split_root, skipping
    any child that is an accepted split, accumulate subtree spectra, and score
    each node's subtree against the remainder of the split root's spectrum.
    """
    # in-scope nodes, parents before children
    scope = []
    stack = [split_root]
    while stack:
        v = stack.pop()
        scope.append(v)
        for c in tree.children[v]:
            if c in stop_set:
                continue
            stack.append(c)

    local = {v: k for k, v in enumerate(scope)}
    spec = [list(own[v]) for v in scope]
    for k in range(len(scope) - 1, 0, -1):
        tgt = spec[local[tree.parent[scope[k]]]]
        src = spec[k]
        for j in range(ncat):
            tgt[j] += src[j]

    root_spec = spec[0]
    root_total = 0
    for v in root_spec:
        root_total += v

    max_chi = 0.0
    max_node = None
    max_below = 0
    for k in range(len(scope)):
        below = spec[k]
        n_below = 0
        for v in below:
            n_below += v
        if n_below < min_mutations:
            continue
        if root_total - n_below < min_mutations:
            continue
        above = [root_spec[j] - below[j] for j in range(ncat)]
        chi = chi2_2xk(below, above)
        if chi > max_chi:
            max_chi = chi
            max_node = scope[k]
            max_below = n_below
    return max_chi, max_node, max_below, root_total - max_below


def run_search(tree, own, ncat, min_chi, min_mutations, max_splits, log, label):
    """Greedy split search; returns the list of accepted splits as records."""
    accepted = {0}
    finalized = set()
    records = []
    iteration = 0
    while len(accepted) > len(finalized) and len(records) < max_splits:
        iteration += 1
        new_splits = set()
        for split_root in sorted(accepted):
            if split_root in finalized:
                continue
            chi, node, n_below, n_above = scan_subtree(
                tree, own, ncat, split_root, accepted, min_mutations)
            if node is not None and chi > min_chi and node not in accepted:
                new_splits.add(node)
                records.append({
                    "label": label,
                    "iteration": iteration,
                    "split_root_id": tree.node_ids[split_root],
                    "node_id": tree.node_ids[node],
                    "chi2": chi,
                    "mutations_below": n_below,
                    "mutations_above": n_above,
                })
                print(f"  [{label}] iter {iteration}: split at {tree.node_ids[node]} "
                      f"chi2={chi:.1f} below={n_below} above={n_above}", file=log)
                if len(records) >= max_splits:
                    break
            else:
                finalized.add(split_root)
                if node is None:
                    print(f"  [{label}] iter {iteration}: no eligible candidate under "
                          f"{tree.node_ids[split_root]} at --min_mutations "
                          f"{min_mutations}", file=log)
        accepted |= new_splits
        if not new_splits:
            break
    return records


# --------------------------------------------------------------------------
# replicate driver
# --------------------------------------------------------------------------

def one_replicate(tree, args, ncat, fold_map, modal, replicate):
    """Score one replicate; replicate 0 means the unpermuted tree."""
    if replicate == 0:
        types = tree.types
        label = "observed"
    else:
        rng = random.Random(args.seed + replicate)
        if args.mode == "types":
            types = permute_types(tree, rng)
        else:
            types = permute_sites(tree, rng, modal)
        label = f"perm{replicate}"
    own = own_spectra(tree, types, ncat, fold_map, args.max_branch_length)
    print(f"Scoring {label} ...", file=sys.stderr)
    return run_search(tree, own, ncat, args.min_chi, args.min_mutations,
                      args.max_splits, sys.stderr, label)


def quantile(sorted_values, p):
    if not sorted_values:
        return float("nan")
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    h = (n - 1) * p
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    return sorted_values[lo] + (h - lo) * (sorted_values[hi] - sorted_values[lo])


def main():
    args = parse_args()

    if args.selftest:
        return selftest()

    fold_map, fold_labels = folded_index()
    if args.fold:
        ncat = 6
        cats = fold_labels
    else:
        ncat = 12
        fold_map = None
        cats = MUT_TYPES

    try:
        import bte
    except ImportError:
        print("error: could not import bte; run this in the environment you use "
              "for spectrumSplits.py", file=sys.stderr)
        return 1

    print(f"Reading {args.input_tree} ...", file=sys.stderr)
    mat = bte.MATree(args.input_tree)
    tree = extract(mat.root)
    del mat
    n_nodes = len(tree.node_ids)
    n_muts = sum(len(t) for t in tree.types)
    print(f"Extracted {n_nodes} nodes, {n_muts} usable SNVs, "
          f"scoring {ncat} categories ({', '.join(cats)})", file=sys.stderr)

    modal = site_modal_types(tree) if args.mode == "sites" else None
    if modal is not None:
        print(f"Site-level null over {len(modal)} distinct sites", file=sys.stderr)

    replicates = list(range(1, args.replicates + 1))
    if args.include_observed:
        replicates = [0] + replicates

    results = []
    use_fork = args.nthreads > 1
    if use_fork:
        try:
            import multiprocessing
            ctx = multiprocessing.get_context("fork")
        except (ImportError, ValueError):
            print("warning: fork start method unavailable; running serially. "
                  "Run several instances with different --seed to parallelize.",
                  file=sys.stderr)
            use_fork = False

    if use_fork:
        with ctx.Pool(processes=args.nthreads) as pool:
            jobs = [pool.apply_async(one_replicate,
                                     (tree, args, ncat, fold_map, modal, r))
                    for r in replicates]
            for job in jobs:
                results.extend(job.get())
    else:
        for r in replicates:
            results.extend(one_replicate(tree, args, ncat, fold_map, modal, r))

    with open(args.output, "w") as fh:
        fh.write("\t".join(["replicate", "iteration", "split_root_id", "node_id",
                            "chi2", "mutations_below", "mutations_above"]) + "\n")
        for rec in results:
            fh.write("\t".join([
                rec["label"], str(rec["iteration"]), rec["split_root_id"],
                rec["node_id"], f"{rec['chi2']:.4f}",
                str(rec["mutations_below"]), str(rec["mutations_above"]),
            ]) + "\n")
    print(f"\nPer-split records written to {args.output}", file=sys.stderr)

    # ---- summary ---------------------------------------------------------
    per_label = defaultdict(list)
    for rec in results:
        per_label[rec["label"]].append(rec["chi2"])

    null_maxima = sorted(max(v) for k, v in per_label.items()
                         if k != "observed" and v)
    out = sys.stdout
    print("=" * 72, file=out)
    print("PERMUTATION NULL SUMMARY", file=out)
    print("=" * 72, file=out)
    print(f"  mode:               {args.mode}", file=out)
    print(f"  categories:         {ncat}", file=out)
    print(f"  min_mutations:      {args.min_mutations}", file=out)
    print(f"  max_branch_length:  {args.max_branch_length}", file=out)
    print(f"  max_splits:         {args.max_splits}", file=out)
    print(f"  replicates scored:  {len(null_maxima)}", file=out)
    print(file=out)

    if not null_maxima:
        print("  No null replicate produced a scoreable split.  Either "
              "--min_mutations", file=out)
        print("  excludes every candidate, or --min_chi is above the null range.", file=out)
    else:
        print("  Maximum chi-squared per null replicate:", file=out)
        print(f"    min    {null_maxima[0]:.1f}", file=out)
        print(f"    median {quantile(null_maxima, 0.5):.1f}", file=out)
        print(f"    p95    {quantile(null_maxima, 0.95):.1f}", file=out)
        print(f"    max    {null_maxima[-1]:.1f}", file=out)
        print(file=out)
        floor = null_maxima[-1]
        print(f"  The null search reached chi-squared {floor:.1f} with no real", file=out)
        print("  clade-specific spectrum present.  --min_chi must exceed that:", file=out)
        print(f"    bare floor         --min_chi {math.ceil(floor)}", file=out)
        print(f"    with 1.5x margin   --min_chi {math.ceil(1.5 * floor)}", file=out)
        print(f"    with 2x margin     --min_chi {math.ceil(2.0 * floor)}", file=out)
        print(file=out)
        print("  Only the replicate maxima matter here; a replicate's later", file=out)
        print("  splits are scored on smaller subtrees and are not comparable.", file=out)

    if "observed" in per_label and per_label["observed"]:
        obs = max(per_label["observed"])
        print(file=out)
        print(f"  Observed (unpermuted) maximum chi-squared: {obs:.1f}", file=out)
        if null_maxima:
            ratio = obs / null_maxima[-1] if null_maxima[-1] > 0 else float("inf")
            print(f"  Ratio to the worst null replicate: {ratio:.1f}x", file=out)
            if ratio < 2:
                print("  That is a thin margin: the strongest real split is close to", file=out)
                print("  what chance alone produces.  Treat any split here as", file=out)
                print("  provisional and check it against the bootstrap.", file=out)

    print(file=out)
    print("  Next: run spectrumSplits.py with this --min_chi, then confirm the", file=out)
    print("  surviving splits with --bootstrap_splits, and inspect which", file=out)
    print("  mutation classes drive each one.", file=out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
