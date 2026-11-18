#!/usr/bin/env python3
"""Compare internal-branch and terminal-branch spectra within each split.

This is step 9 of parameter tuning: a check on whether a spectrum split is
biology or sequencing quality.

The problem it addresses.  On a sparsely sampled tree most mutations sit on
terminal branches -- 63% on the MTBC tree this was written for -- and a
terminal branch carries everything wrong with one genome: mapping error in
repeats, deamination damage, contamination, assembly artifacts.  Those have
their own spectrum, usually enriched for C>T and G>A.  So a clade whose
genomes are individually worse will show a shifted spectrum for reasons that
have nothing to do with its mutational processes, and spectrumSplits.py will
happily split on it.

Internal branches cannot carry per-genome error: a mutation is placed on an
internal branch only because it is shared by the samples below it, so an
artifact of one genome cannot reach there.  (Shared artifacts can -- a repeat
that misaligns the same way in every genome of a lineage -- which is why
masking comes first, not instead.)

The test.  If a split's spectrum difference is real, it holds when only
internal branches are counted.  If it is a quality gradient, it collapses.

What to look at, in order:

  1. Does the BETWEEN-SPLIT ordering survive on internal branches alone?
     This is the test.  Correlations of each class against the divergence
     ratio are printed for both branch classes; if the lineage structure is
     real, the internal-only numbers keep it.

  2. How large is the WITHIN-SPLIT internal-vs-terminal difference?  A large
     chi-squared here is expected on any real tree -- terminal branches
     genuinely differ -- and is not itself evidence of a problem.  It matters
     when it varies wildly between splits, which says the quality gradient is
     uneven across clades.

Usage:
    branch_class_spectra.py --input_tree tree.pb.gz \
        --spectra spectra.tsv --fold --max_branch_length 250
"""

import argparse
import os
import sys
from array import array
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
### Imported rather than reimplemented so the folding, the parse rules and the
### statistic cannot drift from the tool this is meant to check.
from permutation_null import (MUT_TYPES, chi2_2xk, extract, folded_index,
                              parse_mutation)


def parse_args():
    p = argparse.ArgumentParser(
        description="Compare internal-branch and terminal-branch spectra "
                    "within each split.")
    p.add_argument("--input_tree", type=str, default=None,
                   help="Input tree file (protobuf format)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--spectra", type=str, default=None,
                   help="spectrumSplits.py --output_spectrum TSV; split roots "
                        "are read from its Node_ID column")
    g.add_argument("--splits", type=str, default=None,
                   help="File of split root node IDs, one per line")
    p.add_argument("--output", type=str, default=None,
                   help="Output TSV of per-split, per-branch-class spectra")
    p.add_argument("--fold", action="store_true",
                   help="Score the strand-folded 6-type spectrum instead of 12")
    p.add_argument("--max_branch_length", type=int, default=None,
                   help="Match the --max_branch_length used for the search; "
                        "branches above it contribute nothing")
    p.add_argument("--selftest", action="store_true",
                   help="Check group assignment and partitioning, then exit")
    return p.parse_args()


def read_split_roots(args):
    if args.splits:
        with open(args.splits) as fh:
            return [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    with open(args.spectra) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            col = hdr.index("Node_ID")
        except ValueError:
            raise SystemExit("error: --spectra file has no Node_ID column")
        out = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > col and f[col]:
                out.append(f[col])
        return out


def residual_groups(tree, split_roots):
    """Nearest ancestor-or-self split root for every node.

    This reproduces how spectrumSplits.py reports spectra: once a child split
    is carved out, its parent's row describes only what is left.  Parent
    indices are strictly less than child indices, so one forward pass
    suffices.
    """
    index = {nid: i for i, nid in enumerate(tree.node_ids)}
    missing = [s for s in split_roots if s not in index]
    roots = set(index[s] for s in split_roots if s in index)
    group = array("i", [-1] * len(tree.node_ids))
    for i in range(len(tree.node_ids)):
        if i in roots:
            group[i] = i
        else:
            p = tree.parent[i]
            group[i] = group[p] if p >= 0 else -1
    return group, roots, missing


def tabulate(tree, group, roots, ncat, fold_map, cap):
    """Per group, per branch class, counts of each mutation category.

    Returns counts[(group, is_terminal)] -> list of ncat ints, plus tip counts
    per group and the number of branches skipped by the cap.
    """
    counts = {}
    tips = Counter()
    skipped_branches = 0
    skipped_muts = 0
    for i in range(len(tree.node_ids)):
        g = group[i]
        terminal = not tree.children[i]
        if terminal and g >= 0:
            tips[g] += 1
        if g < 0:
            continue
        if cap is not None and tree.raw_len[i] > cap:
            skipped_branches += 1
            skipped_muts += len(tree.types[i])
            continue
        key = (g, terminal)
        row = counts.get(key)
        if row is None:
            row = [0] * ncat
            counts[key] = row
        for t in tree.types[i]:
            row[fold_map[t] if fold_map else t] += 1
    return counts, tips, skipped_branches, skipped_muts


def props(row):
    n = sum(row)
    return [v / n for v in row] if n else [0.0] * len(row)


def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0] * len(v)
    for pos, i in enumerate(order):
        r[i] = pos + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    sx = sum((a - mx) ** 2 for a in x) ** 0.5
    sy = sum((b - my) ** 2 for b in y) ** 0.5
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def spearman(x, y):
    return pearson(rank(x), rank(y))


def selftest():
    """Synthetic tree: group assignment and the internal/terminal partition."""
    class N(object):
        def __init__(self, i, muts, kids=()):
            self.id = i
            self.mutations = list(muts)
            self.children = list(kids)

    ### root -> a -> (leaf1, leaf2);  root -> b(leaf)
    leaf1 = N("leaf1", ["C100T"])
    leaf2 = N("leaf2", ["C101T"])
    a = N("a", ["A200G", "A201G"], [leaf1, leaf2])
    b = N("b", ["C300A"])
    root = N("root", ["G400A"], [a, b])

    tree = extract(root)
    fold_map, labels = folded_index()

    group, roots, missing = residual_groups(tree, ["root", "a"])
    assert not missing, missing
    index = {nid: i for i, nid in enumerate(tree.node_ids)}
    ### a and its leaves belong to a; root and b belong to root.
    for nid, owner in (("root", "root"), ("b", "root"),
                       ("a", "a"), ("leaf1", "a"), ("leaf2", "a")):
        assert tree.node_ids[group[index[nid]]] == owner, nid
    print("ok: residual group assignment follows the nearest split ancestor")

    counts, tips, skipped, skipped_m = tabulate(
        tree, group, roots, 6, fold_map, None)
    assert skipped == 0 and skipped_m == 0
    ### leaf1 + leaf2 are terminal and in group a: two C>T.
    ct = labels.index("CT|GA")
    assert counts[(index["a"], True)][ct] == 2, counts[(index["a"], True)]
    ### a's own branch is internal: two A>G, which folds onto TC|AG.
    tc = labels.index("TC|AG")
    assert counts[(index["a"], False)][tc] == 2, counts[(index["a"], False)]
    ### b is a leaf in group root: one C>A.  root's branch is internal: G>A.
    ca = labels.index("CA|GT")
    assert counts[(index["root"], True)][ca] == 1
    assert counts[(index["root"], False)][ct] == 1
    assert tips[index["a"]] == 2 and tips[index["root"]] == 1
    print("ok: internal and terminal branches partition without overlap")

    ### The cap drops whole branches, mutations and all.
    counts2, _, sk2, skm2 = tabulate(tree, group, roots, 6, fold_map, 1)
    assert sk2 == 1 and skm2 == 2, (sk2, skm2)
    ### A branch dropped by the cap contributes nothing, so no row is created.
    assert counts2.get((index["a"], False), [0] * 6)[tc] == 0
    print("ok: --max_branch_length drops a whole branch, not part of one")

    ### A tree whose root is not a split root leaves those mutations ungrouped.
    group2, _, _ = residual_groups(tree, ["a"])
    assert group2[index["root"]] == -1 and group2[index["b"]] == -1
    print("ok: mutations above the shallowest split root are reported, not "
          "silently assigned")

    assert abs(chi2_2xk([10, 10], [10, 10])) < 1e-12
    print("ok: imported chi-squared agrees with itself on a null table")
    print("all self-tests passed")
    return 0


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    if not args.input_tree:
        print("error: --input_tree is required (or pass --selftest)",
              file=sys.stderr)
        return 2
    if not (args.spectra or args.splits):
        print("error: pass --spectra or --splits to name the split roots",
              file=sys.stderr)
        return 2

    fold_map, fold_labels = folded_index()
    if args.fold:
        ncat, cats = 6, fold_labels
    else:
        ncat, cats, fold_map = 12, MUT_TYPES, None

    split_roots = read_split_roots(args)
    print(f"{len(split_roots)} split roots requested", file=sys.stderr)

    try:
        import bte
    except ImportError:
        print("error: could not import bte; run this in the environment you "
              "use for spectrumSplits.py", file=sys.stderr)
        return 1
    print(f"Reading {args.input_tree} ...", file=sys.stderr)
    mat = bte.MATree(args.input_tree)
    tree = extract(mat.root)
    del mat

    group, roots, missing = residual_groups(tree, split_roots)
    if missing:
        print(f"  WARNING: {len(missing)} split root(s) not in this tree: "
              f"{', '.join(missing[:10])}"
              f"{' ...' if len(missing) > 10 else ''}", file=sys.stderr)
        print("  Node IDs are assigned at read time, so a list from a "
              "differently masked or pruned tree will not match.",
              file=sys.stderr)
    ungrouped = sum(1 for g in group if g < 0)
    if ungrouped:
        print(f"  note: {ungrouped} nodes sit above every split root; their "
              f"mutations are excluded", file=sys.stderr)

    counts, tips, skipped, skipped_m = tabulate(
        tree, group, roots, ncat, fold_map, args.max_branch_length)
    if args.max_branch_length is not None:
        print(f"  --max_branch_length {args.max_branch_length} dropped "
              f"{skipped} branches and {skipped_m} mutations", file=sys.stderr)

    ordered = sorted(roots, key=lambda g: tree.node_ids[g])
    rows = []
    for g in ordered:
        internal = counts.get((g, False), [0] * ncat)
        terminal = counts.get((g, True), [0] * ncat)
        ni, nt = sum(internal), sum(terminal)
        rows.append({
            "node": tree.node_ids[g],
            "tips": tips[g],
            "internal": internal, "terminal": terminal,
            "ni": ni, "nt": nt,
            "chi2": chi2_2xk(internal, terminal),
            "mt": (ni + nt) / tips[g] if tips[g] else 0.0,
            "term_frac": nt / (ni + nt) if (ni + nt) else 0.0,
        })

    out = sys.stdout
    print("\n" + "=" * 72, file=out)
    print("INTERNAL vs TERMINAL BRANCH SPECTRA, PER SPLIT", file=out)
    print("=" * 72, file=out)
    tot_i = [sum(r["internal"][j] for r in rows) for j in range(ncat)]
    tot_t = [sum(r["terminal"][j] for r in rows) for j in range(ncat)]
    print(f"  mutations on internal branches: {sum(tot_i):,}", file=out)
    print(f"  mutations on terminal branches: {sum(tot_t):,} "
          f"({100 * sum(tot_t) / max(1, sum(tot_i) + sum(tot_t)):.1f}%)",
          file=out)
    print(f"\n  whole-tree spectrum by branch class:", file=out)
    print(f"    {'class':8s} {'internal':>10s} {'terminal':>10s} {'ratio':>7s}",
          file=out)
    pi, pt = props(tot_i), props(tot_t)
    for j, c in enumerate(cats):
        ratio = pt[j] / pi[j] if pi[j] else float("nan")
        print(f"    {c:8s} {pi[j]:10.4f} {pt[j]:10.4f} {ratio:7.2f}", file=out)
    print(f"    internal vs terminal chi-squared: "
          f"{chi2_2xk(tot_i, tot_t):,.1f} on {ncat - 1} df", file=out)

    print(f"\n  THE TEST: does the between-split ordering survive on internal "
          f"branches?", file=out)
    mt = [r["mt"] for r in rows]
    print(f"    correlation of each class with the divergence ratio "
          f"(mutations:tips)", file=out)
    print(f"    {'class':8s} {'internal r':>11s} {'internal rho':>13s} "
          f"{'terminal r':>11s} {'terminal rho':>13s}", file=out)
    for j, c in enumerate(cats):
        vi = [props(r["internal"])[j] for r in rows]
        vt = [props(r["terminal"])[j] for r in rows]
        print(f"    {c:8s} {pearson(mt, vi):11.3f} {spearman(mt, vi):13.3f} "
              f"{pearson(mt, vt):11.3f} {spearman(mt, vt):13.3f}", file=out)
    print("    A lineage signal keeps its internal-only correlations; a "
          "quality gradient", file=out)
    print("    shows up in the terminal column and fades in the internal one.",
          file=out)

    print(f"\n  per split, sorted by divergence ratio:", file=out)
    print(f"    {'node':12s} {'tips':>7s} {'M:T':>6s} {'%term':>6s} "
          f"{'internal_n':>11s} {'terminal_n':>11s} {'i_vs_t_chi2':>12s}",
          file=out)
    for r in sorted(rows, key=lambda r: r["mt"]):
        print(f"    {r['node']:12s} {r['tips']:7d} {r['mt']:6.1f} "
              f"{100 * r['term_frac']:6.1f} {r['ni']:11,d} {r['nt']:11,d} "
              f"{r['chi2']:12,.1f}", file=out)

    if args.output:
        with open(args.output, "w") as fh:
            fh.write("\t".join(["Node_ID", "Branch_Class", "Tips",
                                "Mutations"] + list(cats)) + "\n")
            for r in sorted(rows, key=lambda r: r["node"]):
                for label, key, n in (("internal", "internal", r["ni"]),
                                      ("terminal", "terminal", r["nt"])):
                    p = props(r[key])
                    fh.write("\t".join([r["node"], label, str(r["tips"]),
                                        str(n)]
                                       + [f"{v:.6f}" for v in p]) + "\n")
        print(f"\nPer-split spectra written to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
