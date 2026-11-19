#!/usr/bin/env python3
"""Read a cascade permutation null by subtree size, and check it for a hotspot.

This belongs to step 7 of parameter tuning.  permutation_null.py's printed
summary reports each replicate's maximum, which is its first split, over the
whole tree.  Chi-squared scales with the counts in the table, so that number
is far too strict for splits deeper in the cascade, which are scored on
smaller subtrees.  The cascade TSV (--max_splits 25 --min_chi 1) holds the
null for those too; this script reads it.

Three reports:

  1. NULL CEILING BY SUBTREE SIZE.  The whole-tree pass, then the deeper
     null splits grouped into factor-of-2 bins of the size of the table each
     was scored on (mutations_below +
     mutations_above, which is exact: it is the search root minus every split
     already carved out of it).  For each bin: how many null splits, their
     size range, and the maximum and median chi-squared.

  2. OBSERVED SPLITS AGAINST A SIZE-MATCHED NULL.  Each observed split in the
     same file against the highest null chi-squared from subtrees within
     --null_window fold of its size: the same ceiling and window that
     split_confidence.py uses, so the two agree.  Given --min_chi, splits on
     the wrong side of that threshold are flagged: below it but clearing
     their null by --margin (a fixed threshold is rejecting them), or above it
     but not clearing their null at all.

  3. HOTSPOT CHECK.  Under the sites null, each site keeps its branches, so a
     heavily recurrent site inside one clade dumps the same clumped counts
     there in every replicate.  How much that scores depends on which type
     the site draws, so such a clade need not win often; it wins big when it
     draws a rare type.  Either of two signs flags it:
       - one clade is the whole-tree maximum in more than
         --hotspot_fraction of the replicates;
       - one clade holds the --top_maxima largest replicate maxima, so it
         sets the threshold the printed summary recommends.
     Node IDs are comparable across replicates because the permutation
     changes only mutation types, never the tree.  One clade can still show
     up under two IDs, a node and a child holding nearly all of it.  With
     --input_tree, nested maxima whose mutations_below agree within
     --same_clade count as one clade, and each flagged clade's most
     recurrent sites are listed with an annotate_sites.py command for them.
     Masking those sites is cheaper than paying for them in every run.

The observed rows come from the cascade run's --min_chi 1 search, which
accepts every split, so its cascade carves out splits a production run at a
higher --min_chi would not.  Deep observed rows are therefore an
approximation to what a production run would test.  For the production
run's own splits, use split_confidence.py, which models its carving.

Usage:
    null_by_size.py --null null_sites_cascade.tsv
    null_by_size.py --null null_sites_cascade.tsv --min_chi 65 --margin 2 \
        --input_tree tree.pruned.masked.pb.gz \
        --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna

--input_tree needs bte; without it the script reads text only.
--annotation and --reference are used only in the printed command.
"""

import argparse
import collections
import math
import statistics
import sys

from split_confidence import ceiling, parse_null


def parse_cascade(lines):
    """All rows as dicts, observed and null alike."""
    rows = []
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if i == 0 or len(f) < 7:
            continue
        rows.append({
            "replicate": f[0], "iteration": int(f[1]), "root": f[2],
            "node": f[3], "chi2": float(f[4]), "below": int(f[5]),
            "size": int(f[5]) + int(f[6]),
        })
    return rows


def size_bins(null_rows):
    """Group null rows into factor-of-2 bins, counted down from the largest."""
    if not null_rows:
        return []
    top = max(r["size"] for r in null_rows)
    bins = collections.defaultdict(list)
    for r in null_rows:
        k = int(math.floor(math.log2(top / r["size"]))) if r["size"] > 0 else 99
        bins[k].append(r)
    return [bins[k] for k in sorted(bins)]


def fmt_size(n):
    if n >= 1e6:
        return f"{n / 1e6:.2f}M"
    if n >= 1e3:
        return f"{n / 1e3:.0f}k"
    return str(n)


def same_clades(tree, first, same):
    """node -> the clade it counts as: itself, or the nearest ancestor among
    the whole-tree maxima, if its mutations_below is at least `same` times
    that ancestor's (a node and a child holding nearly all of it)."""
    below = {r["node"]: r["below"] for r in first}
    clade = {}
    for node in sorted(below, key=lambda n: -below[n]):
        clade[node] = node
        if tree is None:
            continue
        p = tree.get_node(node).parent
        while p is not None:
            if p.id in below:
                if below[node] >= same * below[p.id]:
                    clade[node] = clade[p.id]
                break
            p = p.parent
    return clade


def hotspot(null_rows, fraction, top_maxima, clade=None):
    """The whole-tree maxima, ranked, and the clades flagged, with why."""
    first = [r for r in null_rows if r["iteration"] == 1]
    clade = clade or {}

    def key(r):
        return clade.get(r["node"], r["node"])

    counts = collections.Counter(key(r) for r in first)
    n = len({r["replicate"] for r in first})
    ranked = sorted(first, key=lambda r: -r["chi2"])
    flagged = {}
    for c, k in counts.most_common():
        if n and k / n > fraction:
            flagged[c] = [f"it is the maximum in {k} of {n} replicates"]
    if top_maxima > 1 and len(ranked) >= top_maxima:
        lead = key(ranked[0])
        if all(key(r) == lead for r in ranked[:top_maxima]):
            rest = [r["chi2"] for r in ranked if key(r) != lead]
            why = (f"it holds the {top_maxima} largest maxima ("
                   + ", ".join(f"{r['chi2']:.1f}" for r in ranked[:top_maxima])
                   + ")" + (f"; no other clade's exceeds {max(rest):.1f}"
                            if rest else ""))
            flagged.setdefault(lead, []).append(why)
    return {"counts": counts, "n": n, "ranked": ranked, "flagged": flagged,
            "key": key}


def subtree(tree, root=None):
    """Every node at and below root (the whole tree if None)."""
    stack = [tree.get_node(root) if root else tree.root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(node.children)


def site_counts(nodes):
    counts = collections.Counter()
    for node in nodes:
        for m in node.mutations:
            if len(m) > 2 and m[0] in "ACGT" and m[-1] in "ACGT":
                counts[int(m[1:-1])] += 1
    return counts


def clade_sites(tree, node, top, whole):
    """(clade's mutations, [(position, in clade, tree-wide)]) for its most
    recurrent sites.  The clade includes the node's own branch, as
    mutations_below does."""
    inside = site_counts(subtree(tree, node))
    return sum(inside.values()), [(pos, k, whole[pos])
                                  for pos, k in inside.most_common(top)]


def report(rows, args, out, tree=None):
    null_rows = [r for r in rows if r["replicate"] != "observed"]
    observed = [r for r in rows if r["replicate"] == "observed"]
    nulls = [(r["size"], r["chi2"]) for r in null_rows]
    n_reps = len({r["replicate"] for r in null_rows})

    print("=" * 72, file=out)
    print("NULL CEILING BY SUBTREE SIZE", file=out)
    print("=" * 72, file=out)
    print(f"  {len(null_rows)} null splits from {n_reps} replicates, in factor-of-2 "
          "bins of the", file=out)
    print("  size of the table each was scored on (mutations_below + "
          "mutations_above).", file=out)
    print("  The whole-tree pass (iteration 1) is listed first; it is what the", file=out)
    print("  printed summary of permutation_null.py reports.", file=out)
    print("    size range\tnull splits\tmax chi2\tmedian chi2\tmax at", file=out)
    whole = [r for r in null_rows if r["iteration"] == 1]
    deeper = [r for r in null_rows if r["iteration"] > 1]
    groups = ([("whole tree", whole)] if whole else []) + [(None, b) for b in size_bins(deeper)]
    for label, b in groups:
        lo = min(r["size"] for r in b)
        hi = max(r["size"] for r in b)
        top = max(b, key=lambda r: r["chi2"])
        rng = fmt_size(lo) if fmt_size(lo) == fmt_size(hi) else f"{fmt_size(lo)}-{fmt_size(hi)}"
        if label:
            rng = f"{rng} ({label})"
        print(f"    {rng}\t{len(b)}\t{top['chi2']:.1f}\t"
              f"{statistics.median(r['chi2'] for r in b):.1f}\t"
              f"{top['replicate']} {top['node']}", file=out)
    print(file=out)

    if observed:
        print("=" * 72, file=out)
        print("OBSERVED SPLITS AGAINST A SIZE-MATCHED NULL", file=out)
        print("=" * 72, file=out)
        print(f"  ceiling = highest null chi2 from subtrees within {args.null_window}x "
              "of the split's", file=out)
        print("  size (split_confidence.py's rule); margin = chi2 / ceiling.", file=out)
        if args.min_chi is not None:
            print(f"  flag: 'rejected' if chi2 <= --min_chi {args.min_chi:g} but "
                  f"margin >= {args.margin:g};", file=out)
            print("        'within null' if chi2 > --min_chi but margin < 1.", file=out)
        print(f"  ('few' marks a ceiling from fewer than {args.min_nulls} null splits.)",
              file=out)
        print("    iter\tsplit_root\tnode\tsize\tchi2\tceiling\tnulls\tmargin\tflag",
              file=out)
        n_rejected = n_within = n_fail = 0
        for r in observed:
            ceil, n = ceiling(nulls, r["size"], args.null_window)
            flags = []
            if ceil is None:
                margin_s, ceil_s = "-", "-"
                margin = None
            else:
                margin = r["chi2"] / ceil
                margin_s, ceil_s = f"{margin:.1f}", f"{ceil:.1f}"
                if margin < 1:
                    n_fail += 1
            if n < args.min_nulls:
                flags.append("few")
            if args.min_chi is not None and margin is not None:
                if r["chi2"] <= args.min_chi and margin >= args.margin:
                    flags.append("rejected")
                    n_rejected += 1
                elif r["chi2"] > args.min_chi and margin < 1:
                    flags.append("within null")
                    n_within += 1
            print(f"    {r['iteration']}\t{r['root']}\t{r['node']}\t{r['size']}\t"
                  f"{r['chi2']:.1f}\t{ceil_s}\t{n}\t{margin_s}\t{','.join(flags)}",
                  file=out)
        print(f"  {n_fail} of {len(observed)} observed splits fall below their "
              "size-matched ceiling.", file=out)
        if args.min_chi is not None:
            print(f"  At --min_chi {args.min_chi:g}: {n_rejected} rejected splits clear "
                  f"their null by >= {args.margin:g}x;", file=out)
            print(f"  {n_within} accepted splits do not clear their null.", file=out)
        print(file=out)

    print("=" * 72, file=out)
    print("HOTSPOT CHECK", file=out)
    print("=" * 72, file=out)
    first = [r for r in null_rows if r["iteration"] == 1]
    clade = same_clades(tree, first, args.same_clade)
    h = hotspot(null_rows, args.hotspot_fraction, args.top_maxima, clade)
    key, n = h["key"], h["n"]
    merged = tree is not None and any(c != node for node, c in clade.items())
    print(f"  The largest whole-tree null maxima (iteration 1), of {n} "
          "replicates:", file=out)
    print("    chi2\treplicate\tnode\tmutations_below\tclade", file=out)
    for r in h["ranked"][:5]:
        print(f"    {r['chi2']:.1f}\t{r['replicate']}\t{r['node']}\t"
              f"{r['below']}\t{key(r)}", file=out)
    print("  Replicates in which each clade is the maximum:", file=out)
    print("    clade\treplicates", file=out)
    for c, k in h["counts"].most_common(5):
        print(f"    {c}\t{k}", file=out)
    if tree is None:
        print("  (Without --input_tree, a node and a child that holds nearly "
              "all of it", file=out)
        print("  count as two clades.  Equal mutations_below above is the "
              "sign.)", file=out)
    elif merged:
        print(f"  (Nested nodes whose mutations_below agree within "
              f"{1 - args.same_clade:.0%} count as one clade.)", file=out)
    print(f"  Flags: a clade that is the maximum in more than "
          f"{args.hotspot_fraction:.0%} of replicates", file=out)
    print(f"  (--hotspot_fraction), or that holds the {args.top_maxima} "
          "largest maxima (--top_maxima).", file=out)
    if not h["flagged"]:
        print("  Verdict: no hotspot.", file=out)
    whole = site_counts(subtree(tree)) if tree is not None and h["flagged"] \
        else None
    for c, why in h["flagged"].items():
        print(f"  Verdict: HOTSPOT at {c}: " + "; ".join(why) + ".", file=out)
        print("  A recurrent site inside that clade is likely setting the "
              "threshold.", file=out)
        if tree is None:
            print("  Rerun with --input_tree to list its most recurrent "
                  "sites.", file=out)
            continue
        total, sites = clade_sites(tree, c, args.hotspot_sites, whole)
        print(f"  {c}'s most recurrent sites ({total} mutations in the "
              "clade):", file=out)
        print("    position\tin clade\ttree-wide\tshare of the clade", file=out)
        for pos, k, w in sites:
            print(f"    {pos}\t{k}\t{w}\t{100 * k / total:.2f}%", file=out)
        hot = [str(pos) for pos, k, w in sites if k >= args.hotspot_min]
        if hot:
            print(f"  Sites with >= {args.hotspot_min} mutations in the clade: "
                  "annotate them, and mask them", file=out)
            print("  if they look like an artifact (step 7):", file=out)
            cmd = "annotate_sites.py"
            if args.annotation:
                cmd += f" --annotation {args.annotation}"
            if args.reference:
                cmd += f" --reference {args.reference}"
            print(f"    {cmd} --positions {','.join(hot)}", file=out)
    print(file=out)


def selftest():
    import io
    lines = ["replicate\titeration\tsplit_root_id\tnode_id\tchi2\t"
             "mutations_below\tmutations_above\n"]
    # Whole tree 1000: observed splits big; null maxima, node_7 wins 3 of 4.
    lines.append("observed\t1\tnode_1\tnode_5\t500\t300\t700\n")
    lines.append("observed\t2\tnode_1\tnode_9\t30\t100\t500\n")      # size 600
    lines.append("observed\t2\tnode_5\tnode_6\t8\t100\t200\n")       # size 300
    for rep, node, chi in (("perm1", "node_7", 40), ("perm2", "node_7", 35),
                           ("perm3", "node_7", 38), ("perm4", "node_8", 33)):
        lines.append(f"{rep}\t1\tnode_1\t{node}\t{chi}\t200\t800\n")
        lines.append(f"{rep}\t2\tnode_1\tnode_x\t{chi / 3:.1f}\t150\t550\n")
        lines.append(f"{rep}\t2\tnode_7\tnode_y\t{chi / 4:.1f}\t50\t250\n")
    rows = parse_cascade(lines)
    assert len(rows) == 15
    null_rows = [r for r in rows if r["replicate"] != "observed"]
    bins = size_bins([r for r in null_rows if r["iteration"] > 1])
    assert [len(b) for b in bins] == [4, 4], [len(b) for b in bins]
    h = hotspot(null_rows, 0.5, 2)
    assert h["n"] == 4 and list(h["flagged"]) == ["node_7"], h["flagged"]
    assert len(h["flagged"]["node_7"]) == 2, h["flagged"]   # both signs
    h = hotspot(null_rows, 0.75, 5)
    assert h["flagged"] == {}, h["flagged"]
    args = argparse.Namespace(null_window=1.5, min_chi=50, margin=2, min_nulls=5,
                              hotspot_fraction=0.5, top_maxima=2,
                              same_clade=0.95, hotspot_sites=3, hotspot_min=3,
                              annotation="ann.gff3", reference="ref.fa")
    buf = io.StringIO()
    report(rows, args, buf)
    text = buf.getvalue()
    # node_9 at 600: the 1.5x window (400-900) holds only the size-700 nulls,
    # whose max is 40/3 -> margin 2.25, below --min_chi 50, so rejected.
    line9 = next(l for l in text.splitlines() if "\tnode_9\t" in l)
    assert line9.split("\t")[5] == "13.3" and "rejected" in line9, line9
    # node_6 at 300: ceiling 10.0 (40/4), chi 8 -> fails its null.
    line6 = next(l for l in text.splitlines() if "\tnode_6\t" in l)
    assert line6.split("\t")[7] == "0.8", line6
    assert "HOTSPOT at node_7: it is the maximum in 3 of 4" in text, text
    assert "1k (whole tree)\t4\t40.0" in text, text

    # A clade that rarely wins but wins big, under two IDs: node_20 and its
    # child node_21, which holds nearly all of it.  It is the maximum in 2 of
    # 6 replicates, under the frequency limit, but holds the 2 largest.
    # Two sites inside it carry most of its mutations.
    class Node:
        def __init__(self, nid, muts, children=()):
            self.id, self.mutations = nid, list(muts)
            self.children, self.parent = list(children), None
            for c in self.children:
                c.parent = self

    class Tree:
        def __init__(self, root):
            self.root, self.nodes, stack = root, {}, [root]
            while stack:
                n = stack.pop()
                self.nodes[n.id] = n
                stack.extend(n.children)

        def get_node(self, nid):
            return self.nodes[nid]

    leaves = [Node(f"leaf{i}", ["A50G", "T52G"] if i % 2 else ["G50A"])
              for i in range(8)]
    n21 = Node("node_21", ["C9T"], leaves)
    n20 = Node("node_20", ["C7T"], [n21])
    tree = Tree(Node("node_1", [], [n20, Node("node_30", ["A50G", "G3T"])]
                     + [Node(f"node_{i}", ["G3T"]) for i in (31, 32, 33)]))
    lines = lines[:1] + [f"perm{i}\t1\tnode_1\t{node}\t{chi}\t{below}\t900\n"
                         for i, (node, chi, below) in enumerate(
                             [("node_20", 45, 100), ("node_21", 31, 98),
                              ("node_30", 24, 300), ("node_31", 22, 200),
                              ("node_32", 21, 200), ("node_33", 20, 200)])]
    rows = parse_cascade(lines)
    h = hotspot(rows, 0.5, 2)
    assert h["flagged"] == {}, "without the tree, two IDs hide the clade"
    first = [r for r in rows if r["iteration"] == 1]
    clade = same_clades(tree, first, 0.95)
    assert clade["node_21"] == "node_20" and clade["node_30"] == "node_30"
    h = hotspot(rows, 0.5, 2, clade)
    assert list(h["flagged"]) == ["node_20"], h["flagged"]
    assert "no other clade's exceeds 24.0" in h["flagged"]["node_20"][0]
    total, sites = clade_sites(tree, "node_20", 3,
                               site_counts(subtree(tree)))
    assert total == 14 and sites[0] == (50, 8, 9) and sites[1] == (52, 4, 4), \
        (total, sites)
    buf = io.StringIO()
    report(rows, args, buf, tree)
    text = buf.getvalue()
    assert "HOTSPOT at node_20: it holds the 2 largest maxima (45.0, 31.0)" \
        in text, text
    assert ("annotate_sites.py --annotation ann.gff3 --reference ref.fa "
            "--positions 50,52") in text, text
    print("ok: size bins, size-matched ceilings, --min_chi flags and the hotspot "
          "check behave on a constructed cascade")
    print("ok: a clade under two IDs that holds the largest maxima is flagged, "
          "with its recurrent sites")
    print("all self-tests passed")
    return 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Read a cascade permutation null by subtree size, and check "
                    "it for a hotspot.")
    p.add_argument("--null", help="permutation_null.py --output TSV from a "
                   "--max_splits > 1 --min_chi 1 run")
    p.add_argument("--null_window", type=float, default=1.5,
                   help="fold-range of subtree size accepted as size-matched, as in "
                        "split_confidence.py. Default 1.5.")
    p.add_argument("--min_chi", type=float, default=None,
                   help="a production --min_chi to flag observed splits against")
    p.add_argument("--margin", type=float, default=2.0,
                   help="margin over the size-matched ceiling for a split below "
                        "--min_chi to be flagged as rejected. Default 2.")
    p.add_argument("--min_nulls", type=int, default=10,
                   help="mark a ceiling drawn from fewer null splits than this. "
                        "Default 10.")
    p.add_argument("--hotspot_fraction", type=float, default=0.5,
                   help="flag a node that is the whole-tree null maximum in more "
                        "than this fraction of replicates. Default 0.5.")
    p.add_argument("--top_maxima", type=int, default=2,
                   help="flag a clade that holds this many of the largest "
                        "whole-tree null maxima. Default 2.")
    p.add_argument("--input_tree", default=None,
                   help="the tree the null was run on (needs bte): counts a "
                        "node and a child holding nearly all of it as one "
                        "clade, and lists a flagged clade's recurrent sites")
    p.add_argument("--same_clade", type=float, default=0.95,
                   help="with --input_tree, a nested node whose "
                        "mutations_below is at least this fraction of its "
                        "ancestor's counts as the same clade. Default 0.95.")
    p.add_argument("--hotspot_sites", type=int, default=10,
                   help="recurrent sites to list per flagged clade. Default 10.")
    p.add_argument("--hotspot_min", type=int, default=20,
                   help="sites with at least this many mutations in a flagged "
                        "clade go into the printed annotate_sites.py command. "
                        "Default 20.")
    p.add_argument("--annotation", default=None,
                   help="used only in the printed annotate_sites.py command")
    p.add_argument("--reference", default=None,
                   help="used only in the printed annotate_sites.py command")
    p.add_argument("--selftest", action="store_true",
                   help="run the built-in checks and exit")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    if not args.null:
        print("error: --null is required", file=sys.stderr)
        return 2
    with open(args.null) as fh:
        rows = parse_cascade(fh.readlines())
    if not any(r["iteration"] > 1 for r in rows if r["replicate"] != "observed"):
        print("warning: no null split beyond iteration 1; this looks like a "
              "--max_splits 1 run, so there is nothing to bin by size",
              file=sys.stderr)
    tree = None
    if args.input_tree:
        import bte
        print(f"Loading {args.input_tree}", file=sys.stderr)
        tree = bte.MATree(args.input_tree)
    report(rows, args, sys.stdout, tree)
    return 0


if __name__ == "__main__":
    sys.exit(main())
