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
     there in every replicate, and that clade can win the whole-tree pass
     again and again.  If one node is the whole-tree maximum in more than
     --hotspot_fraction of the replicates, a single hotspot is setting the
     threshold, and masking it is cheaper than paying for it in every run.
     Node IDs are comparable across replicates because the permutation
     changes only mutation types, never the tree.

The observed rows come from the cascade run's --min_chi 1 search, which
accepts every split, so its cascade carves out splits a production run at a
higher --min_chi would not.  Deep observed rows are therefore an
approximation to what a production run would test.  For the production
run's own splits, use split_confidence.py, which models its carving.

Usage:
    null_by_size.py --null null_sites_cascade.tsv
    null_by_size.py --null null_sites_cascade.tsv --min_chi 65 --margin 2
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
            "node": f[3], "chi2": float(f[4]),
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


def hotspot(null_rows, fraction):
    """(counts, n_replicates, flagged) for the whole-tree maximum's node."""
    first = [r for r in null_rows if r["iteration"] == 1]
    counts = collections.Counter(r["node"] for r in first)
    n = len({r["replicate"] for r in first})
    flagged = [(node, c) for node, c in counts.most_common() if n and c / n > fraction]
    return counts, n, flagged


def report(rows, args, out):
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
    counts, n, flagged = hotspot(null_rows, args.hotspot_fraction)
    print(f"  node that is the whole-tree maximum (iteration 1), across {n} "
          "replicates:", file=out)
    print("    node\treplicates", file=out)
    for node, c in counts.most_common(5):
        print(f"    {node}\t{c}", file=out)
    if flagged:
        for node, c in flagged:
            print(f"  Verdict: HOTSPOT.  {node} is the maximum in {c} of {n} "
                  "replicates", file=out)
            print(f"  (> {args.hotspot_fraction:.0%}).  A recurrent site inside "
                  "that clade is setting", file=out)
            print("  the threshold: find the top recurrent sites whose mutations "
                  "fall in it", file=out)
            print("  (profile_tree.py's SITE RECURRENCE, annotate_sites.py) "
                  "and mask them.", file=out)
    elif counts:
        node, c = counts.most_common(1)[0]
        print(f"  Verdict: no hotspot.  The most frequent node, {node}, is the "
              f"maximum in {c} of {n}", file=out)
        print(f"  replicates (limit: more than {args.hotspot_fraction:.0%}).", file=out)
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
    counts, n, flagged = hotspot(null_rows, 0.5)
    assert n == 4 and flagged == [("node_7", 3)], (n, flagged)
    _, _, flagged = hotspot(null_rows, 0.75)
    assert flagged == []
    args = argparse.Namespace(null_window=1.5, min_chi=50, margin=2, min_nulls=5,
                              hotspot_fraction=0.5)
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
    assert "HOTSPOT.  node_7 is the maximum in 3 of 4" in text
    assert "1k (whole tree)\t4\t40.0" in text, text
    print("ok: size bins, size-matched ceilings, --min_chi flags and the hotspot "
          "check behave on a constructed cascade")
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
    report(rows, args, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
