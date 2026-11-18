#!/usr/bin/env python3
"""What is each split made of?

This belongs to step 11 of parameter tuning.  Steps 9 and 10 ask whether a
split is robust; neither asks what is inside it, and a tier A split can still
be a third laboratory strain.  This script lists, for every split worth
interpreting, the clades of the samples its spectrum actually covers.

That is the split's RESIDUAL: its samples minus those of every split nested
inside it, because spectrumSplits.py carves nested splits out of their
parent's row.  matUtils extract -I lists every sample below a node, so a
parent counted that way borrows its children's clades (all of lineage 2's
node would report the Beijing samples that belong to the Beijing split's own
row).  Nesting comes from sample membership, as in split_confidence.py.

Per split: tier, residual samples, the number of distinct clades, and the
--top clades with their shares.  Flags:

  mixed   --mixed_share or more of the residual lies outside the top
          clade's own lineage, so the spectrum is an average over more than
          one population.  "Own lineage" goes by the dotted names: the top
          clade, its sublineages and its ancestors (lineage2.2.1 owns
          lineage2.2.1.1 and lineage2.2), since a label at another level of
          resolution is the same population, not a second one
  lab     only with --lab_pattern, a regex for the clade labels of
          laboratory strains in YOUR tree (there is no default: which
          labels those are depends on the organism and the annotation,
          e.g. 'BCG' for the MTBC vaccine strain).  A clade matching it
          holds --lab_share or more, so part of the spectrum accumulated in
          culture rather than in hosts

Inputs:
  --clade_samples  node<TAB>sample, from matUtils extract -I per split (step 10)
  --clades         sample<TAB>clade, from matUtils summary --sample-clades
  --confidence     split_confidence.tsv (step 10), for tiers and order

Usage:
  split_composition.py --clade_samples split_clade_samples.tsv \\
      --clades clades.txt --confidence split_confidence.tsv \\
      --lab_pattern BCG --output split_composition.tsv

  split_composition.py --selftest
"""

import argparse
import collections
import re
import sys

from split_confidence import nesting_from_samples, read_clade_samples

NO_CLADE = "(none)"


def read_clades(lines):
    """sample -> clade; 'None' and a header line are skipped."""
    out = {}
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or (i == 0 and f[0] == "sample"):
            continue
        if f[1] and f[1] != "None":
            out[f[0]] = f[1]
    return out


def read_confidence(lines):
    """node -> tier, in file order (strongest chi-squared first)."""
    tiers = {}
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if i == 0 or len(f) < 2:
            continue
        tiers[f[0]] = f[-1]
    return tiers


def residuals(sets):
    """node -> its samples minus those of every split nested inside it."""
    inside = nesting_from_samples(sets)
    out = {}
    for node, samples in sets.items():
        nested = set()
        for m in inside.get(node, ()):
            nested |= sets[m]
        out[node] = samples - nested
    return out


def related(a, b):
    """True if one dotted clade name is the other or an ancestor of it."""
    return a == b or a.startswith(b + ".") or b.startswith(a + ".")


def compose(samples, clades):
    counts = collections.Counter(clades.get(s, NO_CLADE) for s in samples)
    return counts.most_common()


def build_rows(sets, clades, tiers, args):
    res = residuals(sets)
    lab = re.compile(args.lab_pattern) if args.lab_pattern else None
    want = set(args.tiers.split(",")) if args.tiers != "all" else None
    order = list(tiers) + sorted(n for n in sets if n not in tiers)
    rows = []
    for node in order:
        if node not in sets:
            continue
        tier = tiers.get(node, "root" if node == "node_1" else "")
        if want is not None and tier not in want:
            continue
        samples = res[node]
        comp = compose(samples, clades)
        n = len(samples)
        flags = []
        named = [(c, k) for c, k in comp if c != NO_CLADE]
        other = sum(k for c, k in named[1:] if not related(c, named[0][0])) \
            if named else 0
        if n and other / n >= args.mixed_share:
            flags.append("mixed")
        lab_n = sum(k for c, k in comp if lab and c != NO_CLADE
                    and lab.search(c))
        if lab and n and lab_n / n >= args.lab_share:
            flags.append("lab")
        rows.append({"node": node, "tier": tier, "samples": len(sets[node]),
                     "residual": n, "clades": len(named), "other": other,
                     "top": comp[:args.top], "lab": lab_n,
                     "flags": ",".join(flags)})
    return rows


def fmt_top(top, n):
    return ", ".join(f"{c} {100 * k / n:.0f}%" for c, k in top) if n else ""


def report(rows, args, out):
    want = "all tiers" if args.tiers == "all" else f"tiers {args.tiers}"
    print(f"SPLIT COMPOSITION ({want}; the samples each split's own spectrum "
          f"covers, nested splits removed)", file=out)
    print(f"  mixed: >= {args.mixed_share:.0%} of the residual outside the top "
          f"clade's own lineage ('other')", file=out)
    if not args.lab_pattern:
        print("  lab:   not tested (give --lab_pattern)", file=out)
    else:
        print(f"  lab:   clades matching /{args.lab_pattern}/ >= "
              f"{args.lab_share:.0%} of the residual", file=out)
    print(file=out)
    table = [["node", "tier", "samples", "residual", "clades", "other",
              f"top {args.top} clades", "flags"]]
    for r in rows:
        other = f"{100 * r['other'] / r['residual']:.0f}%" if r["residual"] \
            else ""
        table.append([r["node"], r["tier"], f"{r['samples']:,}",
                      f"{r['residual']:,}", str(r["clades"]), other,
                      fmt_top(r["top"], r["residual"]), r["flags"]])
    widths = [max(len(t[i]) for t in table) for i in range(len(table[0]))]
    for t in table:
        print("  " + "  ".join(v.ljust(w) for v, w in zip(t, widths)).rstrip(),
              file=out)
    flagged = [r for r in rows if r["flags"]]
    print(file=out)
    print(f"  {len(rows)} splits, {len(flagged)} flagged"
          + (": " + " ".join(f"{r['node']} ({r['flags']})" for r in flagged)
             if flagged else ""), file=out)


def write_tsv(rows, path):
    with open(path, "w") as fh:
        fh.write("node\ttier\tsamples\tresidual\tclades\tother_lineage\t"
                 "lab_samples\tflags\tcomposition\n")
        for r in rows:
            comp = ",".join(f"{c}:{k}" for c, k in r["top"])
            fh.write(f"{r['node']}\t{r['tier']}\t{r['samples']}\t"
                     f"{r['residual']}\t{r['clades']}\t{r['other']}\t"
                     f"{r['lab']}\t"
                     f"{r['flags']}\t{comp}\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Clade composition of each split's residual samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clade_samples", help="node<TAB>sample (step 10)")
    p.add_argument("--clades", help="sample<TAB>clade from matUtils summary "
                   "--sample-clades")
    p.add_argument("--confidence", default=None,
                   help="split_confidence.tsv, for tiers and order")
    p.add_argument("--tiers", default="A,B",
                   help="Comma-separated tiers to list, or 'all' (default A,B; "
                        "needs --confidence unless 'all')")
    p.add_argument("--top", type=int, default=3)
    p.add_argument("--mixed_share", type=float, default=0.2)
    p.add_argument("--lab_pattern", default=None,
                   help="Regex for the clade labels of laboratory strains in "
                        "this tree, e.g. 'BCG' for MTBC; no default, and the "
                        "lab flag is not tested without it")
    p.add_argument("--lab_share", type=float, default=0.01)
    p.add_argument("--output", default=None)
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    if not args.selftest:
        if not (args.clade_samples and args.clades):
            p.error("--clade_samples and --clades are required")
        if args.tiers != "all" and not args.confidence:
            p.error("--confidence is needed to select tiers; or --tiers all")
    return args


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    with open(args.clade_samples) as fh:
        sets = read_clade_samples(fh)
    with open(args.clades) as fh:
        clades = read_clades(fh)
    tiers = {}
    if args.confidence:
        with open(args.confidence) as fh:
            tiers = read_confidence(fh)
    rows = build_rows(sets, clades, tiers, args)
    if args.output:
        write_tsv(rows, args.output)
    report(rows, args, sys.stdout)
    return 0


def selftest():
    ok = True

    def check(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    ### node_1 holds node_2, which holds node_3.  node_3 is all La3; node_2's
    ### residual is La1 plus some BCG; node_1's residual is lineage 4.
    sets = {"node_1": {f"s{i}" for i in range(20)},
            "node_2": {f"s{i}" for i in range(10)},
            "node_3": {f"s{i}" for i in range(3)}}
    clades = read_clades(["sample\tannotation_1\n", "s0\tLa3\n", "s1\tLa3\n",
                          "s2\tLa3\n", "s3\tLa1\n", "s4\tLa1\n", "s5\tLa1\n",
                          "s6\tLa1\n", "s7\tLa1\n", "s8\tLa1.2.BCG\n",
                          "s9\tLa1.2.BCG\n", "s10\tNone\n"]
                         + [f"s{i}\tlineage4\n" for i in range(11, 20)])
    check(clades.get("s10") is None and "sample" not in clades,
          "header and None skipped")
    res = residuals(sets)
    check(res["node_2"] == {f"s{i}" for i in range(3, 10)},
          "nested split removed from its parent")
    check(len(res["node_1"]) == 10, "root keeps only the unsplit samples")
    tiers = read_confidence(["node\ttier\n", "node_2\tA\n", "node_3\tC\n"])
    args = parse_args(["--selftest"])
    rows = build_rows(sets, clades, tiers, args)
    check([r["node"] for r in rows] == ["node_2"],
          "tier C and the root skipped by default")
    check(rows[0]["flags"] == "" and rows[0]["lab"] == 0,
          "no lab flag without --lab_pattern")
    args.lab_pattern = "BCG"
    rows = build_rows(sets, clades, tiers, args)
    r = rows[0]
    check(r["top"][0] == ("La1", 5) and r["lab"] == 2 and r["clades"] == 2,
          f"composition {r['top']}")
    check(r["flags"] == "lab" and r["other"] == 0,
          f"BCG under La1 is La1's own lineage, so lab not mixed: "
          f"'{r['flags']}'")
    check(related("lineage2.2.1", "lineage2.2.1.1")
          and related("lineage6", "lineage6.1.1")
          and not related("lineage2.2.1", "lineage2.2.10")
          and not related("La1.8.1", "La1.7.1"), "dotted-name lineage")
    args.tiers = "all"
    rows = build_rows(sets, clades, tiers, args)
    root = [r for r in rows if r["node"] == "node_1"][0]
    check(root["tier"] == "root" and root["top"][0] == ("lineage4", 9)
          and (NO_CLADE, 1) in root["top"], "root row and unlabelled samples")
    import io
    buf = io.StringIO()
    report(rows, args, buf)
    check("flagged: node_2 (lab)" in buf.getvalue(), "report renders")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
