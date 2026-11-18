#!/usr/bin/env python3
"""Are the highlighted-gene hits on split branches more than chance?

This belongs to step 8 of parameter tuning, after split_drivers.py's loop has
run annotate_sites.py --highlight on each split's defining branch.  A branch
with hundreds of mutations will hit a gene list covering 1-2% of the genome
whatever the mechanism, so a hit means something only against its
expectation.  This script computes that expectation.

By default only hits that change the protein count: nonsynonymous changes,
stop gains included.  A synonymous change cannot disable a repair gene, so
counting it adds noise and, worse, can promote a mutational hotspot into a
"recurrent hit".  The expected rate is scaled to match: the genes' share of
the searchable genome times the share of coding mutations on these branches
that are nonsynonymous.  --count all counts every hit, unscaled.

Every result comes with a chance count: how many of the tests made would
reach its P by chance, from each test's own expectation (a short branch or a
short gene cannot reach a small P at all, so it adds almost nothing).  Read
that column first.  Near or above 1, the result is what the best of that
many draws looks like, and its rank among the others means nothing.

Reports:

  1. THE CHANCE RATE, from the highlighted genes' length outside the mask
     over the genome's length outside the mask.  Masked sites cannot carry a
     mutation, so the whole genome is the wrong denominator.

  2. ALL BRANCHES.  Hits against expectation summed over every branch, and
     the nonsynonymous share of all hits beside the same share for every
     coding mutation on the branches.

  3. PER BRANCH, branches with counted hits.  On a short branch any one hit
     gives P of about (mutations x rate), so the P of such a branch mostly
     reflects how short it is.

  4. ONE GENE, ONE BRANCH.  Genes with several counted hits on one branch,
     with the mutations.  Hits a few bases apart may be a single multi-base
     change or a local misalignment; check them in the reads.

  5. ONE GENE, SEVERAL BRANCHES.  Independent recurrence, the evidence a
     single branch cannot give.

A hit is a site INSIDE a highlighted gene.  annotate_sites.py --highlight_flank
also labels nearby sites as name@+N; those are not counted.

Usage:
  branch_hits.py --branches node_*_branch.tsv \\
      --annotation NC_012345.1.gff3.gz --highlight repair_genes.txt \\
      --mask selection_loci.bed,alignment_loci.bed --reference NC_012345.1.fna

  branch_hits.py --selftest
"""

import argparse
import math
import os
import re
import sys
from collections import defaultdict

from annotate_sites import (LANDMARK_TYPES, load_annotation, load_fasta,
                            merge_intervals, read_bed_intervals,
                            resolve_region_names)


def poisson_tail(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    term = math.exp(-lam)
    below = term
    for i in range(1, k):
        term *= lam / i
        below += term
    return max(0.0, 1.0 - below)


def chance_count(p, lams):
    """Expected number of tests, of Poisson means lams, reaching P <= p.

    Each test's counts are discrete, so its chance of reaching P <= p is the
    tail at the smallest count whose tail is <= p, not p itself.
    """
    total = 0.0
    for lam in lams:
        k = 0
        while poisson_tail(k, lam) > p:
            k += 1
            if k > 10000:
                break
        total += poisson_tail(k, lam)
    return total


def overlap(intervals, mask):
    """Bases of merged 1-based intervals that fall inside merged mask."""
    n = 0
    for s, e in intervals:
        for ms, me in mask:
            if me < s:
                continue
            if ms > e:
                break
            n += min(e, me) - max(s, ms) + 1
    return n


def length(intervals):
    return sum(e - s + 1 for s, e in intervals)


def read_names(path):
    names = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                names.append(line.split()[0])
    return names


def gene_intervals(features, names):
    """name -> merged intervals, plus the names the annotation lacks."""
    features = [f for f in features if f.ftype not in LANDMARK_TYPES]
    matched, missing = resolve_region_names(features, names)
    out = {}
    for name, feats in matched.items():
        out[name] = merge_intervals([(f.start, f.end) for f in feats])
    return out, missing


def read_branch(path):
    """Summarise one annotate_sites.py output file.

    hits: (gene, mutation, aa_change, nonsynonymous) per site inside a
    highlighted gene.
    """
    n = coding = nonsyn = 0
    hits = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header) - 1:
                continue
            n += 1
            aa = f[col["aa_change"]]
            is_nonsyn = "(nonsyn)" in aa
            if is_nonsyn or "(syn)" in aa:
                coding += 1
                nonsyn += is_nonsyn
            hl = f[col["highlight"]] if col["highlight"] < len(f) else ""
            if hl and "@" not in hl:
                mut = f[col["mutation"]] or f[col["position"]]
                hits.append((hl, mut, aa.split(" ")[0], is_nonsyn))
    return {"mutations": n, "coding": coding, "nonsyn": nonsyn, "hits": hits}


def branch_name(path):
    return re.sub(r"_branch\.tsv$", "", os.path.basename(path))


def analyse(branches, gene_ivs, mask, genome_size, count="nonsyn"):
    """All the numbers the report prints, as a dict."""
    masked = length([(s, min(e, genome_size)) for s, e in mask
                     if s <= genome_size])
    searchable = genome_size - masked
    union = merge_intervals([iv for ivs in gene_ivs.values() for iv in ivs])
    hl_bp = length(union)
    hl_free = hl_bp - overlap(union, mask)
    gene_free = {g: length(ivs) - overlap(ivs, mask)
                 for g, ivs in gene_ivs.items()}

    all_hits = sum(len(d["hits"]) for d in branches.values())
    hit_nonsyn = sum(h[3] for d in branches.values() for h in d["hits"])
    coding = sum(d["coding"] for d in branches.values())
    nonsyn = sum(d["nonsyn"] for d in branches.values())
    ### P(a mutation is counted | it lands in a highlighted gene).
    frac = 1.0 if count == "all" else (nonsyn / coding if coding else 0.0)
    rate = hl_free / searchable * frac

    counted = {b: [h for h in d["hits"] if count == "all" or h[3]]
               for b, d in branches.items()}
    lams = {b: d["mutations"] * rate for b, d in branches.items()}
    rows = []
    for b, d in branches.items():
        k = len(counted[b])
        genes = defaultdict(list)
        for h in counted[b]:
            genes[h[0]].append(h)
        rows.append({"branch": b, "mutations": d["mutations"], "hits": k,
                     "expected": lams[b], "p": poisson_tail(k, lams[b]),
                     "genes": dict(genes)})
    for r in rows:
        r["chance"] = chance_count(r["p"], lams.values()) if r["hits"] \
            else None
    rows.sort(key=lambda r: (r["p"], -r["hits"]))

    n_mut = sum(d["mutations"] for d in branches.values())
    n_counted = sum(len(v) for v in counted.values())
    tot = {"mutations": n_mut, "hits": n_counted, "expected": n_mut * rate,
           "all_hits": all_hits, "hit_nonsyn": hit_nonsyn, "coding": coding,
           "nonsyn": nonsyn}
    tot["p"] = poisson_tail(n_counted, tot["expected"])

    ### One gene, one branch.
    pair_lams = {(b, g): d["mutations"] * frac * bp / searchable
                 for b, d in branches.items() for g, bp in gene_free.items()}
    cluster = []
    for r in rows:
        for g, hs in r["genes"].items():
            if len(hs) < 2 or g not in gene_free:
                continue
            lam = pair_lams[(r["branch"], g)]
            p = poisson_tail(len(hs), lam)
            cluster.append({"branch": r["branch"], "gene": g, "hits": hs,
                            "expected": lam, "p": p,
                            "chance": chance_count(p, pair_lams.values())})
    cluster.sort(key=lambda r: r["p"])

    ### One gene, several branches.
    gene_lams = {g: sum(1 - math.exp(-d["mutations"] * frac * bp / searchable)
                        for d in branches.values())
                 for g, bp in gene_free.items()}
    recur = []
    for g in gene_free:
        hit_by = [(r["branch"], r["genes"][g]) for r in rows if g in r["genes"]]
        if len(hit_by) < 2:
            continue
        p = poisson_tail(len(hit_by), gene_lams[g])
        recur.append({"gene": g, "bp": gene_free[g], "hit_by": hit_by,
                      "expected": gene_lams[g], "p": p,
                      "chance": chance_count(p, gene_lams.values())})
    recur.sort(key=lambda r: r["p"])
    return {"genome": genome_size, "masked": masked, "searchable": searchable,
            "hl_bp": hl_bp, "hl_free": hl_free, "frac": frac, "rate": rate,
            "count": count, "rows": rows, "total": tot, "cluster": cluster,
            "recur": recur, "n_genes": len(gene_free)}


def pct(a, b):
    return f"{100 * a / b:.0f}%" if b else "n/a"


def hit_label(h):
    return f"{h[0]} {h[2]}" if h[2] else f"{h[0]} {h[1]}"


def print_table(table, out):
    widths = [max(len(t[i]) for t in table) for i in range(len(table[0]))]
    for t in table:
        print("  " + "  ".join(v.ljust(w) for v, w in zip(t, widths)).rstrip(),
              file=out)


def report(res, names_file, n_names, missing, out):
    t = res["total"]
    what = "nonsynonymous hits" if res["count"] == "nonsyn" else "hits"
    print(f"BRANCH HITS: genes from {names_file}; counting {what}", file=out)
    print(f"  {res['n_genes']} of {n_names} names found in the annotation"
          + (f"; not found: {', '.join(missing)}" if missing else ""), file=out)
    print(f"  highlighted genes: {res['hl_bp']:,} bp, {res['hl_free']:,} "
          f"outside the mask", file=out)
    print(f"  genome: {res['genome']:,} bp, {res['masked']:,} masked, "
          f"{res['searchable']:,} searchable", file=out)
    line = (f"  chance rate per mutation: {res['hl_free']:,} / "
            f"{res['searchable']:,}")
    if res["count"] == "nonsyn":
        line += (f" x {t['nonsyn']:,}/{t['coding']:,} nonsynonymous "
                 f"coding mutations")
    print(f"{line} = {100 * res['rate']:.3f}%", file=out)
    print(file=out)
    print("  'by chance' is how many of the tests made would reach that P by "
          "chance. Read it", file=out)
    print("  first: near or above 1, the result is what the best of that many "
          "draws looks", file=out)
    print("  like, and its rank means nothing.", file=out)
    print(file=out)

    print("ALL BRANCHES", file=out)
    print(f"  {len(res['rows'])} branches, {t['mutations']:,} mutations, "
          f"{t['hits']} {what}, {t['expected']:.1f} expected, "
          f"P(>={t['hits']}) = {t['p']:.3g}", file=out)
    print(f"  nonsynonymous share: {t['hit_nonsyn']} of {t['all_hits']} hits "
          f"({pct(t['hit_nonsyn'], t['all_hits'])}), against {t['nonsyn']:,} "
          f"of {t['coding']:,} ({pct(t['nonsyn'], t['coding'])}) coding "
          f"mutations on these branches", file=out)
    print(file=out)

    n = len(res["rows"])
    print(f"PER BRANCH ('by chance' out of {n} branches)", file=out)
    table = [["branch", "by chance", "P", what.split()[0] if
              res["count"] == "nonsyn" else "hits", "expected", "mutations",
              "genes"]]
    for r in res["rows"]:
        if not r["hits"]:
            continue
        genes = ", ".join(hit_label(h) for hs in r["genes"].values()
                          for h in hs)
        table.append([r["branch"], f"{r['chance']:.2f}", f"{r['p']:.3g}",
                      str(r["hits"]), f"{r['expected']:.2f}",
                      str(r["mutations"]), genes])
    if len(table) > 1:
        print_table(table, out)
    else:
        print("  none", file=out)
    print(file=out)

    pairs = n * res["n_genes"]
    print(f"ONE GENE, ONE BRANCH ('by chance' out of {pairs:,} branch-gene "
          f"pairs)", file=out)
    if res["cluster"]:
        table = [["branch", "gene", "by chance", "P", "hits", "expected",
                  "mutations"]]
        for r in res["cluster"]:
            table.append([r["branch"], r["gene"], f"{r['chance']:.3g}",
                          f"{r['p']:.3g}", str(len(r["hits"])),
                          f"{r['expected']:.3f}",
                          ", ".join(f"{h[1]} {h[2]}" for h in r["hits"])])
        print_table(table, out)
    else:
        print("  none", file=out)
    print(file=out)

    print(f"ONE GENE, SEVERAL BRANCHES ('by chance' out of {res['n_genes']} "
          f"genes)", file=out)
    if res["recur"]:
        table = [["gene", "by chance", "P", "branches", "expected", "bp",
                  "hits"]]
        for r in res["recur"]:
            hits = "; ".join(f"{b} " + ", ".join(h[2] or h[1] for h in hs)
                             for b, hs in r["hit_by"])
            table.append([r["gene"], f"{r['chance']:.2f}", f"{r['p']:.3g}",
                          str(len(r["hit_by"])), f"{r['expected']:.2f}",
                          f"{r['bp']:,}", hits])
        print_table(table, out)
    else:
        print("  none", file=out)


def write_tsv(res, path):
    with open(path, "w") as fh:
        fh.write("branch\tmutations\thits\texpected\tp\tby_chance\tgenes\n")
        for r in res["rows"]:
            vals = [r["branch"], str(r["mutations"]), str(r["hits"]),
                    f"{r['expected']:.4g}", f"{r['p']:.4g}",
                    "" if r["chance"] is None else f"{r['chance']:.4g}",
                    ",".join(f"{h[0]}:{h[1]}" for hs in r["genes"].values()
                             for h in hs)]
            fh.write("\t".join(vals) + "\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Test highlighted-gene hits on split branches against "
                    "chance.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--branches", nargs="+",
                   help="annotate_sites.py --highlight outputs, one per "
                        "branch, named node_N_branch.tsv")
    p.add_argument("--annotation", help="The annotation annotate_sites.py used")
    p.add_argument("--highlight", help="The gene list given to --highlight")
    p.add_argument("--mask", default="",
                   help="Comma-separated BED files of the masked regions "
                        "(annotate_sites.py --mask_bed outputs)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--reference", help="Reference FASTA, for the genome length")
    g.add_argument("--genome_size", type=int)
    p.add_argument("--seqid", default=None)
    p.add_argument("--count", choices=("nonsyn", "all"), default="nonsyn",
                   help="Count only protein-changing hits (default), or all")
    p.add_argument("--output", help="Per-branch TSV")
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    if not args.selftest:
        need = [a for a in ("branches", "annotation", "highlight")
                if not getattr(args, a)]
        if need or not (args.reference or args.genome_size):
            p.error("--branches, --annotation, --highlight and one of "
                    "--reference or --genome_size are required")
    return args


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    if args.genome_size:
        genome = args.genome_size
    else:
        seqs = load_fasta(args.reference)
        genome = len(seqs[args.seqid]) if args.seqid else \
            max(len(s) for s in seqs.values())
    names = read_names(args.highlight)
    gene_ivs, missing = gene_intervals(
        load_annotation(args.annotation, args.seqid), names)
    mask = read_bed_intervals([m for m in args.mask.split(",") if m]) \
        if args.mask else []
    branches = {branch_name(p): read_branch(p) for p in args.branches}
    res = analyse(branches, gene_ivs, mask, genome, args.count)
    if args.output:
        write_tsv(res, args.output)
    report(res, args.highlight, len(names), missing, sys.stdout)
    return 0


def selftest():
    ok = True

    def check(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    check(abs(poisson_tail(8, 3.16) - 0.0158) < 0.001,
          f"P(X>=8 | 3.16) = {poisson_tail(8, 3.16):.4f}")
    check(poisson_tail(0, 2.0) == 1.0, "P(X>=0) is 1")
    check(chance_count(0.01, [0.0001]) < 0.01, "a short branch adds ~0")
    check(abs(chance_count(0.5, [5.0]) - poisson_tail(6, 5.0)) < 1e-12,
          "discrete tail used, not p")
    check(overlap([(10, 20)], [(15, 30)]) == 6, "overlap is inclusive")

    import io
    import tempfile
    d = tempfile.mkdtemp()
    hdr = ("position\tcount\tmutation\tcontext\tfeature_type\tgene\t"
           "locus_tag\tstrand\tcodon\tcodon_pos\tref_codon\tref_aa\t"
           "aa_change\tdistance\tproduct\thighlight\n")
    rows = ["1\t\tA1G\tCDS\tCDS\tx\tR1\t+\t1\t1\tATG\tM\tM1V (nonsyn)\t0\tp\trecA\n",
            "2\t\tA2G\tCDS\tCDS\tx\tR1\t+\t1\t2\tATG\tM\tM1I (nonsyn)\t0\tp\trecA\n",
            "3\t\tA3G\tCDS\tCDS\tx\tR1\t+\t1\t3\tATG\tM\tM1M (syn)\t0\tp\trecA\n",
            "4\t\tA4G\tCDS\tCDS\ty\tR2\t+\t2\t1\tATG\tM\tM2V (nonsyn)\t0\tp\t\n",
            "5\t\tA5G\tpromoter\t\t\t\t\t\t\t\t\t\t-5\t\tmutY@-5\n"]
    path = os.path.join(d, "node_7_branch.tsv")
    with open(path, "w") as fh:
        fh.write(hdr + "".join(rows))
    b = read_branch(path)
    check(branch_name(path) == "node_7", "branch name from file name")
    check(b["mutations"] == 5 and len(b["hits"]) == 3 and b["coding"] == 4
          and b["nonsyn"] == 3, "hits inside the gene only; flank not counted")

    gene_ivs = {"recA": [(1, 100)], "mutY": [(201, 300)]}
    syn_only = {"mutations": 10, "coding": 10, "nonsyn": 5,
                "hits": [("mutY", "C250T", "L5L", False)]}
    branches = {"node_7": b, "node_9": {"mutations": 10, "coding": 10,
                                        "nonsyn": 5,
                                        "hits": [("recA", "A9G", "K3E", True),
                                                 ("mutY", "C251T", "L6L",
                                                  False)]},
                "node_80": syn_only}
    res = analyse(branches, gene_ivs, [(51, 100), (901, 1000)], 1000)
    check(res["searchable"] == 850 and res["hl_free"] == 150,
          f"mask removed from both: {res['hl_free']}/{res['searchable']}")
    check(abs(res["frac"] - 13 / 24) < 1e-12, "nonsyn fraction pooled")
    check(abs(res["rate"] - 150 / 850 * 13 / 24) < 1e-12, "rate scaled")
    check(res["total"]["hits"] == 3, "synonymous hits not counted")
    check([r["gene"] for r in res["recur"]] == ["recA"],
          "a synonymous-only gene does not recur")
    check([(r["branch"], r["gene"], len(r["hits"])) for r in res["cluster"]]
          == [("node_7", "recA", 2)], "clustering counts nonsyn only")
    res_all = analyse(branches, gene_ivs, [(51, 100), (901, 1000)], 1000,
                      "all")
    check(res_all["total"]["hits"] == 6 and res_all["frac"] == 1.0
          and {r["gene"] for r in res_all["recur"]} == {"recA", "mutY"},
          "--count all counts every hit, unscaled")
    buf = io.StringIO()
    report(res, "genes.txt", 2, [], buf)
    check("ONE GENE, SEVERAL BRANCHES" in buf.getvalue(), "report renders")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
