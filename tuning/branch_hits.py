#!/usr/bin/env python3
"""Are the highlighted-gene hits on split branches more than chance?

This belongs to step 8 of parameter tuning, after split_drivers.py's loop has
run annotate_sites.py --highlight on each split's defining branch.  A branch
with hundreds of mutations will hit a gene list covering 1-2% of the genome
whatever the mechanism, so a hit means something only against its
expectation.  This script computes that expectation and reports:

  1. THE CHANCE RATE.  The highlighted genes' length outside the mask, over
     the genome's length outside the mask.  Masked sites cannot carry a
     mutation, so the whole genome is the wrong denominator.

  2. ALL BRANCHES.  Mutations, hits and expected hits summed over every
     branch file, with the Poisson tail P(>= observed).  Also the share of
     hits that change an amino acid, beside the same share for every coding
     mutation on the branches: a mechanism that disables genes should raise
     the first above the second.

  3. PER BRANCH.  The same per branch, sorted by P.  chance_branches is how
     many of the branches tested would reach that P by chance, computed from
     each branch's own length (a short branch cannot reach a small P at all).
     A value near or above 1 means the branch looks like the best of that
     many draws, not like a mechanism.

  4. GENES HIT MORE THAN ONCE ON ONE BRANCH.  Several hits in one gene on
     one branch is either a strong signal or one event: hits a few bases
     apart are often a single multi-base change or a local misalignment, so
     the positions are listed.  chance_pairs is the correction across every
     branch and gene pair.

  5. GENES HIT ON MORE THAN ONE BRANCH.  Independent recurrence is the
     evidence a single branch cannot give.  For each gene hit on two or more
     branches: the branches expected to hit it by chance, from the gene's
     length and each branch's mutation count, the Poisson tail, and
     chance_genes, the same correction across every gene in the list.

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
    """Summarise one annotate_sites.py output file."""
    n = coding = nonsyn = 0
    hits, hit_nonsyn = 0, 0
    genes = defaultdict(int)
    positions = defaultdict(list)
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
                hits += 1
                hit_nonsyn += is_nonsyn
                genes[hl] += 1
                positions[hl].append(f[col["mutation"]] or f[col["position"]])
    return {"mutations": n, "coding": coding, "nonsyn": nonsyn,
            "hits": hits, "hit_nonsyn": hit_nonsyn, "genes": dict(genes),
            "positions": dict(positions)}


def branch_name(path):
    base = os.path.basename(path)
    return re.sub(r"_branch\.tsv$", "", base)


def analyse(branches, gene_ivs, mask, genome_size):
    """All the numbers the report prints, as a dict."""
    masked = length([(s, min(e, genome_size)) for s, e in mask
                     if s <= genome_size])
    searchable = genome_size - masked
    union = merge_intervals([iv for ivs in gene_ivs.values() for iv in ivs])
    hl_bp = length(union)
    hl_free = hl_bp - overlap(union, mask)
    rate = hl_free / searchable
    gene_free = {g: length(ivs) - overlap(ivs, mask)
                 for g, ivs in gene_ivs.items()}

    lams = {b: d["mutations"] * rate for b, d in branches.items()}
    rows = []
    for b, d in branches.items():
        p = poisson_tail(d["hits"], lams[b])
        rows.append(dict(d, branch=b, expected=lams[b], p=p))
    for r in rows:
        r["chance"] = chance_count(r["p"], lams.values()) \
            if r["hits"] else None
    rows.sort(key=lambda r: (r["p"], -r["hits"]))

    tot = {k: sum(d[k] for d in branches.values())
           for k in ("mutations", "coding", "nonsyn", "hits", "hit_nonsyn")}
    tot["expected"] = tot["mutations"] * rate
    tot["p"] = poisson_tail(tot["hits"], tot["expected"])

    ### Recurrence: per gene, how many branches hit it, against how many
    ### would by chance given each branch's length and the gene's.
    gene_lams = {}
    for g, bp in gene_free.items():
        gene_lams[g] = sum(1 - math.exp(-d["mutations"] * bp / searchable)
                           for d in branches.values())
    recur = []
    for g in gene_free:
        hit_by = sorted(b for b, d in branches.items() if g in d["genes"])
        if len(hit_by) < 2:
            continue
        p = poisson_tail(len(hit_by), gene_lams[g])
        recur.append({"gene": g, "bp": gene_free[g], "branches": hit_by,
                      "expected": gene_lams[g], "p": p,
                      "chance": chance_count(p, gene_lams.values())})
    recur.sort(key=lambda r: r["p"])

    ### Clustering: one gene hit more than once on one branch.
    pair_lams = {(b, g): d["mutations"] * bp / searchable
                 for b, d in branches.items() for g, bp in gene_free.items()}
    cluster = []
    for b, d in branches.items():
        for g, k in d["genes"].items():
            if k < 2 or g not in gene_free:
                continue
            p = poisson_tail(k, pair_lams[(b, g)])
            cluster.append({"branch": b, "gene": g, "hits": k,
                            "mutations": d["positions"].get(g, []),
                            "expected": pair_lams[(b, g)], "p": p,
                            "chance": chance_count(p, pair_lams.values())})
    cluster.sort(key=lambda r: r["p"])
    return {"cluster": cluster, "genome": genome_size, "masked": masked, "searchable": searchable,
            "hl_bp": hl_bp, "hl_free": hl_free, "rate": rate, "rows": rows,
            "total": tot, "recur": recur, "n_genes": len(gene_free)}


def pct(a, b):
    return f"{100 * a / b:.0f}%" if b else "n/a"


def report(res, names_file, n_names, missing, out):
    t = res["total"]
    print(f"BRANCH HITS: genes from {names_file}", file=out)
    print(f"  {res['n_genes']} of {n_names} names found in the annotation"
          + (f"; not found: {', '.join(missing)}" if missing else ""), file=out)
    print(f"  highlighted genes: {res['hl_bp']:,} bp, {res['hl_free']:,} "
          f"outside the mask", file=out)
    print(f"  genome: {res['genome']:,} bp, {res['masked']:,} masked, "
          f"{res['searchable']:,} searchable", file=out)
    print(f"  chance rate per mutation: {res['hl_free']:,} / "
          f"{res['searchable']:,} = {100 * res['rate']:.3f}%", file=out)
    print(file=out)
    print("ALL BRANCHES", file=out)
    print(f"  {len(res['rows'])} branches, {t['mutations']:,} mutations, "
          f"{t['hits']} hits, {t['expected']:.1f} expected, "
          f"P(>={t['hits']}) = {t['p']:.3g}", file=out)
    print(f"  hits that change an amino acid: {t['hit_nonsyn']} of {t['hits']} "
          f"({pct(t['hit_nonsyn'], t['hits'])}), against {t['nonsyn']:,} of "
          f"{t['coding']:,} ({pct(t['nonsyn'], t['coding'])}) for every coding "
          f"mutation on these branches", file=out)
    print(file=out)
    print("PER BRANCH, branches with hits, by P (chance_branches: how many of "
          f"the {len(res['rows'])} would", file=out)
    print("reach this P by chance)", file=out)
    head = ["branch", "mutations", "hits", "nonsyn", "expected", "P",
            "chance_branches", "genes"]
    table = [head]
    for r in res["rows"]:
        if not r["hits"]:
            continue
        genes = ", ".join(f"{g}={k}" if k > 1 else g
                          for g, k in sorted(r["genes"].items()))
        table.append([r["branch"], str(r["mutations"]), str(r["hits"]),
                      str(r["hit_nonsyn"]), f"{r['expected']:.2f}",
                      f"{r['p']:.3g}", f"{r['chance']:.2f}", genes])
    print_table(table, out)
    print(file=out)
    print("GENES HIT MORE THAN ONCE ON ONE BRANCH (chance_pairs: how many of "
          f"the {len(res['rows']) * res['n_genes']:,}", file=out)
    print("branch-gene pairs would reach this P by chance)", file=out)
    if res["cluster"]:
        table = [["branch", "gene", "hits", "expected", "P", "chance_pairs",
                  "mutations"]]
        for r in res["cluster"]:
            table.append([r["branch"], r["gene"], str(r["hits"]),
                          f"{r['expected']:.3f}", f"{r['p']:.3g}",
                          f"{r['chance']:.3g}", " ".join(r["mutations"])])
        print_table(table, out)
    else:
        print("  none", file=out)
    print(file=out)
    print("GENES HIT ON MORE THAN ONE BRANCH (chance_genes: how many of the "
          f"{res['n_genes']} genes", file=out)
    print("would reach this P by chance)", file=out)
    if not res["recur"]:
        print("  none", file=out)
        return
    table = [["gene", "bp", "branches", "expected", "P", "chance_genes",
              "hit on"]]
    for r in res["recur"]:
        table.append([r["gene"], f"{r['bp']:,}", str(len(r["branches"])),
                      f"{r['expected']:.2f}", f"{r['p']:.3g}",
                      f"{r['chance']:.2f}", " ".join(r["branches"])])
    print_table(table, out)


def print_table(table, out):
    widths = [max(len(t[i]) for t in table) for i in range(len(table[0]))]
    for t in table:
        print("  " + "  ".join(v.ljust(w) for v, w in zip(t, widths)).rstrip(),
              file=out)


def write_tsv(res, path):
    cols = ["branch", "mutations", "coding", "nonsyn", "hits", "hit_nonsyn",
            "expected", "p", "chance"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols + ["genes"]) + "\n")
        for r in res["rows"]:
            vals = [r[c] for c in cols]
            vals = ["" if v is None else (f"{v:.4g}" if isinstance(v, float)
                                          else str(v)) for v in vals]
            genes = ",".join(f"{g}={k}" for g, k in sorted(r["genes"].items()))
            fh.write("\t".join(vals + [genes]) + "\n")


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
    res = analyse(branches, gene_ivs, mask, genome)
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
    ### A test that cannot reach p adds nothing; one that always can adds 1.
    check(chance_count(0.01, [0.0001]) < 0.01, "a short branch adds ~0")
    check(abs(chance_count(0.5, [5.0]) - poisson_tail(6, 5.0)) < 1e-12,
          "discrete tail used, not p")
    check(overlap([(10, 20)], [(15, 30)]) == 6, "overlap is inclusive")

    import tempfile
    d = tempfile.mkdtemp()
    hdr = ("position\tcount\tmutation\tcontext\tfeature_type\tgene\t"
           "locus_tag\tstrand\tcodon\tcodon_pos\tref_codon\tref_aa\t"
           "aa_change\tdistance\tproduct\thighlight\n")
    rows = ["1\t\tA1G\tCDS\tCDS\tx\tR1\t+\t1\t1\tATG\tM\tM1V (nonsyn)\t0\tp\trecA\n",
            "2\t\tA2G\tCDS\tCDS\tx\tR1\t+\t1\t2\tATG\tM\tM1M (syn)\t0\tp\trecA\n",
            "3\t\tA3G\tCDS\tCDS\ty\tR2\t+\t1\t3\tATG\tM\tM1I (nonsyn)\t0\tp\t\n",
            "4\t\tA4G\tpromoter\t\t\t\t\t\t\t\t\t\t-5\t\tmutY@-5\n"]
    path = os.path.join(d, "node_7_branch.tsv")
    with open(path, "w") as fh:
        fh.write(hdr + "".join(rows))
    b = read_branch(path)
    check(branch_name(path) == "node_7", "branch name from file name")
    check(b["mutations"] == 4 and b["hits"] == 2 and b["hit_nonsyn"] == 1
          and b["coding"] == 3 and b["nonsyn"] == 2,
          "hits inside the gene only; flank labels not counted")

    gene_ivs = {"recA": [(1, 100)], "mutY": [(201, 300)]}
    branches = {"node_7": b,
                "node_9": {"mutations": 10, "coding": 10, "nonsyn": 5,
                           "hits": 1, "hit_nonsyn": 1, "genes": {"recA": 1},
                           "positions": {"recA": ["A5G"]}}}
    res = analyse(branches, gene_ivs, [(51, 100), (901, 1000)], 1000)
    check(res["searchable"] == 850 and res["hl_free"] == 150,
          f"mask removed from both: {res['hl_free']}/{res['searchable']}")
    check(abs(res["rate"] - 150 / 850) < 1e-12, "rate")
    check([r["gene"] for r in res["recur"]] == ["recA"],
          "recurrence across branches")
    check([(r["branch"], r["gene"], r["mutations"]) for r in res["cluster"]]
          == [("node_7", "recA", ["A1G", "A2G"])], "clustering on one branch")
    import io
    buf = io.StringIO()
    report(res, "genes.txt", 2, [], buf)
    check("GENES HIT ON MORE THAN ONE BRANCH" in buf.getvalue(),
          "report renders")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
