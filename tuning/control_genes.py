#!/usr/bin/env python3
"""Draw a length-matched random control gene list for gene_hit_contrasts.py.

Step 12's synonymous control shows the comparison itself is unbiased.  It
cannot show that a nonsynonymous excess belongs to the genes tested rather
than to nonsynonymous changes anywhere.  This draws ordinary genes to run the
same test on: if they show the same excess, the candidate list's signal is
weaker than it looks; if they stay flat, the candidates stand out against
them as well as against chance.

For each gene in --match, --per_gene protein-coding genes whose CDS length is
within --tolerance of it are drawn at random (the closest ones if too few
qualify), never the same gene twice.  A gene is not eligible if:

  - it is named in --match or any --exclude list (any name or locus tag);
  - more than --max_masked of its CDS lies in the --mask BED intervals,
    since a masked site cannot carry a hit;
  - its product matches --exclude_product, a regular expression for
    functions that should not be in a control (DNA metabolism, say).

The output is a gene list in the usual format, one locus tag per line, each
with a comment naming the gene it was matched to; give it to
gene_hit_contrasts.py --genes.  --seed makes the draw reproducible.

Usage:
  control_genes.py --annotation NC_012345.1.gff3.gz \\
      --match repair_genes.txt --exclude resistance_genes.txt \\
      --mask selection_loci.bed,alignment_loci.bed \\
      --exclude_product '(?i)(DNA|repair|replicat)' \\
      --output control_genes.txt

  control_genes.py --selftest
"""

import argparse
import bisect
import random
import re
import sys

from annotate_sites import (LANDMARK_TYPES, load_annotation,
                            resolve_region_names)


def read_names(path):
    names = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                names.append(line.split()[0])
    return names


def read_bed(paths):
    """Merged 1-based inclusive intervals, all chromosomes pooled."""
    ivs = []
    for path in paths:
        with open(path) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) >= 3 and not line.startswith(("#", "track")):
                    ivs.append((int(f[1]) + 1, int(f[2])))
    ivs.sort()
    out = []
    for s, e in ivs:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def masked_bases(ivs, starts, s, e):
    n = 0
    i = max(bisect.bisect_right(starts, s) - 1, 0)
    while i < len(ivs) and ivs[i][0] <= e:
        a, b = ivs[i]
        n += max(0, min(b, e) - max(a, s) + 1)
        i += 1
    return n


def cds_loci(features):
    """locus key -> (first CDS feature, CDS length), protein-coding only."""
    out = {}
    for f in features:
        if f.ftype != "CDS" or f.ftype in LANDMARK_TYPES:
            continue
        key = f.locus_tag or f.gene or f"{f.seqid}:{f.start}-{f.end}"
        if key not in out:
            out[key] = (f, f.end - f.start + 1)
    return out


def draw(features, match, exclude, mask, exclude_product, per_gene,
         tolerance, max_masked, seed):
    """[(control locus, its length, product, matched name, matched length)]"""
    loci = cds_loci(features)
    starts = [s for s, _e in mask]
    pattern = re.compile(exclude_product) if exclude_product else None

    def keys_for(names):
        matched, _missing = resolve_region_names(features, names)
        out = {}
        for name, feats in matched.items():
            for f in feats:
                if f.ftype == "CDS":
                    out.setdefault(name, set()).add(
                        f.locus_tag or f.gene or
                        f"{f.seqid}:{f.start}-{f.end}")
        return out

    targets = keys_for(match)
    barred = set()
    for ks in targets.values():
        barred |= ks
    for ks in keys_for(exclude).values():
        barred |= ks
    pool = []
    for key, (f, length) in loci.items():
        if key in barred:
            continue
        if pattern and f.product and pattern.search(f.product):
            continue
        if mask and masked_bases(mask, starts, f.start, f.end) \
                > max_masked * length:
            continue
        pool.append((key, length, f.product or ""))
    pool.sort()
    rng = random.Random(seed)
    used = set()
    out = []
    for name in match:
        if name not in targets:
            continue
        key = sorted(targets[name])[0]
        tlen = loci[key][1]
        free = [p for p in pool if p[0] not in used]
        near = [p for p in free if abs(p[1] - tlen) <= tolerance * tlen]
        if len(near) >= per_gene:
            chosen = rng.sample(near, per_gene)
        else:
            chosen = sorted(free, key=lambda p: (abs(p[1] - tlen), p[0]))[
                :per_gene]
        for key_c, length, product in chosen:
            used.add(key_c)
            out.append((key_c, length, product, name, tlen))
    return out, len(pool)


def write_list(rows, args, n_pool, fh):
    fh.write(f"# Length-matched random control genes, drawn by "
             f"control_genes.py --seed {args.seed}\n")
    fh.write(f"# {args.per_gene} per gene of {args.match}, CDS length within "
             f"{args.tolerance:.0%}, from {n_pool:,} eligible genes\n")
    for key, length, product, name, tlen in rows:
        fh.write(f"{key}\t# for {name} ({tlen:,} bp): {length:,} bp, "
                 f"{product}\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Draw length-matched random control genes.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotation", help="GFF3 or GenBank, as for "
                   "annotate_sites.py")
    p.add_argument("--seqid", default=None)
    p.add_argument("--match", help="gene list to match lengths to")
    p.add_argument("--exclude", default="",
                   help="Comma-separated further gene lists to keep out")
    p.add_argument("--mask", default="",
                   help="Comma-separated BED files of masked intervals")
    p.add_argument("--max_masked", type=float, default=0.1,
                   help="Largest masked fraction of a CDS (default 0.1)")
    p.add_argument("--exclude_product", default=None,
                   help="Regular expression; genes whose product matches are "
                        "not eligible")
    p.add_argument("--per_gene", type=int, default=3)
    p.add_argument("--tolerance", type=float, default=0.2,
                   help="Largest relative CDS length difference (default "
                        "0.2)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--output", default=None, help="default stdout")
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    if not args.selftest and not (args.annotation and args.match):
        p.error("--annotation and --match are required")
    return args


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    features = load_annotation(args.annotation, args.seqid)
    exclude = []
    for path in [x for x in args.exclude.split(",") if x]:
        exclude += read_names(path)
    mask = read_bed([x for x in args.mask.split(",") if x])
    rows, n_pool = draw(features, read_names(args.match), exclude, mask,
                        args.exclude_product, args.per_gene, args.tolerance,
                        args.max_masked, args.seed)
    if args.output:
        with open(args.output, "w") as fh:
            write_list(rows, args, n_pool, fh)
    else:
        write_list(rows, args, n_pool, sys.stdout)
    far = [r for r in rows if abs(r[1] - r[4]) > args.tolerance * r[4]]
    print(f"{len(rows)} control genes for {len(set(r[3] for r in rows))} "
          f"genes, from {n_pool:,} eligible"
          + (f"; {len(far)} outside the length tolerance" if far else ""),
          file=sys.stderr)
    return 0


def selftest():
    from annotate_sites import Feature
    ok = True

    def check(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    feats = [Feature("c", "region", 1, 100000, 1, None, None, None)]
    ### Target T1 of 1,000 bp; candidates of assorted lengths.
    spec = [("T1", 1000, "target"), ("A", 1050, "kinase"),
            ("B", 950, "transporter"), ("C", 1100, "DNA helicase"),
            ("D", 990, "hypothetical"), ("E", 3000, "synthase"),
            ("F", 1010, "masked one"), ("X", 1000, "excluded")]
    pos = 1
    for tag, length, product in spec:
        feats.append(Feature("c", "gene", pos, pos + length - 1, 1, tag, tag,
                             None))
        feats.append(Feature("c", "CDS", pos, pos + length - 1, 1, tag, tag,
                             product))
        if tag == "F":
            fstart = pos
        pos += length + 100
    mask = [(fstart, fstart + 500)]
    rows, n_pool = draw(feats, ["T1"], ["X"], mask, "(?i)DNA", 2, 0.2, 0.1,
                        1)
    chosen = {r[0] for r in rows}
    check(n_pool == 4, f"T1, X, the DNA product and the masked gene are not "
                       f"eligible: pool {n_pool}")
    check(len(rows) == 2 and chosen <= {"A", "B", "D"},
          f"two drawn within 20% of 1,000 bp: {sorted(chosen)}")
    again, _ = draw(feats, ["T1"], ["X"], mask, "(?i)DNA", 2, 0.2, 0.1, 1)
    check(rows == again, "the same seed gives the same draw")
    rows, _ = draw(feats, ["T1"], ["X"], mask, "(?i)DNA", 4, 0.2, 0.1, 1)
    check(len(rows) == 4 and "E" in {r[0] for r in rows},
          "too few within tolerance: the closest fill in")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
