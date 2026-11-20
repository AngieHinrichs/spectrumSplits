#!/usr/bin/env python3
"""Check a gene-name list against an annotation before using it.

annotate_sites.py --regions and --highlight, and split_drivers.py --highlight,
read a list of gene names and match each against the annotation's gene, Name,
gene_synonym, Alias, locus_tag and old_locus_tag fields.  A name the
annotation spells differently is reported as NOT FOUND and silently left out:
harmless for a mask you can see is smaller, easy to miss for a highlight list,
where the gene you most wanted can drop out unnoticed.  Annotations differ in
exactly these ways: NCBI's H37Rv annotation names nucS only by its locus tag
Rv1321, calls ku "mku", and keeps dinB1 only as a synonym of dinX.

This script resolves every entry the way annotate_sites.py does and reports:

  ok         one locus, matched as written.  If the annotation displays it
             under another name (a synonym or old locus tag matched), that
             name is shown, since output tables will use it.
  duplicate  the same locus as an earlier entry (mutM and fpg).  Harmless
             in a --regions list, where a name plus its locus tag keeps a
             mask working on annotations that name the gene differently.
             Not in a --highlight list: split_drivers.py counts each name as
             a gene of its own, which inflates its 'by chance' counts, and
             a site's label then depends on which name sorts first.
  several    more than one locus: probably a short or reused name.
  not found  with candidates, best first:
               - a name within --near characters that contains the entry or
                 is contained in it (tag -> tagA, ku -> mku, fpg1 -> fpg,
                 dinB -> dinB1, dinB2);
               - a locus whose product mentions the entry as a word
                 ("endonuclease NucS" for nucS), proposed by locus tag when
                 the annotation gives it no other name.
             A candidate already in the list under another entry is marked,
             since the fix may then be to drop the entry, but only if it is
             the gene meant: endMS is not the end gene.  An entry with no candidate
             needs a search by function: grep the product lines.

Exits 1 if any entry is not found or matches several loci, so it can guard
a pipeline.

Candidates are suggestions.  Check each against its product before editing
the list: a locus tag in the list resolves exactly, so a wrong one misleads
quietly.

Usage:
  check_gene_names.py --names repair_genes.txt \\
      --annotation NC_012345.1.gff3.gz

  check_gene_names.py --selftest
"""

import argparse
import collections
import re
import sys

from annotate_sites import LANDMARK_TYPES, MASKABLE_TYPES, load_annotation


ACCESSION = re.compile(r"^[A-Z]{2}_\d+(\.\d+)?$")


def read_names(path):
    """Entries in file order, comments stripped, as annotate_sites.py does."""
    names = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                names.append(line.split()[0])
    return names


class Locus(object):
    """Every feature sharing a locus tag (gene, CDS, ...), merged."""

    def __init__(self, key):
        self.key = key
        self.names = {}          # lowercased -> as written
        self.display = None
        self.locus_tag = None
        self.product = None
        self.start = self.end = None

    def add(self, f):
        for n in (f.gene, f.locus_tag) + f.aliases:
            if n:
                self.names.setdefault(n.lower(), n)
        if f.ftype == "gene" or not self.display:
            self.display = f.label()
        self.locus_tag = self.locus_tag or f.locus_tag
        self.product = self.product or f.product
        self.start = f.start if self.start is None else min(self.start, f.start)
        self.end = f.end if self.end is None else max(self.end, f.end)

    def tag_only(self):
        """True if the annotation names this locus only by locus tag."""
        return bool(self.locus_tag) and self.display == self.locus_tag

    def describe(self):
        ### Protein accessions (NCBI puts NP_215837.1 in a CDS's Name) match
        ### but would only clutter the description.
        alias = sorted(v for k, v in self.names.items()
                       if v not in (self.display, self.locus_tag)
                       and not ACCESSION.match(v))
        parts = [self.display]
        if self.locus_tag and self.locus_tag != self.display:
            parts.append(self.locus_tag)
        out = " ".join(parts)
        if alias:
            out += " (also " + ", ".join(alias) + ")"
        if self.product:
            out += f": {self.product}"
        return out


def build_loci(features):
    loci = collections.OrderedDict()
    for f in features:
        if f.ftype in LANDMARK_TYPES or f.ftype not in MASKABLE_TYPES:
            continue
        key = f.locus_tag or f.gene or f"{f.seqid}:{f.start}-{f.end}"
        loci.setdefault(key, Locus(key)).add(f)
    by_name = collections.defaultdict(list)
    for locus in loci.values():
        for n in locus.names:
            by_name[n].append(locus)
    return list(loci.values()), by_name


def candidates(entry, loci, near):
    """(locus, reason) for an entry the annotation lacks, best first."""
    e = entry.lower()
    word = re.compile(r"(?<![A-Za-z0-9])" + re.escape(entry)
                      + r"(?![A-Za-z0-9])", re.IGNORECASE)
    out, seen = [], set()
    for locus in loci:
        for n, written in locus.names.items():
            ### A name inside the entry must have 3 characters, or every
            ### two-letter name (G2 in fpg2) would be a candidate.
            if n != e and abs(len(n) - len(e)) <= near \
                    and (e in n or (n in e and len(n) >= 3)):
                if locus.key not in seen:
                    out.append((0, abs(len(n) - len(e)), written, locus,
                                f"name {written}"))
                    seen.add(locus.key)
    for locus in loci:
        if locus.key not in seen and locus.product \
                and word.search(locus.product):
            use = locus.locus_tag if locus.tag_only() else locus.display
            out.append((1, 0, use, locus, "product mentions it"))
            seen.add(locus.key)
    out.sort(key=lambda c: c[:3])
    return [(c[2], c[3], c[4]) for c in out]


def check(names, loci, by_name, near):
    rows = []
    covered = {}                 # locus key -> first entry resolving to it
    for entry in names:
        hits = by_name.get(entry.lower(), [])
        row = {"entry": entry, "loci": hits, "candidates": []}
        if len(hits) > 1:
            row["status"] = "several"
        elif hits:
            first = covered.get(hits[0].key)
            row["status"] = "duplicate" if first else "ok"
            row["of"] = first
        else:
            row["status"] = "not found"
            row["candidates"] = candidates(entry, loci, near)
        for locus in hits:
            covered.setdefault(locus.key, entry)
        rows.append(row)
    for row in rows:
        row["candidates"] = [(use, locus, why, covered.get(locus.key))
                             for use, locus, why in row["candidates"]]
    return rows


def report(rows, args, out):
    count = collections.Counter(r["status"] for r in rows)
    print(f"GENE NAME CHECK: {args.names} against {args.annotation}", file=out)
    print(f"  {len(rows)} entries: " + ", ".join(
        f"{count[s]} {s}" for s in ("ok", "duplicate", "several", "not found")
        if count[s]), file=out)
    renamed = [r for r in rows if r["status"] == "ok"
               and r["loci"][0].display.lower() != r["entry"].lower()]
    if renamed:
        print(file=out)
        print("SHOWN UNDER ANOTHER NAME (matched a synonym or locus tag; output "
              "tables use the second name)", file=out)
        for r in renamed:
            print(f"  {r['entry']}  ->  {r['loci'][0].describe()}", file=out)
    dups = [r for r in rows if r["status"] == "duplicate"]
    if dups:
        print(file=out)
        print("SAME LOCUS AS AN EARLIER ENTRY (harmless in a --regions list; "
              "in a --highlight list keep one, since split_drivers.py counts "
              "each as a gene)", file=out)
        for r in dups:
            print(f"  {r['entry']}  =  {r['of']}: {r['loci'][0].describe()}",
                  file=out)
    problems = [r for r in rows if r["status"] in ("several", "not found")]
    if problems:
        print(file=out)
        print("TO FIX", file=out)
    for r in problems:
        if r["status"] == "several":
            print(f"  {r['entry']}: matches {len(r['loci'])} loci; use the "
                  f"locus tag of the one meant", file=out)
            for locus in r["loci"]:
                print(f"      {locus.describe()}", file=out)
        else:
            cands = r["candidates"][:args.max_candidates]
            if not cands:
                print(f"  {r['entry']}: not found, no candidates; search by "
                      f"function: grep -i '<function>' on the product lines",
                      file=out)
                continue
            print(f"  {r['entry']}: not found; candidates:", file=out)
            for use, locus, why, by in cands:
                note = f"; already in the list as {by}" if by else ""
                print(f"      {use}  ({why}) {locus.describe()}{note}",
                      file=out)
            more = len(r["candidates"]) - len(cands)
            if more > 0:
                print(f"      ... {more} more (--max_candidates)", file=out)
    if not problems:
        print(file=out)
        print("  Every entry resolves to exactly one locus.", file=out)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Check a gene-name list against an annotation.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--names", help="list for --regions or --highlight")
    p.add_argument("--annotation", help="GFF3 or GenBank, as for "
                   "annotate_sites.py")
    p.add_argument("--seqid", default=None,
                   help="Only this sequence of the annotation")
    p.add_argument("--near", type=int, default=2,
                   help="Largest length difference for a name candidate "
                        "(default 2)")
    p.add_argument("--max_candidates", type=int, default=5)
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    if not args.selftest and not (args.names and args.annotation):
        p.error("--names and --annotation are required")
    return args


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    loci, by_name = build_loci(load_annotation(args.annotation, args.seqid))
    rows = check(read_names(args.names), loci, by_name, args.near)
    report(rows, args, sys.stdout)
    return 1 if any(r["status"] in ("several", "not found") for r in rows) \
        else 0


def selftest():
    import io
    import os
    import tempfile
    ok = True

    def expect(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    gff = [
        "##gff-version 3",
        "chr\tx\tregion\t1\t9000\t.\t+\t.\tID=chr",
        "chr\tx\tgene\t100\t400\t.\t+\t.\tName=fpg;gene_synonym=mutM;"
        "locus_tag=X0001",
        "chr\tx\tCDS\t100\t400\t.\t+\t0\tlocus_tag=X0001;"
        "product=formamidopyrimidine-DNA glycosylase",
        "chr\tx\tgene\t500\t900\t.\t+\t.\tName=X0002;locus_tag=X0002",
        "chr\tx\tCDS\t500\t900\t.\t+\t0\tlocus_tag=X0002;"
        "product=endonuclease NucS",
        "chr\tx\tgene\t1000\t1400\t.\t-\t.\tName=tagA;locus_tag=X0003",
        "chr\tx\tgene\t1500\t1900\t.\t+\t.\tName=mku;locus_tag=X0004",
        "chr\tx\tgene\t2000\t2400\t.\t+\t.\tName=dinX;gene_synonym=dinB1;"
        "locus_tag=X0005",
        "chr\tx\tgene\t2500\t2900\t.\t+\t.\tName=dinP;gene_synonym=dinB2;"
        "locus_tag=X0006",
        "chr\tx\tgene\t3000\t3400\t.\t+\t.\tName=Rv0001;locus_tag=X0007;"
        "product=nucleoside triphosphate pyrophosphohydrolase",
        "chr\tx\tgene\t3500\t3900\t.\t+\t.\tName=dup;locus_tag=X0008",
        "chr\tx\tgene\t4000\t4400\t.\t+\t.\tName=dup;locus_tag=X0009",
    ]
    with tempfile.NamedTemporaryFile("w", suffix=".gff3", delete=False) as fh:
        fh.write("\n".join(gff) + "\n")
        path = fh.name
    try:
        loci, by_name = build_loci(load_annotation(path))
    finally:
        os.unlink(path)
    expect(len(loci) == 9, f"region line skipped, gene+CDS merged: {len(loci)}")
    names = ["fpg", "mutM", "fpg1", "nucS", "tag", "ku", "dinB", "dinB1",
             "mazG", "dup"]
    rows = check(names, loci, by_name, 2)
    by = {r["entry"]: r for r in rows}
    expect(by["fpg"]["status"] == "ok", "exact name ok")
    expect(by["mutM"]["status"] == "duplicate" and by["mutM"]["of"] == "fpg",
           "synonym of an earlier entry is a duplicate")
    c = by["fpg1"]["candidates"]
    expect(c and c[0][0] == "fpg" and c[0][3] == "fpg",
           "fpg1 -> fpg, already covered")
    c = by["nucS"]["candidates"]
    expect(c and c[0][0] == "X0002" and c[0][2] == "product mentions it",
           f"nucS -> its locus tag via the product: {c[:1]}")
    expect(by["tag"]["candidates"][0][0] == "tagA", "tag -> tagA")
    expect(by["ku"]["candidates"][0][0] == "mku", "ku -> mku")
    expect({u for u, *_ in by["dinB"]["candidates"]} == {"dinB1", "dinB2"},
           f"dinB -> dinB1, dinB2: {[u for u, *_ in by['dinB']['candidates']]}")
    expect(by["dinB"]["candidates"][0][3] == "dinB1"
           or by["dinB"]["candidates"][1][3] == "dinB1",
           "dinB1 marked as covered by its own entry")
    expect(by["dinB1"]["status"] == "ok", "synonym match ok")
    expect(by["mazG"]["status"] == "not found"
           and not by["mazG"]["candidates"], "no candidate for mazG")
    expect(by["dup"]["status"] == "several", "a name on two loci")
    args = parse_args(["--selftest"])
    args.names, args.annotation = "list.txt", "test.gff3"
    buf = io.StringIO()
    report(rows, args, buf)
    text = buf.getvalue()
    expect("dinB1  ->  dinX X0005" in text,
           "synonym match shown under the displayed name")
    expect("search by function" in text and "use the locus tag" in text,
           "report renders")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
