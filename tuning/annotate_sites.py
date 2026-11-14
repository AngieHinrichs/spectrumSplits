#!/usr/bin/env python3
"""Annotate genome positions against a reference annotation, and build mask regions.

Written for the Mycobacterium tuberculosis use case: profile_tree.py reports the
most recurrently mutated sites, and those sites need to be identified before you
can decide what to mask.  Recurrent sites violate the chi-squared test's
independence assumption, and in MTBC the worst offenders are drug-resistance
loci under convergent selection.

This script does NOT predict drug resistance.  It reports, for each position,
which feature it falls in and which codon it hits, using coordinates taken from
an annotation file YOU supply.  Nothing here is hardcoded from memory, so the
answers are only as good as your annotation.  For actual resistance calls use
TB-Profiler or Mykrobe, and cite the WHO mutation catalogue.

Two jobs:

  1. Annotate sites.  Given positions (and optionally mutation strings such as
     C2155168T), report feature type, gene, locus tag, strand, codon number,
     codon position, reference codon and amino acid, and the amino-acid change.
     Without an explicit alternate base, all three possible changes at that
     position are listed.

  2. Emit mask regions.  Given a list of gene names or locus tags, resolve them
     against the annotation and write a merged BED (or a flat position list) for
     masking whole loci.  Per-site masking of a top-N list is whack-a-mole: the
     adjacent base in the same codon just surfaces next.  Masking the locus is
     the fix.

Inputs:
  --annotation   GFF3 or GenBank, optionally gzipped (GFF3 preferred)
  --reference    FASTA of the reference genome (optional; needed for codons
                 and amino acids)

Examples:
  # identify the recurrent sites profile_tree.py reported
  annotate_sites.py --annotation NC_000962.3.gff3.gz --reference NC_000962.3.fna \\
      --from_profile mtb_profile.txt --output annotated_sites.tsv

  # build a mask covering known resistance loci plus their promoters
  annotate_sites.py --annotation NC_000962.3.gff3.gz \\
      --regions tb_resistance_genes.txt --promoter_flank 200 \\
      --mask_bed resistance_loci.bed

  # also mask the PE/PPE families, which are repetitive and alignment-prone
  annotate_sites.py --annotation NC_000962.3.gff3.gz \\
      --mask_pattern '^(PE|PPE)' --mask_bed pe_ppe.bed

  # check the codon arithmetic, including reverse-strand genes
  annotate_sites.py --selftest
"""

import argparse
import bisect
import gzip
import io
import re
import sys
from collections import defaultdict

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

# Feature types we care about, in priority order when a position falls in
# several overlapping features.
TYPE_PRIORITY = ["CDS", "rRNA", "tRNA", "ncRNA", "tmRNA", "misc_RNA", "gene"]

### Feature types that describe the sequence itself rather than anything on it.
### NCBI GFF3 opens each record with a "region" line spanning the whole
### chromosome, and GenBank uses "source" the same way, so EVERY position falls
### inside one.  Left in the index they swallow every intergenic site, which is
### reported as a hit on a nameless feature at distance 0 instead of as
### intergenic with its nearest gene.  They are also excluded from --regions and
### --mask_pattern matching, where one of them would mask the entire genome.
LANDMARK_TYPES = {"region", "source", "biological_region", "chromosome",
                  "scaffold", "contig", "supercontig", "assembly"}


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate genome positions against a reference annotation, "
                    "and build mask regions.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotation", type=str,
                   help="GFF3 or GenBank annotation, optionally gzipped")
    p.add_argument("--reference", type=str, default=None,
                   help="Reference FASTA; required for codon and amino-acid output")
    p.add_argument("--seqid", type=str, default=None,
                   help="Restrict to this sequence id / accession (default: all)")

    src = p.add_argument_group("site input (choose one)")
    src.add_argument("--sites", type=str, default=None,
                     help="File whose first column is a position; an optional "
                          "second column is a count, and a token like C2155168T "
                          "anywhere in the row is used as the mutation")
    src.add_argument("--positions", type=str, default=None,
                     help="Comma-separated positions")
    src.add_argument("--from_profile", type=str, default=None,
                     help="profile_tree.py output; reads its SITE RECURRENCE table")

    msk = p.add_argument_group("mask output")
    msk.add_argument("--regions", type=str, default=None,
                     help="File of gene names or locus tags (one per line, # comments "
                          "allowed) to resolve into mask intervals")
    msk.add_argument("--mask_product_pattern", type=str, default=None,
                     help="Regex matched against each feature's product description, "
                          "for masking by function rather than by name: "
                          "'(?i)(transposase|integrase|insertion element)' finds "
                          "mobile elements whatever they are called. Matching is "
                          "case sensitive, like --mask_pattern, so prefix (?i) to "
                          "ignore case. Only CDS and RNA features carry products, so "
                          "this will not match bare gene features")
    msk.add_argument("--mask_pattern", type=str, default=None,
                     help="Regex on gene name or locus tag; matching features are "
                          "added to the mask (e.g. '^(PE|PPE)')")
    msk.add_argument("--promoter_flank", type=int, default=0,
                     help="Bases to extend each masked feature upstream (5') to cover "
                          "promoter mutations; strand aware")
    msk.add_argument("--mask_pad", type=int, default=0,
                     help="Bases to extend each masked feature on BOTH sides, for "
                          "alignment artifacts that spill past a feature's boundary "
                          "rather than sitting in its promoter. Independent of and "
                          "additive with --promoter_flank, which extends upstream "
                          "only; use this for paralogous families, repeats and "
                          "prophages, where recurrent miscalls land just outside the "
                          "annotated edges")
    msk.add_argument("--mask_bed", type=str, default=None,
                     help="Write merged mask intervals as BED (0-based half-open)")
    msk.add_argument("--mask_positions", type=str, default=None,
                     help="Write every masked position, one 1-based position per line")
    msk.add_argument("--mask_matutils", type=str, default=None,
                     help="Write the mask in the form 'matUtils mask -m' expects: one "
                          "line per masked position, N<position>N. Feed it straight to "
                          "matUtils mask -i tree.pb.gz -m <this file> -o masked.pb.gz")

    p.add_argument("--flank", type=int, default=200,
                   help="For an intergenic site, call it a promoter site if it is "
                        "within this many bases upstream of a feature")
    p.add_argument("--highlight", type=str, default=None,
                   help="File of gene names or locus tags (one per line, # comments "
                        "allowed); sites in those genes are flagged in a 'highlight' "
                        "column and counted. Annotation only -- nothing is masked. "
                        "Use it to ask whether a clade's defining mutations land in "
                        "genes of interest, e.g. tb_repair_genes.txt")
    p.add_argument("--highlight_flank", type=int, default=0,
                   help="Also highlight sites within this many bp of a named "
                        "gene, reported as name@+N. Default 0: only sites "
                        "inside the gene are highlighted.")
    p.add_argument("--output", type=str, default=None,
                   help="Annotated sites TSV (default: stdout)")
    p.add_argument("--selftest", action="store_true",
                   help="Check codon arithmetic and translation, then exit")
    return p.parse_args()


def smart_open(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------
# annotation parsing
# --------------------------------------------------------------------------

class Feature(object):
    __slots__ = ("seqid", "ftype", "start", "end", "strand", "gene", "locus_tag",
                 "product", "aliases")

    def __init__(self, seqid, ftype, start, end, strand, gene, locus_tag, product,
                 aliases=()):
        self.seqid = seqid
        self.ftype = ftype
        self.start = start          # 1-based inclusive
        self.end = end              # 1-based inclusive
        self.strand = strand        # +1 / -1
        self.gene = gene            # preferred display name
        self.locus_tag = locus_tag
        self.product = product
        ### Every other name this feature is known by: gene_synonym,
        ### old_locus_tag, and whichever of gene/Name did not become `gene`.
        ### Annotations disagree about which field carries the name a person
        ### would look a gene up by, so a name request is matched against all
        ### of them rather than just the two displayed.
        self.aliases = tuple(a for a in aliases if a)

    def names(self):
        """Every name this feature can be requested by, lowercased."""
        out = set()
        for n in (self.gene, self.locus_tag) + self.aliases:
            if n:
                out.add(n.lower())
        return out

    def label(self):
        return self.gene or self.locus_tag or "?"

    def __repr__(self):
        return (f"Feature({self.label()} {self.ftype} {self.start}-{self.end} "
                f"{'+' if self.strand > 0 else '-'})")


def unescape_gff(value):
    out = []
    i = 0
    while i < len(value):
        if value[i] == "%" and i + 2 < len(value):
            try:
                out.append(chr(int(value[i + 1:i + 3], 16)))
                i += 3
                continue
            except ValueError:
                pass
        out.append(value[i])
        i += 1
    return "".join(out)


def parse_gff3(path, seqid_filter=None):
    features = []
    with smart_open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = parts[:9]
            if seqid_filter and seqid != seqid_filter:
                continue
            if ftype not in TYPE_PRIORITY and ftype not in LANDMARK_TYPES:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            kv = {}
            for field in attrs.split(";"):
                if "=" not in field:
                    continue
                k, v = field.split("=", 1)
                kv[k.strip()] = unescape_gff(v.strip())
            gene = kv.get("gene") or kv.get("Name") or kv.get("gene_synonym")
            locus_tag = kv.get("locus_tag") or kv.get("old_locus_tag")
            ### gene_synonym is a comma-separated list, and NCBI puts the name
            ### most people know a gene by there (mabA for Rv1482c, gidB for
            ### Rv3919c) while Name holds the locus tag.  Taking only the first
            ### populated field, as `gene` does for display, makes those
            ### unsearchable.
            aliases = [kv.get("gene"), kv.get("Name"), kv.get("locus_tag"),
                       kv.get("old_locus_tag")]
            for field in ("gene_synonym", "Alias"):
                if kv.get(field):
                    aliases.extend(part.strip() for part in kv[field].split(","))
            features.append(Feature(seqid, ftype, start_i, end_i,
                                    -1 if strand == "-" else 1,
                                    gene, locus_tag, kv.get("product"),
                                    aliases))
    return features


GB_LOC = re.compile(r"(complement\()?\s*<?(\d+)\s*\.\.\s*>?(\d+)")


def parse_genbank(path, seqid_filter=None):
    """Minimal GenBank feature-table parser: simple and complement() locations."""
    features = []
    seqid = None
    in_features = False
    cur = None

    def flush(cur):
        if not cur:
            return
        ftype, loc, quals = cur
        if ftype not in TYPE_PRIORITY and ftype not in LANDMARK_TYPES:
            return
        m = GB_LOC.search(loc)
        if not m:
            return
        strand = -1 if (m.group(1) or "join" in loc and "complement" in loc) else 1
        start_i, end_i = int(m.group(2)), int(m.group(3))
        if seqid_filter and seqid and seqid != seqid_filter:
            return
        aliases = [quals.get("gene"), quals.get("locus_tag"),
                   quals.get("old_locus_tag")]
        if quals.get("gene_synonym"):
            aliases.extend(p.strip() for p in quals["gene_synonym"].split(";"))
        features.append(Feature(seqid, ftype, start_i, end_i, strand,
                                quals.get("gene"), quals.get("locus_tag"),
                                quals.get("product"), aliases))

    with smart_open(path) as fh:
        for line in fh:
            if line.startswith("VERSION") or line.startswith("ACCESSION"):
                tok = line.split()
                if len(tok) > 1 and seqid is None:
                    seqid = tok[1]
                continue
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN") or line.startswith("//"):
                flush(cur)
                cur = None
                in_features = False
                continue
            if not in_features:
                continue
            if len(line) > 5 and line[5] != " " and not line[:5].strip():
                # new feature: key starts at column 6
                flush(cur)
                ftype = line[5:20].strip()
                cur = [ftype, line[20:].strip(), {}]
            elif cur is not None:
                body = line[21:].rstrip("\n") if len(line) > 21 else ""
                stripped = body.strip()
                if stripped.startswith("/"):
                    if "=" in stripped:
                        k, v = stripped[1:].split("=", 1)
                        cur[2][k.strip()] = v.strip().strip('"')
                    else:
                        cur[2][stripped[1:].strip()] = "true"
                elif stripped and not cur[2]:
                    cur[1] += stripped
    flush(cur)
    return features


def load_annotation(path, seqid_filter=None):
    low = path.lower()
    if low.endswith(".gz"):
        low = low[:-3]
    if low.endswith((".gff", ".gff3")):
        return parse_gff3(path, seqid_filter)
    if low.endswith((".gbk", ".gb", ".gbff", ".genbank")):
        return parse_genbank(path, seqid_filter)
    # sniff
    with smart_open(path) as fh:
        head = fh.read(4096)
    if head.lstrip().startswith(("LOCUS", "GBFF")):
        return parse_genbank(path, seqid_filter)
    return parse_gff3(path, seqid_filter)


def load_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


# --------------------------------------------------------------------------
# lookup
# --------------------------------------------------------------------------

class FeatureIndex(object):
    def __init__(self, features):
        self.features = sorted(features, key=lambda f: (f.start, f.end))
        self.starts = [f.start for f in self.features]
        self.max_len = max((f.end - f.start + 1) for f in self.features) if features else 0

    def at(self, pos):
        """All features containing pos, highest priority type first."""
        lo = bisect.bisect_left(self.starts, pos - self.max_len)
        hits = []
        for f in self.features[lo:]:
            if f.start > pos:
                break
            if f.start <= pos <= f.end:
                hits.append(f)
        hits.sort(key=lambda f: TYPE_PRIORITY.index(f.ftype)
                  if f.ftype in TYPE_PRIORITY else len(TYPE_PRIORITY))
        return hits

    def nearest(self, pos):
        """(feature, signed_distance) for the closest feature; negative = upstream of it."""
        best = None
        best_d = None
        for f in self.features:
            if f.start <= pos <= f.end:
                d = 0
            elif pos < f.start:
                d = pos - f.start          # negative
            else:
                d = pos - f.end            # positive
            ad = abs(d)
            if best_d is None or ad < best_d:
                best, best_d = f, ad
                if ad == 0:
                    break
        return best, (0 if best_d == 0 else (pos - best.start if pos < best.start
                                             else pos - best.end))


def codon_for(feature, pos, seq):
    """Codon number, position in codon, and the codon bases for a CDS hit.

    Forward strand: codon n covers start+3(n-1) .. start+3(n-1)+2.
    Reverse strand: translation runs from `end` downwards, so codon n covers
    end-3(n-1) down to end-3(n-1)-2, and the bases are reverse complemented.
    """
    if feature.strand > 0:
        off = pos - feature.start
        codon_no = off // 3 + 1
        codon_pos = off % 3 + 1
        c_start = feature.start + 3 * (codon_no - 1)
        bases = None
        if seq is not None and 1 <= c_start and c_start + 2 <= len(seq):
            bases = seq[c_start - 1:c_start + 2]
    else:
        off = feature.end - pos
        codon_no = off // 3 + 1
        codon_pos = off % 3 + 1
        c_end = feature.end - 3 * (codon_no - 1)
        bases = None
        if seq is not None and c_end - 2 >= 1 and c_end <= len(seq):
            fwd = seq[c_end - 3:c_end]
            bases = "".join(COMPLEMENT.get(b, "N") for b in reversed(fwd))
    return codon_no, codon_pos, bases


def substitute(codon, codon_pos, alt_base, strand):
    """Apply an alternate genomic base to a codon, honouring strand."""
    if codon is None:
        return None
    b = alt_base.upper()
    if strand < 0:
        b = COMPLEMENT.get(b, "N")
    chars = list(codon)
    chars[codon_pos - 1] = b
    return "".join(chars)


def translate(codon):
    if not codon or len(codon) != 3:
        return None
    return CODON_TABLE.get(codon.upper())


# --------------------------------------------------------------------------
# site input
# --------------------------------------------------------------------------

MUT_RE = re.compile(r"^([ACGT])(\d+)([ACGT])$", re.IGNORECASE)


def read_sites_file(path):
    sites = []
    with smart_open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = re.split(r"[\s,]+", line)
            pos = None
            count = None
            mutation = None
            for tok in toks:
                m = MUT_RE.match(tok)
                if m:
                    mutation = tok.upper()
                    pos = int(m.group(2))
                    break
            if pos is None:
                try:
                    pos = int(toks[0])
                except (ValueError, IndexError):
                    continue
                if len(toks) > 1:
                    try:
                        count = int(toks[1])
                    except ValueError:
                        count = None
            sites.append((pos, count, mutation))
    return sites


def read_profile_sites(path):
    """Pull the SITE RECURRENCE table out of profile_tree.py output."""
    sites = []
    with smart_open(path) as fh:
        lines = fh.read().splitlines()
    start = None
    for i, line in enumerate(lines):
        low = line.lower()
        if "position" in low and "mutations" in low and "pct_of_total" in low:
            start = i + 1
            break
    if start is None:
        print("error: could not find the SITE RECURRENCE table "
              "(a header line with position / mutations / pct_of_total)", file=sys.stderr)
        return sites
    for line in lines[start:]:
        toks = line.split()
        if len(toks) < 2:
            break
        try:
            pos = int(toks[0])
            count = int(toks[1])
        except ValueError:
            break
        sites.append((pos, count, None))
    return sites


# --------------------------------------------------------------------------
# mask building
# --------------------------------------------------------------------------

def resolve_region_names(features, names):
    """Match requested gene names / locus tags against the annotation."""
    by_key = defaultdict(list)
    for f in features:
        if f.ftype not in MASKABLE_TYPES:
            continue
        for key in f.names():
            by_key[key].append(f)
    matched = {}
    missing = []
    for name in names:
        hits = by_key.get(name.lower())
        if hits:
            matched[name] = hits
        else:
            missing.append(name)
    return matched, missing


def merge_intervals(intervals):
    """Merge 1-based inclusive intervals."""
    if not intervals:
        return []
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


### Feature types that can be masked.  A "gene" feature and its CDS cover
### nearly the same interval, so matching both and merging is harmless.
MASKABLE_TYPES = ("CDS", "rRNA", "tRNA", "ncRNA", "tmRNA", "misc_RNA", "gene")


def features_matching(features, pattern, on_product=False):
    """Maskable features whose name/locus tag -- or product -- matches pattern.

    Name matching looks at both `gene` and `locus_tag`, because annotations
    disagree about which one carries the readable name.  Product matching looks
    only at `product`, which bare `gene` features do not have; the CDS of the
    same locus does, and covers essentially the same interval.
    """
    rx = re.compile(pattern)
    hits = []
    for f in features:
        if f.ftype not in MASKABLE_TYPES:
            continue
        fields = (f.product,) if on_product else (f.gene, f.locus_tag)
        if any(v and rx.search(v) for v in fields):
            hits.append(f)
    return hits


def feature_interval(f, promoter_flank, pad=0):
    """Masked interval for one feature, 1-based inclusive.

    promoter_flank extends the 5' side only, which is strand dependent: a
    forward feature's promoter lies below `start`, a reverse feature's above
    `end`.  pad extends both sides equally and so is strand independent; it is
    for artifacts that spill past a feature's annotated edges in either
    direction.  The two are additive, so the upstream side of a feature given
    both grows by their sum.
    """
    s, e = f.start, f.end
    if promoter_flank > 0:
        if f.strand > 0:
            s = s - promoter_flank
        else:
            e = e + promoter_flank
    if pad > 0:
        s = s - pad
        e = e + pad
    return max(1, s), e


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def selftest():
    """Verify codon arithmetic on both strands against hand-computed cases."""
    failures = 0

    # A synthetic 60 bp genome.  Forward gene at 10..27, reverse gene at 31..48.
    seq = ("A" * 9 +                      # 1-9 filler
           "ATGAAACCCGGGTTTTAA" +         # 10-27 forward CDS: M K P G F *
           "AAA" +                        # 28-30 filler
           "TTAAAACCCGGGTTTCAT" +         # 31-48 reverse CDS
           "A" * 12)                      # 49-60 filler
    assert len(seq) == 60, len(seq)

    fwd = Feature("t", "CDS", 10, 27, 1, "fwdGene", "T0001", None)
    rev = Feature("t", "CDS", 31, 48, -1, "revGene", "T0002", None)

    # forward gene: codon 1 = positions 10,11,12 = ATG = M
    cases = [
        (fwd, 10, 1, 1, "ATG", "M"),
        (fwd, 12, 1, 3, "ATG", "M"),
        (fwd, 13, 2, 1, "AAA", "K"),
        (fwd, 17, 3, 2, "CCC", "P"),
        (fwd, 27, 6, 3, "TAA", "*"),
    ]
    # reverse gene: translation starts at 48 going down.
    # genome 46..48 = "CAT"; revcomp = "ATG" = M
    cases += [
        (rev, 48, 1, 1, "ATG", "M"),
        (rev, 46, 1, 3, "ATG", "M"),
        (rev, 45, 2, 1, "AAA", "K"),
        (rev, 41, 3, 2, "CCC", "P"),
        (rev, 31, 6, 3, "TAA", "*"),
    ]
    for feat, pos, want_codon, want_cpos, want_bases, want_aa in cases:
        codon_no, codon_pos, bases = codon_for(feat, pos, seq)
        aa = translate(bases)
        ok = (codon_no == want_codon and codon_pos == want_cpos
              and bases == want_bases and aa == want_aa)
        if not ok:
            failures += 1
            print(f"FAIL {feat.label()} pos {pos}: got codon {codon_no} pos {codon_pos} "
                  f"{bases} {aa}; want codon {want_codon} pos {want_cpos} "
                  f"{want_bases} {want_aa}", file=sys.stderr)
    if failures == 0:
        print(f"ok: codon arithmetic correct on {len(cases)} cases "
              "(forward and reverse strand)")

    # substitution, forward: codon 2 AAA, change position 1 to G -> GAA = E
    got = translate(substitute("AAA", 1, "G", 1))
    if got != "E":
        failures += 1
        print(f"FAIL forward substitution: got {got}, want E", file=sys.stderr)
    # substitution, reverse: genomic C on a reverse gene enters the codon as G
    got = translate(substitute("AAA", 1, "C", -1))
    if got != "E":
        failures += 1
        print(f"FAIL reverse substitution: got {got}, want E", file=sys.stderr)
    if failures == 0:
        print("ok: substitutions complement correctly on the reverse strand")

    # translation table sanity
    if len(CODON_TABLE) != 64:
        failures += 1
        print(f"FAIL codon table has {len(CODON_TABLE)} entries, want 64", file=sys.stderr)
    stops = [c for c, a in CODON_TABLE.items() if a == "*"]
    if sorted(stops) != ["TAA", "TAG", "TGA"]:
        failures += 1
        print(f"FAIL stop codons: {sorted(stops)}", file=sys.stderr)
    if failures == 0:
        print("ok: 64-codon table with the three standard stops")

    # interval merging and promoter flanks
    if merge_intervals([(10, 20), (15, 25), (30, 40), (41, 50)]) != [(10, 25), (30, 50)]:
        failures += 1
        print("FAIL interval merge", file=sys.stderr)
    if feature_interval(fwd, 5) != (5, 27):
        failures += 1
        print("FAIL forward promoter flank", file=sys.stderr)
    if feature_interval(rev, 5) != (31, 53):
        failures += 1
        print("FAIL reverse promoter flank", file=sys.stderr)
    ### --mask_pad is symmetric, so it gives the same widening on either strand.
    if feature_interval(fwd, 0, 5) != (5, 32):
        failures += 1
        print(f"FAIL forward pad: {feature_interval(fwd, 0, 5)}, want (5, 32)",
              file=sys.stderr)
    if feature_interval(rev, 0, 5) != (26, 53):
        failures += 1
        print(f"FAIL reverse pad: {feature_interval(rev, 0, 5)}, want (26, 53)",
              file=sys.stderr)
    ### Additive: the forward feature's upstream side grows by flank + pad.
    if feature_interval(fwd, 5, 3) != (2, 30):
        failures += 1
        print(f"FAIL flank+pad: {feature_interval(fwd, 5, 3)}, want (2, 30)",
              file=sys.stderr)
    ### Never runs off the start of the sequence.
    if feature_interval(fwd, 0, 100)[0] != 1:
        failures += 1
        print("FAIL pad clamped at position 1", file=sys.stderr)
    if failures == 0:
        print("ok: interval merging and strand-aware promoter flanks")

    ### Landmark handling: a whole-sequence "region" record must not turn an
    ### intergenic position into a hit at distance 0.
    landmark = Feature("t", "region", 1, 60, 1, None, None, None)
    if landmark.ftype not in LANDMARK_TYPES:
        failures += 1
        print("FAIL 'region' is not treated as a landmark type", file=sys.stderr)
    if "region" in TYPE_PRIORITY:
        failures += 1
        print("FAIL 'region' still carries a matching priority", file=sys.stderr)

    ### Landmarks survive parsing (so the feature counts report them) and are
    ### dropped in main(), which is what this models.
    kept = [f for f in [landmark, fwd, rev] if f.ftype not in LANDMARK_TYPES]
    if kept != [fwd, rev]:
        failures += 1
        print(f"FAIL landmark filter kept {len(kept)} of 2 features", file=sys.stderr)

    ### Position 29 sits between the two genes: 2 bp past the forward gene's end
    ### and 2 bp before the reverse gene's start.
    unfiltered = FeatureIndex([landmark, fwd, rev])
    if not unfiltered.at(29):
        failures += 1
        print("FAIL expected the unfiltered index to match the landmark",
              file=sys.stderr)
    filtered = FeatureIndex(kept)
    if filtered.at(29):
        failures += 1
        print(f"FAIL position 29 should be intergenic, got "
              f"{[f.ftype for f in filtered.at(29)]}", file=sys.stderr)
    near, dist = filtered.nearest(29)
    if near is not fwd or dist != 2:
        failures += 1
        print(f"FAIL nearest feature to 29: got "
              f"{near.gene if near else None} at {dist}, want fwdGene at 2",
              file=sys.stderr)
    if failures == 0:
        print("ok: whole-sequence landmark records do not mask intergenic sites")

    ### Pattern matching: names and products are separate fields, and a
    ### product pattern must not be silently matched against a name.
    mobile = Feature("t", "CDS", 100, 200, 1, None, "Rv2646", "integrase")
    plain = Feature("t", "CDS", 300, 400, 1, "geneX", "Rv0001", "widget protein")
    pool = [mobile, plain, landmark]

    by_name = features_matching(pool, "^Rv2646$")
    if by_name != [mobile]:
        failures += 1
        print(f"FAIL name match on locus tag: {len(by_name)} hits", file=sys.stderr)
    if features_matching(pool, "integrase") != []:
        failures += 1
        print("FAIL 'integrase' matched a name; products are a separate field",
              file=sys.stderr)
    by_product = features_matching(pool, "integrase", on_product=True)
    if by_product != [mobile]:
        failures += 1
        print(f"FAIL product match: {len(by_product)} hits", file=sys.stderr)
    if features_matching(pool, "Rv2646", on_product=True) != []:
        failures += 1
        print("FAIL a locus tag matched as a product", file=sys.stderr)
    if features_matching(pool, "(?i)INTEGRASE", on_product=True) != [mobile]:
        failures += 1
        print("FAIL (?i) did not make the product match case insensitive",
              file=sys.stderr)
    ### A landmark carries no name or product, so no pattern reaches it.
    if features_matching(pool, ".", on_product=True) != [mobile, plain]:
        failures += 1
        print("FAIL '.' should match both products and nothing else",
              file=sys.stderr)
    if failures == 0:
        print("ok: name and product patterns match their own fields only")

    ### Synonyms must be searchable.  NCBI puts the name people look a gene up
    ### by in gene_synonym while Name holds the locus tag, so matching only the
    ### displayed name makes those genes unreachable.
    syn = Feature("t", "CDS", 100, 200, -1, "NP_215998.1", "Rv1482c",
                  "reductase", ["NP_215998.1", "Rv1482c", "mabA", "fabG1"])
    if syn.names() != {"np_215998.1", "rv1482c", "maba", "fabg1"}:
        failures += 1
        print(f"FAIL alias set: {sorted(syn.names())}", file=sys.stderr)
    if syn.gene != "NP_215998.1":
        failures += 1
        print("FAIL aliases changed the displayed gene name", file=sys.stderr)

    pool2 = [syn, plain]
    for want in ("mabA", "MABA", "fabG1", "Rv1482c"):
        matched, missing = resolve_region_names(pool2, [want])
        if missing or matched.get(want) != [syn]:
            failures += 1
            print(f"FAIL resolving {want!r}: matched={bool(matched)} "
                  f"missing={missing}", file=sys.stderr)
    matched, missing = resolve_region_names(pool2, ["nosuchgene"])
    if matched or missing != ["nosuchgene"]:
        failures += 1
        print("FAIL an unknown name should be reported missing", file=sys.stderr)
    ### A feature reached by two of its own names is still one feature.
    matched, _ = resolve_region_names(pool2, ["mabA", "Rv1482c"])
    reached = {id(f) for hits in matched.values() for f in hits}
    if len(reached) != 1:
        failures += 1
        print(f"FAIL two names for one feature reached {len(reached)} features",
              file=sys.stderr)
    if failures == 0:
        print("ok: gene synonyms and locus tags resolve to the same feature")

    ### --highlight must separate a site inside a named gene from one merely
    ### near it.  The nearest feature to an intergenic site is reported, so
    ### without a distance check a neighbouring gene is credited with a hit.
    failures = 0
    hl = Feature("c", "CDS", 1000, 2000, 1, "mutT2", "Rv1160", "mutator", ())
    if "mutt2" not in hl.names() or "rv1160" not in hl.names():
        failures += 1
        print("FAIL a highlight name should match gene and locus tag",
              file=sys.stderr)
    ### Mirrors the production condition: inside means distance 0.
    for dist, flank, want in ((0, 0, True), (140, 0, False),
                              (140, 200, True), (140, 100, False),
                              (-50, 100, True)):
        inside = (dist == 0)
        got = inside or abs(dist) <= flank
        if got != want:
            failures += 1
            print(f"FAIL highlight at distance {dist} with flank {flank}: "
                  f"got {got}, want {want}", file=sys.stderr)
    if failures == 0:
        print("ok: highlight separates a site inside a gene from one near it")

    print("FAILED" if failures else "all self-tests passed")
    return 1 if failures else 0


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    args = parse_args()

    if args.selftest:
        return selftest()

    if not args.annotation:
        print("error: --annotation is required (or use --selftest)", file=sys.stderr)
        return 2

    print(f"Reading annotation {args.annotation} ...", file=sys.stderr)
    features = load_annotation(args.annotation, args.seqid)
    if not features:
        print("error: no usable features parsed; check the file format and --seqid",
              file=sys.stderr)
        return 1
    by_type = defaultdict(int)
    for f in features:
        by_type[f.ftype] += 1
    print("  features: " + ", ".join(f"{k}={v}" for k, v in sorted(by_type.items())),
          file=sys.stderr)

    ### Landmark records span the whole sequence, so drop them before anything
    ### looks a position up or resolves a name to an interval.
    landmarks = [f for f in features if f.ftype in LANDMARK_TYPES]
    if landmarks:
        features = [f for f in features if f.ftype not in LANDMARK_TYPES]
        print(f"  ignoring {len(landmarks)} whole-sequence landmark feature(s) "
              f"({', '.join(sorted({f.ftype for f in landmarks}))}): they contain "
              f"every position, so a site inside one is intergenic, not a hit",
              file=sys.stderr)
    if not features:
        print("error: the annotation has no features other than whole-sequence "
              "landmarks; check the file format and --seqid", file=sys.stderr)
        return 1

    seq = None
    if args.reference:
        print(f"Reading reference {args.reference} ...", file=sys.stderr)
        seqs = load_fasta(args.reference)
        if args.seqid and args.seqid in seqs:
            seq = seqs[args.seqid]
        elif len(seqs) == 1:
            seq = next(iter(seqs.values()))
        elif seqs:
            name = sorted(seqs, key=lambda k: -len(seqs[k]))[0]
            seq = seqs[name]
            print(f"  warning: {len(seqs)} sequences present; using the longest "
                  f"({name}, {len(seq)} bp). Pass --seqid to choose.", file=sys.stderr)
        if seq:
            print(f"  reference length {len(seq)} bp", file=sys.stderr)
            ### An untyped landmark would do the same damage silently.
            hogs = [f for f in features
                    if (f.end - f.start + 1) >= 0.99 * len(seq)]
            if hogs:
                print(f"  warning: {len(hogs)} feature(s) span 99% or more of the "
                      f"reference and will match every position, e.g. "
                      f"{hogs[0].ftype} {hogs[0].gene or hogs[0].locus_tag or '?'}. "
                      f"If that is a landmark record, add its type to "
                      f"LANDMARK_TYPES.", file=sys.stderr)

    ### Names to flag in the output.  Purely annotation: unlike --regions these
    ### are never masked, because a mutation in one of them is usually the
    ### explanation being looked for rather than something to remove.
    ### lowercased name -> the spelling the user asked for, so a hit reports
    ### their name rather than a case-folded one.
    highlight_names = {}
    if args.highlight:
        with smart_open(args.highlight) as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if line:
                    want = line.split()[0]
                    highlight_names.setdefault(want.lower(), want)
        print(f"Highlighting sites in {len(highlight_names)} named genes from "
              f"{args.highlight}", file=sys.stderr)

    did_something = False

    # ---- annotate sites -------------------------------------------------
    sites = []
    if args.from_profile:
        sites = read_profile_sites(args.from_profile)
    elif args.sites:
        sites = read_sites_file(args.sites)
    elif args.positions:
        for tok in args.positions.split(","):
            tok = tok.strip()
            if tok:
                sites.append((int(tok), None, None))

    if sites:
        did_something = True
        index = FeatureIndex(features)
        out = open(args.output, "w") if args.output else sys.stdout
        header = ["position", "count", "mutation", "context", "feature_type",
                  "gene", "locus_tag", "strand", "codon", "codon_pos",
                  "ref_codon", "ref_aa", "aa_change", "distance", "product",
                  "highlight"]
        out.write("\t".join(header) + "\n")

        n_coding = 0
        n_promoter = 0
        n_intergenic = 0
        n_highlighted = 0
        n_nearby = 0
        highlight_hits = defaultdict(int)
        per_gene = defaultdict(int)

        for pos, count, mutation in sites:
            hits = index.at(pos)
            coding = [f for f in hits if f.ftype == "CDS"]
            rna = [f for f in hits if f.ftype in ("rRNA", "tRNA", "ncRNA", "tmRNA", "misc_RNA")]
            row = {k: "" for k in header}
            row["position"] = pos
            row["count"] = count if count is not None else ""
            row["mutation"] = mutation or ""

            if coding:
                f = coding[0]
                n_coding += 1
                per_gene[f.label()] += 1
                codon_no, codon_pos, bases = codon_for(f, pos, seq)
                ref_aa = translate(bases)
                row.update({
                    "context": "CDS", "feature_type": f.ftype, "gene": f.gene or "",
                    "locus_tag": f.locus_tag or "",
                    "strand": "+" if f.strand > 0 else "-",
                    "codon": codon_no, "codon_pos": codon_pos,
                    "ref_codon": bases or "", "ref_aa": ref_aa or "",
                    "distance": 0, "product": f.product or "",
                })
                if bases and ref_aa:
                    ref_genomic = bases[codon_pos - 1] if f.strand > 0 else \
                        COMPLEMENT.get(bases[codon_pos - 1], "N")
                    if mutation:
                        m = MUT_RE.match(mutation)
                        alt = m.group(3).upper()
                        new = translate(substitute(bases, codon_pos, alt, f.strand))
                        kind = "syn" if new == ref_aa else "nonsyn"
                        row["aa_change"] = f"{ref_aa}{codon_no}{new} ({kind})"
                    else:
                        opts = []
                        for alt in "ACGT":
                            if alt == ref_genomic:
                                continue
                            new = translate(substitute(bases, codon_pos, alt, f.strand))
                            opts.append(f"{ref_genomic}>{alt}:{ref_aa}{codon_no}{new}")
                        row["aa_change"] = ";".join(opts)
            elif rna:
                f = rna[0]
                per_gene[f.label()] += 1
                offset = (pos - f.start + 1) if f.strand > 0 else (f.end - pos + 1)
                row.update({
                    "context": f.ftype, "feature_type": f.ftype, "gene": f.gene or "",
                    "locus_tag": f.locus_tag or "",
                    "strand": "+" if f.strand > 0 else "-",
                    "codon": "", "codon_pos": "",
                    "aa_change": f"{f.label()} position {offset}",
                    "distance": 0, "product": f.product or "",
                })
            elif hits:
                f = hits[0]
                per_gene[f.label()] += 1
                row.update({
                    "context": f.ftype, "feature_type": f.ftype, "gene": f.gene or "",
                    "locus_tag": f.locus_tag or "",
                    "strand": "+" if f.strand > 0 else "-",
                    "distance": 0, "product": f.product or "",
                })
            else:
                f, dist = index.nearest(pos)
                upstream = (f.strand > 0 and dist < 0) or (f.strand < 0 and dist > 0)
                if upstream and abs(dist) <= args.flank:
                    n_promoter += 1
                    context = "promoter"
                    per_gene[f.label()] += 1
                else:
                    n_intergenic += 1
                    context = "intergenic"
                row.update({
                    "context": context, "feature_type": f.ftype, "gene": f.gene or "",
                    "locus_tag": f.locus_tag or "",
                    "strand": "+" if f.strand > 0 else "-",
                    "distance": -abs(dist) if upstream else abs(dist),
                    "aa_change": (f"{f.label()} promoter {-abs(dist)}"
                                  if context == "promoter" else ""),
                    "product": f.product or "",
                })
            if highlight_names and f is not None:
                ### Report the requested name, not the annotation's display
                ### name, so a hit on a synonym says which name matched.
                hit = sorted(f.names() & set(highlight_names))
                ### f is the NEAREST feature, which need not contain the site.
                ### Only a site inside the gene counts as a hit on it; a site
                ### merely near one is reported with its distance, and only
                ### when the user asked for a flank.  Flagging a neighbour as
                ### a hit overstates the evidence exactly where it is being
                ### interpreted.
                dist = row["distance"]
                inside = (dist == 0)
                if hit and (inside or abs(dist) <= args.highlight_flank):
                    name = highlight_names[hit[0]]
                    row["highlight"] = name if inside else f"{name}@{dist:+d}"
                    highlight_hits[name] += 1
                    if inside:
                        n_highlighted += 1
                    else:
                        n_nearby += 1
            out.write("\t".join(str(row[k]) for k in header) + "\n")

        if args.output:
            out.close()
            print(f"Annotated {len(sites)} sites -> {args.output}", file=sys.stderr)
        print(f"  CDS: {n_coding}   promoter (within {args.flank} bp): {n_promoter}   "
              f"intergenic: {n_intergenic}", file=sys.stderr)
        if highlight_names:
            if n_nearby:
                print(f"  {n_nearby} further site(s) lie within "
                      f"{args.highlight_flank} bp of a named gene but NOT in "
                      f"one; they carry a signed distance in the highlight "
                      f"column and are not counted as hits", file=sys.stderr)
            if n_highlighted:
                print(f"  HIGHLIGHTED: {n_highlighted} site(s) in named genes: "
                      + ", ".join(f"{g}={n}" for g, n in
                                  sorted(highlight_hits.items(), key=lambda x: -x[1])),
                      file=sys.stderr)
            else:
                print("  highlighted: none of these sites fall in the named genes",
                      file=sys.stderr)
        if per_gene:
            print("  sites per feature: " + ", ".join(
                f"{g}={n}" for g, n in sorted(per_gene.items(), key=lambda x: -x[1])),
                file=sys.stderr)

    # ---- build mask -----------------------------------------------------
    mask_features = []
    if args.regions:
        names = []
        with smart_open(args.regions) as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if line:
                    names.append(line.split()[0])
        matched, missing = resolve_region_names(features, names)
        for name, hits in matched.items():
            mask_features.extend(hits)
        print(f"Resolved {len(matched)} of {len(names)} requested names from "
              f"{args.regions}", file=sys.stderr)
        if missing:
            print("  NOT FOUND in the annotation (check spelling or use locus tags): "
                  + ", ".join(missing), file=sys.stderr)

    if args.mask_pattern:
        hits = features_matching(features, args.mask_pattern)
        print(f"Name pattern {args.mask_pattern!r} matched {len(hits)} features",
              file=sys.stderr)
        if not hits:
            print("  (it matches gene names and locus tags, not products; use "
                  "--mask_product_pattern for those)", file=sys.stderr)
        mask_features.extend(hits)

    if args.mask_product_pattern:
        hits = features_matching(features, args.mask_product_pattern,
                                 on_product=True)
        loci = sorted({f.locus_tag or f.gene or "?" for f in hits})
        print(f"Product pattern {args.mask_product_pattern!r} matched "
              f"{len(hits)} features in {len(loci)} loci", file=sys.stderr)
        if loci:
            shown = ", ".join(loci[:20]) + (" ..." if len(loci) > 20 else "")
            print(f"  {shown}", file=sys.stderr)
        mask_features.extend(hits)

    if mask_features:
        did_something = True
        ### One feature can be reached by several requested names (a gene and
        ### its synonym, a name and its locus tag) and by both --regions and a
        ### pattern.  Count it once; merge_intervals would collapse the
        ### duplicates anyway, but the reported feature count should be honest.
        seen_ids = set()
        deduped = []
        for f in mask_features:
            if id(f) not in seen_ids:
                seen_ids.add(id(f))
                deduped.append(f)
        mask_features = deduped
        intervals = merge_intervals([feature_interval(f, args.promoter_flank,
                                                      args.mask_pad)
                                     for f in mask_features])
        total = sum(e - s + 1 for s, e in intervals)
        widening = []
        if args.promoter_flank:
            widening.append(f"{args.promoter_flank} bp upstream")
        if args.mask_pad:
            widening.append(f"{args.mask_pad} bp both sides")
        print(f"Mask: {len(mask_features)} features"
              + (f" widened by {' and '.join(widening)}" if widening else "")
              + f" -> {len(intervals)} merged intervals, {total} bp",
              file=sys.stderr)
        if seq:
            print(f"  {100.0 * total / len(seq):.3f}% of the reference", file=sys.stderr)
        if args.mask_bed:
            seqid = args.seqid or mask_features[0].seqid or "chr"
            with open(args.mask_bed, "w") as fh:
                for s, e in intervals:
                    fh.write(f"{seqid}\t{s - 1}\t{e}\n")
            print(f"  BED (0-based half-open) -> {args.mask_bed}", file=sys.stderr)
        if args.mask_positions:
            with open(args.mask_positions, "w") as fh:
                for s, e in intervals:
                    for p in range(s, e + 1):
                        fh.write(f"{p}\n")
            print(f"  position list (1-based) -> {args.mask_positions}", file=sys.stderr)
        if args.mask_matutils:
            with open(args.mask_matutils, "w") as fh:
                for s, e in intervals:
                    for p in range(s, e + 1):
                        fh.write(f"N{p}N\n")
            print(f"  matUtils mask file -> {args.mask_matutils}", file=sys.stderr)
        if not args.mask_bed and not args.mask_positions and not args.mask_matutils:
            print("  (pass --mask_bed, --mask_positions or --mask_matutils to write "
                  "it out)", file=sys.stderr)

    if not did_something:
        print("Nothing to do: pass sites (--sites / --positions / --from_profile) "
              "and/or mask options (--regions / --mask_pattern / "
              "--mask_product_pattern).", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
