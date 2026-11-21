#!/usr/bin/env python3
"""Do the splits' spectra show in the resistance mutations they acquire?

Resistance genes were masked before the search, because selection makes their
mutations recur all over the tree.  With the splits found, this script asks
the reverse question: given a split's mutation spectrum, does it acquire
resistance mutations at a different rate, or different ones, than the rest
of the tree, and how much of that does its spectrum account for?

Resistance mutations come from a catalogue of graded variants: the WHO
catalogue of M. tuberculosis complex mutations (3rd edition, Web Annex B, an
.xlsx), with confidence groups chosen by --groups (default 1 and 2: associated
with resistance, and associated - interim).  Only what one single-nucleotide
substitution can make is visible on the tree: missense, stop_gained,
start_lost, stop_lost, upstream (c.-N) and rRNA (n.N) variants, and coding
c.N changes.  Frameshifts, in-frame indels and gene deletions are not, and
the report says how many graded rows that leaves out per drug.  Two rules
reach past the listed variants: a graded GENE_LoF row covers every
stop_gained and start_lost substitution in the gene, and a graded RRDR row
every missense substitution in rpoB codons 426-452.  A variant graded
explicitly (in any group) takes that grade over a rule.

Each mutation on the --input_tree (the UNMASKED tree: the masked one has had
these sites removed) that falls in a catalogue gene or its upstream region
is named against its PARENT's genome, codon by codon, in the catalogue's
style (rpoB_p.Ser450Leu, katG_p.Trp328*, fabG1_c.-15C>T, rrs_n.1401A>G), and
is a resistance event for each drug whose graded variants include that name.
Naming against the parent means a lineage whose background codon differs from
H37Rv gets the name the catalogue would give its isolates.

Splits come from --clade_samples (one line per split and sample, as made in
step 10 from the tree the search ran on), and are placed on this tree at the
most recent common ancestor of their samples; a split whose MRCA's clade has
extra samples is reported.  Branches belong to residual groups, as in
spectra.tsv: a split's group is its clade minus the splits nested in it.

Each group's spectrum is its trinucleotide-context counts on this tree at
unmasked sites (--mask, the masks the search's tree was built with; masked
flanks are not called), shrunk towards the whole tree's spectrum by --prior
pseudo-mutations.  A context's rate in a group is its share of the group's
spectrum over its opportunity (occurrences of the trinucleotide at unmasked
sites of --reference, the root genome).  A resistance substitution's rate is
the rate of its own context, from the reference flanks.

Two questions per split and drug, on all branches and again on internal
branches only (no per-genome error, and no resistance acquired in the patient
sampled):

  Rate (descriptive).  O, the split's resistance events, against two
  expectations that each sum to the drug's events over the tree:
    E_global = events x (split's mutations / tree's mutations), and
    E_spectrum, the same with each group's mutations weighted by the summed
    rate of the drug's resistance substitutions under the group's spectrum.
  spectrum_factor = E_spectrum / E_global is how much more (or less) the
  split's spectrum alone would make resistance mutations available; O /
  E_spectrum is what is left for treatment history and selection to
  explain.  These carry no p value: how often a lineage is treated, and with
  what, swamps the mutational supply, and no null models that.

  Composition.  Among the split's resistance events for the drug, are the
  variants different from those of the rest of the tree?
    G       the 2 x V likelihood-ratio statistic of the split's variant
            counts against the rest's: different, for any reason.
    llr     whether its spectrum predicts the difference.  q_v, the rest's
            share of variant v (with 0.5 added to each), is reweighted by
            rho_v = rate under the split's spectrum / rate under the rest's,
            and renormalised to q'_v.  llr = sum over the split's events of
            log(q'_v / q_v): positive when its events lean towards the
            variants its spectrum favours.
    llr_expected  the llr expected if its events followed q' exactly:
            n * KL(q' || q).  Near zero means the spectrum cannot change
            which variants arise enough to see, whatever the data.
  G is tested by shuffling the drug's events' split labels (the drug's
  variants and each split's event count held fixed).  It is significant
  almost everywhere: which drugs a lineage was treated with, and where,
  differ, so it is a description, not a finding about spectra.
  llr is tested differently, because a composition that differs for those
  reasons and happens to lean the spectrum's way would beat label shuffles
  without the spectrum having anything to do with it (the self-test shows
  shuffling calling a spectrum responsible in about a third of such cases).
  Instead each split's counts are held as observed, and the spectrum ratios
  rho are permuted across the drug's variants (one permutation for all
  splits): does the composition lean towards the variants this spectrum
  favours more than towards a random reassignment of the same ratios?
  Each test is run 'permutations' times; p = (1 + permutations at or above
  the observed) / (1 + permutations), and q is Benjamini-Hochberg over the
  splits tested for that drug.  A split is tested with --min_events or more
  events.  Per drug, the sum of llr over the tested splits is tested on the
  same ratio permutations, and G over the whole splits x variants table on
  the same label shuffles.  A drug with few variants (amikacin's handful)
  has few distinct permutations, and its llr p values cannot get small.

'any drug' pools the drugs: each distinct event (branch, site, base) once,
over every drug's variants.

Outputs, named from --prefix (default resistance_spectra):
  PREFIX_events.tsv   one row per resistance event and drug
  PREFIX.tsv          one row per reading, drug and split
  PREFIX_drugs.tsv    one row per reading and drug: totals and the pooled
                      tests
  PREFIX_variants.tsv per reading, drug and tested split, the shares
                      behind llr: each variant's share of the rest (q), as
                      the split's spectrum reweights it (q'), and of the
                      split's own events; for plot_resistance_spectra.py
The report goes to stdout.

Usage:
    resistance_spectra.py --input_tree tree.pruned.pb.gz \\
        --clade_samples split_clade_samples.tsv \\
        --catalogue WHO-catalogue-3rd-edition-annex-B.xlsx \\
        --annotation NC_000962.3.gff3 --reference root_genome.fa \\
        --mask selection_loci.bed,alignment_loci.bed --fold \\
        --max_branch_length 250 > resistance_spectra.txt

    resistance_spectra.py --selftest

Needs numpy, bte, and openpyxl to read the catalogue.
"""

import argparse
import math
import os
import re
import sys
from array import array
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trinucleotide as tn
from annotate_sites import load_annotation, translate
from split_contexts import bh, fmt_p, read_fasta

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
       "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
       "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
       "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
       "*": "Ter"}
### Effects one substitution can produce, as the catalogue names them.
SNV_EFFECTS = {"missense_variant", "stop_gained", "start_lost", "stop_lost",
               "upstream_gene_variant", "non_coding_transcript_exon_variant",
               "synonymous_variant", "initiator_codon_variant"}
LOF_EFFECTS = {"stop_gained", "start_lost"}
RRDR_GENE, RRDR_CODONS = "rpoB", (426, 452)
ANY = "any drug"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Test each split's resistance mutations against its "
                    "mutation spectrum.")
    p.add_argument("--input_tree",
                   help="The UNMASKED tree (protobuf): resistance sites must "
                        "still carry their mutations")
    p.add_argument("--clade_samples",
                   help="split<TAB>sample per line, every sample of every "
                        "split's clade (step 10's split_clade_samples.tsv)")
    p.add_argument("--catalogue",
                   help="WHO catalogue Web Annex B (.xlsx), with the sheets "
                        "Graded_Mutation_List and Gene_Coordinates")
    p.add_argument("--groups", default="1,2",
                   help="Confidence groups counted as resistance "
                        "(default 1,2)")
    p.add_argument("--annotation",
                   help="GFF3 or GenBank of the reference, for each "
                        "catalogue gene's strand and type")
    p.add_argument("--reference",
                   help="FASTA of the tree's ROOT genome")
    p.add_argument("--mask", default=None,
                   help="Comma-separated BED files: the masks the search's "
                        "tree was built with.  Spectra and opportunity are "
                        "counted outside them")
    p.add_argument("--fold", action="store_true",
                   help="96 folded context classes, as for a --fold search; "
                        "without it, 192 forward-strand ones")
    p.add_argument("--max_branch_length", type=int, default=None,
                   help="Branches with more mutations than this count for "
                        "nothing (match the search)")
    p.add_argument("--prior", type=float, default=100.0,
                   help="Pseudo-mutations shrinking each group's spectrum "
                        "towards the whole tree's (default 100)")
    p.add_argument("--min_events", type=int, default=5,
                   help="Events a split needs for its composition to be "
                        "tested (default 5)")
    p.add_argument("--permutations", type=int, default=10000,
                   help="Label shuffles per drug and reading (default 10000)")
    p.add_argument("--max_q", type=float, default=0.05,
                   help="q below which a split is listed (default 0.05)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--prefix", default="resistance_spectra",
                   help="Prefix of the output TSVs (default "
                        "resistance_spectra)")
    p.add_argument("--selftest", action="store_true",
                   help="Run the self-tests and exit")
    return p.parse_args(argv)


# --------------------------------------------------------------------------
# catalogue
# --------------------------------------------------------------------------

class Gene(object):
    __slots__ = ("name", "locus", "start", "end", "strand", "rna")

    def __init__(self, name, locus, start, end, strand, rna):
        self.name, self.locus = name, locus
        self.start, self.end = start, end
        self.strand, self.rna = strand, rna

    def coding_index(self, pos):
        """1-based index from the gene's first base, on its strand."""
        return (pos - self.start + 1 if self.strand > 0
                else self.end - pos + 1)

    def genomic(self, k):
        return self.start + k - 1 if self.strand > 0 else self.end - k + 1

    def upstream(self, pos):
        """Distance upstream of the gene's first base, or 0."""
        d = self.start - pos if self.strand > 0 else pos - self.end
        return d if d > 0 else 0


class Catalogue(object):
    """Graded variants, rules and gene coordinates.

    explicit[(drug, gene, name)] = (group, route)   every graded row
    lof[(drug, gene)] = (group, route)              GENE_LoF rows
    rrdr[drug] = (group, route)                     RRDR rows
    """

    def __init__(self, rows, coords, groups):
        self.groups = set(groups)
        self.explicit = {}
        self.lof = {}
        self.rrdr = {}
        self.coords = coords                 # gene -> (locus, left, right)
        self.listed = defaultdict(Counter)   # drug -> effect -> rows in groups
        self.listed_names = defaultdict(set)  # drug -> {(gene, name)} SNV-type
        self.max_up = Counter()              # gene -> furthest listed c.-N
        for r in rows:
            drug, gene, var = r["drug"], r["gene"], r["variant"]
            grade, route, effect = r["grade"], r["route"], r["effect"]
            if grade is None:
                continue
            if var == "RRDR" or gene == "RRDR":
                self.rrdr[drug] = (grade, route or "Rule: RRDR")
                if grade in self.groups:
                    self.listed[drug]["RRDR rule"] += 1
                continue
            name = var[len(gene) + 1:] if var.startswith(gene + "_") \
                else var.split("_", 1)[-1]
            if name == "LoF":
                self.lof[(drug, gene)] = (grade, route or "Rule: LoF")
                if grade in self.groups:
                    self.listed[drug]["LoF rule"] += 1
                continue
            self.explicit[(drug, gene, name)] = (grade, route)
            if grade not in self.groups:
                continue
            self.listed[drug][effect] += 1
            if effect in SNV_EFFECTS:
                self.listed_names[drug].add((gene, name))
            m = re.match(r"c\.-(\d+)", name)
            if m:
                self.max_up[gene] = max(self.max_up[gene], int(m.group(1)))

    def drugs(self):
        out = set(d for d, _g, _n in self.explicit
                  if self.explicit[(d, _g, _n)][0] in self.groups)
        out |= set(d for d, _g in self.lof if self.lof[(d, _g)][0]
                   in self.groups)
        out |= set(d for d in self.rrdr if self.rrdr[d][0] in self.groups)
        return sorted(out)

    def genes_of(self, drug):
        """Genes with something graded in --groups for the drug."""
        out = set(g for d, g, _n in self.explicit
                  if d == drug and self.explicit[(d, g, _n)][0] in self.groups)
        out |= set(g for d, g in self.lof
                   if d == drug and self.lof[(d, g)][0] in self.groups)
        if drug in self.rrdr and self.rrdr[drug][0] in self.groups:
            out.add(RRDR_GENE)
        return out

    def classify(self, drug, gene, name, effect, codon):
        """(group, route) if this named substitution is resistance for the
        drug in --groups, else None.  An explicit grade beats a rule."""
        hit = self.explicit.get((drug, gene, name))
        if hit is not None:
            grade, route = hit
            if grade not in self.groups:
                return None
            return grade, route or "listed"
        if effect in LOF_EFFECTS:
            rule = self.lof.get((drug, gene))
            if rule and rule[0] in self.groups:
                return rule[0], "Rule: LoF (not listed)"
        if effect == "missense_variant" and gene == RRDR_GENE \
                and RRDR_CODONS[0] <= codon <= RRDR_CODONS[1]:
            rule = self.rrdr.get(drug)
            if rule and rule[0] in self.groups:
                return rule[0], "Rule: RRDR (not listed)"
        return None


def _grade(value):
    m = re.match(r"\s*(\d)\)", str(value or ""))
    return int(m.group(1)) if m else None


def read_catalogue(path, groups):
    try:
        import openpyxl
    except ImportError:
        raise SystemExit("error: reading the catalogue needs openpyxl")
    import warnings
    warnings.simplefilter("ignore")
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)

    def table(sheet, first):
        it = wb[sheet].iter_rows(values_only=True)
        for r in it:
            if r and r[0] == first:
                hdr = [str(x) if x is not None else "" for x in r]
                break
        else:
            raise SystemExit(f"error: no header row in {sheet}")
        for r in it:
            if r and r[0]:
                yield dict(zip(hdr, r))

    rows = []
    for r in table("Graded_Mutation_List", "drug"):
        rows.append({
            "drug": str(r["drug"]), "gene": str(r["Gene"]),
            "variant": str(r["variant"]), "effect": str(r["effect"]),
            "grade": _grade(r.get("FINAL CONFIDENCE GRADING")),
            "route": r.get("SUPPLEMENTARY GRADING CONSIDERATIONS")})
    coords = {}
    for r in table("Gene_Coordinates", "Gene"):
        try:
            coords[str(r["Gene"])] = (str(r["Rv ID"]), int(r["Left"]),
                                      int(r["Right"]))
        except (TypeError, ValueError):
            continue
    return Catalogue(rows, coords, groups)


def build_genes(cat, features, names, warn):
    """Gene objects for the named catalogue genes, strand and type from the
    annotation's feature with the catalogue's locus tag."""
    by_locus = defaultdict(list)
    for f in features:
        if f.locus_tag:
            by_locus[f.locus_tag.lower()].append(f)
    out = {}
    for g in sorted(names):
        if g not in cat.coords:
            warn(f"{g}: no coordinates in the catalogue; left out")
            continue
        locus, left, right = cat.coords[g]
        fs = by_locus.get(locus.lower(), [])
        if not fs:
            warn(f"{g} ({locus}): not in the annotation; left out")
            continue
        rna = any(f.ftype == "rRNA" for f in fs)
        coding = [f for f in fs if f.ftype in ("CDS", "rRNA")] or fs
        if not any(f.start == left and f.end == right for f in coding):
            warn(f"{g} ({locus}): catalogue {left}-{right}, annotation "
                 f"{coding[0].start}-{coding[0].end}; using the catalogue's")
        out[g] = Gene(g, locus, left, right, coding[0].strand, rna)
    return out


# --------------------------------------------------------------------------
# naming
# --------------------------------------------------------------------------

def snv_names(gene, pos, base_of, alt):
    """[(name, effect, codon)] for a substitution to `alt` at `pos`, named
    against the genome base_of(p) gives, in the catalogue's style."""
    comp = (lambda b: b) if gene.strand > 0 else \
        (lambda b: COMPLEMENT.get(b, "N"))
    if gene.start <= pos <= gene.end:
        k = gene.coding_index(pos)
        r, a = comp(base_of(pos)), comp(alt)
        if gene.rna:
            return [(f"n.{k}{r}>{a}", "non_coding_transcript_exon_variant",
                     0)]
        out = [(f"c.{k}{r}>{a}", "nucleotide", 0)]
        n, cp = (k - 1) // 3 + 1, (k - 1) % 3
        first = 3 * (n - 1) + 1
        if gene.coding_index(gene.end if gene.strand > 0 else gene.start) \
                < first + 2:
            return out
        codon = "".join(comp(base_of(gene.genomic(first + j)))
                        for j in range(3))
        mutant = codon[:cp] + a + codon[cp + 1:]
        ra, aa = translate(codon), translate(mutant)
        if ra is None or aa is None:
            return out
        if n == 1:
            if mutant != "ATG":
                out.append(("p.Met1?", "start_lost", 1))
        elif ra == aa:
            out.append((f"p.{AA3[ra]}{n}=", "synonymous_variant", n))
        elif aa == "*":
            out.append((f"p.{AA3[ra]}{n}*", "stop_gained", n))
        elif ra == "*":
            out.append((f"p.Ter{n}{AA3[aa]}ext*?", "stop_lost", n))
        else:
            out.append((f"p.{AA3[ra]}{n}{AA3[aa]}", "missense_variant", n))
        return out
    d = gene.upstream(pos)
    if d:
        return [(f"c.-{d}{comp(base_of(pos))}>{comp(alt)}",
                 "upstream_gene_variant", 0)]
    return []


def restores_reference(gene, pos, base_of, alt, ref_of):
    """True if the substitution gives back the reference's amino acid (in a
    coding gene's codon) or base (anywhere else): a reversion towards
    susceptibility, whatever name it gets against the parent."""
    if gene.rna or not gene.start <= pos <= gene.end:
        return alt == ref_of(pos)
    comp = (lambda b: b) if gene.strand > 0 else \
        (lambda b: COMPLEMENT.get(b, "N"))
    k = gene.coding_index(pos)
    first = 3 * ((k - 1) // 3) + 1
    cod = [gene.genomic(first + j) for j in range(3)]
    child = "".join(comp(alt if g == pos else base_of(g)) for g in cod)
    refc = "".join(comp(ref_of(g)) for g in cod)
    if first == 1:
        return child == refc
    ca, ra = translate(child), translate(refc)
    return ca is not None and ca == ra


def window(gene, max_up):
    if gene.strand > 0:
        return gene.start - max_up, gene.end
    return gene.start, gene.end + max_up


class Targets(object):
    """Which drugs' resistance each substitution is, by position."""

    def __init__(self, cat, genes):
        self.cat = cat
        self.genes = genes
        self.drugs = cat.drugs()
        self.drugs_of = defaultdict(list)    # gene -> drugs
        for d in self.drugs:
            for g in cat.genes_of(d):
                if g in genes:
                    self.drugs_of[g].append(d)
        self.at = defaultdict(list)          # pos -> genes
        for g, gene in genes.items():
            if g not in self.drugs_of:
                continue
            lo, hi = window(gene, cat.max_up[g])
            for p in range(lo, hi + 1):
                self.at[p].append(gene)

    def hits(self, pos, base_of, alt):
        """[(drug, gene, name, effect, group, route)], one per drug."""
        out, seen = [], set()
        for gene in self.at.get(pos, ()):
            for name, effect, codon in snv_names(gene, pos, base_of, alt):
                for d in self.drugs_of[gene.name]:
                    if d in seen:
                        continue
                    c = self.cat.classify(d, gene.name, name, effect, codon)
                    if c:
                        seen.add(d)
                        out.append((d, gene.name, name, effect) + c)
        return out

    def alternatives(self, ref):
        """drug -> {(pos, alt): (gene, name)} reachable from the reference,
        plus the listed SNV-type names no single substitution reaches."""
        alts = defaultdict(dict)
        reached = defaultdict(set)
        base_of = lambda p: ref[p - 1] if 1 <= p <= len(ref) else "N"
        for pos in sorted(self.at):
            r = base_of(pos)
            for a in "ACGT":
                if a == r or r not in COMPLEMENT:
                    continue
                for d, g, name, _e, _gr, _ro in self.hits(pos, base_of, a):
                    alts[d].setdefault((pos, a), (g, name))
                    reached[d].add((g, name))
        unreached = {d: sorted(self.cat.listed_names[d] - reached[d])
                     for d in self.drugs}
        return alts, unreached


# --------------------------------------------------------------------------
# tree
# --------------------------------------------------------------------------

def find_events(tree, targets, ref, forward=None):
    """[(node, pos, parent_base, alt, drug, gene, name, effect, group,
    route)] for every resistance substitution on the tree, named against
    the parent's genome (reconstructed down the tree), and per drug the
    back-mutations of its resistance substitutions: a return to the
    reference base where the parent carries a substitution that `forward`
    ({(pos, base): drugs}) lists.  Substitutions that give back the
    reference's amino acid or base are never events, whatever name they get
    against the parent (rpoB Leu450Ser on an S450L background is in the
    RRDR, and no resistance)."""
    reversions = Counter()
    ref_of = lambda p: ref[p - 1] if 1 <= p <= len(ref) else "N"
    n = len(tree.node_ids)
    L = len(ref)
    watched = targets.at
    state, undo = {}, {}
    events = []
    roots = [i for i in range(n) if tree.parent[i] < 0]
    stack = [(r, False) for r in reversed(roots)]

    def base_of(p):
        b = state.get(p)
        return b if b is not None else (ref[p - 1] if 1 <= p <= L else "N")

    while stack:
        i, leaving = stack.pop()
        if leaving:
            for p, old in reversed(undo.pop(i)):
                if old is None:
                    del state[p]
                else:
                    state[p] = old
            continue
        changes = []
        mine = [(p, tn.MUT_TYPES[t][1]) for p, t in
                zip(tree.positions[i], tree.types[i]) if p in watched]
        for p, alt in mine:
            pb = base_of(p)
            if pb == alt:
                continue
            if forward and alt == ref_of(p):
                for d in forward.get((p, pb), ()):
                    reversions[d] += 1
            for h in targets.hits(p, base_of, alt):
                if not restores_reference(targets.genes[h[1]], p, base_of,
                                          alt, ref_of):
                    events.append((i, p, pb, alt) + h)
        for p, alt in mine:
            changes.append((p, state.get(p)))
            state[p] = alt
        undo[i] = changes
        stack.append((i, True))
        for c in reversed(tree.children[i]):
            stack.append((c, False))
    return events, reversions


def subtree_sizes(tree):
    n = len(tree.node_ids)
    size = array("i", [1] * n)
    for i in range(n - 1, 0, -1):
        p = tree.parent[i]
        if p >= 0:
            size[p] += size[i]
    return size


def place_splits(tree, clade_samples):
    """split -> (node index, extra samples in its MRCA's clade, samples not
    on this tree).  Relies on extract()'s preorder: a node's subtree is the
    index range [i, i + size)."""
    size = subtree_sizes(tree)
    leaf = {tree.node_ids[i]: i for i in range(len(tree.node_ids))
            if not tree.children[i]}
    leaves_below = array("i", [0] * len(tree.node_ids))
    for i in range(len(tree.node_ids) - 1, -1, -1):
        if not tree.children[i]:
            leaves_below[i] += 1
        p = tree.parent[i]
        if p >= 0:
            leaves_below[p] += leaves_below[i]
    out = {}
    for split, samples in clade_samples.items():
        idx = [leaf[s] for s in samples if s in leaf]
        missing = len(samples) - len(idx)
        if not idx:
            out[split] = (None, 0, missing)
            continue
        lo, hi = min(idx), max(idx)
        m = lo
        while not (m <= hi < m + size[m]):
            m = tree.parent[m]
        out[split] = (m, leaves_below[m] - len(idx), missing)
    return out


def read_clade_samples(path):
    out = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and f[1]:
                out[f[0]].add(f[1])
    return out


def groups_of(tree, roots):
    """Nearest ancestor-or-self split root per node, -1 outside every split."""
    group = array("i", [-1] * len(tree.node_ids))
    for i in range(len(tree.node_ids)):
        if i in roots:
            group[i] = i
        else:
            p = tree.parent[i]
            group[i] = group[p] if p >= 0 else -1
    return group


def counted(tree, i, cap, internal_only):
    if cap is not None and tree.raw_len[i] > cap:
        return False
    return not internal_only or bool(tree.children[i])


def group_spectra(np, tree, contexts, masked, group, cap, internal_only,
                  n_ctx):
    """group -> context counts at unmasked sites."""
    out = defaultdict(lambda: np.zeros(n_ctx))
    for i in range(len(tree.node_ids)):
        g = group[i]
        if not counted(tree, i, cap, internal_only):
            continue
        c = [x for p, x in zip(tree.positions[i], contexts[i])
             if x >= 0 and not masked[p]]
        if c:
            out[g] += np.bincount(c, minlength=n_ctx)
    return out


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------

def rates(np, counts, total, opp_ctx, prior):
    """Per-context rate (share over opportunity), the share shrunk towards
    `total`'s by `prior` pseudo-mutations."""
    base = (total + 0.5) / (total + 0.5).sum()
    share = (counts + prior * base) / (counts.sum() + prior)
    return share / opp_ctx


def snv_context(sch, ref, pos, rbase, alt):
    L = len(ref)
    left = ref[pos - 2] if pos >= 2 else "N"
    right = ref[pos] if pos < L else "N"
    f = tn.context_index if sch.fold else tn.stranded_context_index
    return f(left, rbase, right, alt)


def _xlogx(np, a):
    a = np.asarray(a, dtype=np.float64)
    out = np.zeros_like(a)
    pos = a > 0
    out[pos] = a[pos] * np.log(a[pos])
    return out


def composition_stats(np, C, logrho, tested):
    """Per row of C (groups x variants): llr, G against the other rows, and
    the drug's total llr and whole-table G.  logrho is groups x variants."""
    tot = C.sum(0)
    N = tot.sum()
    n = C.sum(1)
    V = C.shape[1]
    rest = tot[None, :] - C
    q = (rest + 0.5) / (rest.sum(1, keepdims=True) + 0.5 * V)
    w = q * np.exp(logrho)
    qp = w / w.sum(1, keepdims=True)
    llr = (C * (np.log(qp) - np.log(q))).sum(1)
    xl = _xlogx(np, C).sum(1) + _xlogx(np, rest).sum(1)
    G = 2 * (xl - _xlogx(np, n) - _xlogx(np, N - n)
             - _xlogx(np, tot).sum() + _xlogx(np, N))
    Gall = 2 * (_xlogx(np, C).sum() - _xlogx(np, n).sum()
                - _xlogx(np, tot).sum() + _xlogx(np, N))
    expected = n * (qp * (np.log(qp) - np.log(q))).sum(1)
    return llr, G, float(llr[tested].sum()), float(Gall), expected


def shuffle_test(np, labels, cats, S, V, logrho, tested, permutations, rng):
    """Observed statistics, and the G tests' p values from shuffling the
    events' split labels."""
    def table(lab):
        return np.bincount(lab * V + cats, minlength=S * V).reshape(S, V) \
            .astype(np.float64)
    obs = composition_stats(np, table(labels), logrho, tested)
    ge_G = np.zeros(S)
    ge_all = 0
    for _ in range(permutations):
        _l, G, _t, Gall, _e = composition_stats(np, table(rng.permutation(
            labels)), logrho, tested)
        ge_G += G >= obs[1] - 1e-9
        ge_all += Gall >= obs[3] - 1e-9
    P = permutations + 1
    return obs, (ge_G + 1) / P, (ge_all + 1) / P


def rho_test(np, C, logrho, tested, permutations, rng, chunk=32):
    """p values for llr, per tested row and summed over them, from permuting
    the spectrum ratios across the variants: the same permutation for every
    split.  Each split's variant counts stay as observed, so a composition
    that differs from the rest's for any other reason (treatment, lineage)
    is held fixed, and only its alignment with the spectrum is tested."""
    Ct = C[tested]
    L = logrho[tested]
    n = Ct.sum(1)
    V = C.shape[1]
    tot = C.sum(0)
    rest = tot[None, :] - Ct
    q = (rest + 0.5) / (rest.sum(1, keepdims=True) + 0.5 * V)

    def llr_of(Lk):                          # (S, k, V) -> (S, k)
        a = np.einsum("sv,skv->sk", Ct, Lk)
        b = np.log(np.einsum("sv,skv->sk", q, np.exp(Lk)))
        return a - n[:, None] * b

    obs = llr_of(L[:, None, :])[:, 0]
    ge = np.zeros(len(n))
    ge_tot = 0
    done = 0
    while done < permutations:
        k = min(chunk, permutations - done)
        idx = np.stack([rng.permutation(V) for _ in range(k)])
        null = llr_of(L[:, idx])
        ge += (null >= obs[:, None] - 1e-9).sum(1)
        ge_tot += int((null.sum(0) >= obs.sum() - 1e-9).sum())
        done += k
    P = permutations + 1
    p = np.ones(C.shape[0])
    p[tested] = (ge + 1) / P
    return p, (ge_tot + 1) / P


def analyse_drug(np, drug, events, alts, groups, gid, exposure, rate_of,
                 sch, ref, args, rng):
    """Rows per split, and the drug's totals, for one drug and reading.

    events: [(group index, pos, parent base, alt, name)]; alts: every
    resistance substitution reachable from the reference, for the rate
    expectation; gid: group -> row; exposure: per row, its mutations;
    rate_of(row, ctx) -> rate."""
    S = len(groups)
    N = len(events)
    # rate expectations
    actx = [snv_context(sch, ref, p, ref[p - 1], a) for p, a in alts]
    actx = np.array([c for c in actx if c >= 0], dtype=int)
    R = np.array([rate_of(s, actx).sum() for s in range(S)])
    T = exposure
    e_glob = N * T / T.sum() if T.sum() else np.zeros(S)
    tr = T * R
    e_spec = N * tr / tr.sum() if tr.sum() else np.zeros(S)
    obs = np.zeros(S)
    for g, *_r in events:
        obs[g] += 1
    # composition
    cat_of, names = {}, []
    cats, labels, cctx = [], [], []
    for g, p, pb, a, name in events:
        key = (p, a)
        if key not in cat_of:
            cat_of[key] = len(names)
            names.append(name)
            cctx.append(snv_context(sch, ref, p, pb, a))
        cats.append(cat_of[key])
        labels.append(g)
    V = len(names)
    rows = []
    tested = obs >= args.min_events
    res = None
    if V >= 2 and tested.sum() >= 1:
        cctx = np.array(cctx)
        ok = cctx >= 0
        allr = np.stack([rate_of(s, np.where(ok, cctx, 0)) for s in range(S)])
        allr[:, ~ok] = 1.0
        sumT = exposure[:, None] * allr
        logrho = np.zeros((S, V))
        for s in range(S):
            rest = (sumT.sum(0) - sumT[s]) / max(T.sum() - T[s], 1e-12)
            logrho[s] = np.log(allr[s]) - np.log(rest)
        res = shuffle_test(np, np.array(labels), np.array(cats), S, V,
                           logrho, tested, args.permutations, rng)
        C = np.bincount(np.array(labels) * V + np.array(cats),
                        minlength=S * V).reshape(S, V).astype(np.float64)
        rho_p, rho_tot = rho_test(np, C, logrho, tested, args.permutations,
                                  rng)
    llr_q = bh([float(rho_p[s]) if res is not None and tested[s] else None
                for s in range(S)])
    g_q = bh([float(res[1][s]) if res is not None and tested[s] else None
              for s in range(S)])
    for s in range(S):
        row = {"drug": drug, "split": groups[s], "events": int(obs[s]),
               "mutations": float(T[s]), "e_global": float(e_glob[s]),
               "e_spectrum": float(e_spec[s]),
               "spectrum_factor": (float(e_spec[s] / e_glob[s])
                                   if e_glob[s] else None),
               "obs_over_spectrum": (float(obs[s] / e_spec[s])
                                     if e_spec[s] else None),
               "llr": None, "llr_expected": None, "llr_p": None,
               "llr_q": llr_q[s], "g": None, "g_p": None, "g_q": g_q[s],
               "top_variant": ""}
        if res is not None and tested[s]:
            stats = res[0]
            row.update(llr=float(stats[0][s]), g=float(stats[1][s]),
                       llr_expected=float(stats[4][s]),
                       llr_p=float(rho_p[s]), g_p=float(res[1][s]))
            row["top_variant"] = top_variant(np, s, labels, cats, names, V)
        rows.append(row)
    drug_row = {"drug": drug, "events": N, "variants_seen": V,
                "substitutions": len(alts), "splits_tested": int(tested.sum()),
                "llr_total": None, "llr_total_p": None, "g_all": None,
                "g_all_p": None}
    if res is not None:
        drug_row.update(llr_total=res[0][2], llr_total_p=float(rho_tot),
                        g_all=res[0][3], g_all_p=float(res[2]))
    var_rows = variant_rows(np, drug, groups, C, logrho, tested, names) \
        if res is not None else []
    return rows, drug_row, var_rows


def variant_rows(np, drug, groups, C, logrho, tested, names, min_share=0.005):
    """Per tested split, the shares behind its llr: q (the rest's, with 0.5
    added to each variant), q' (q reweighted by the split's spectrum) and
    the split's own, for each variant the split has or that makes up
    min_share or more of q or q'."""
    out = []
    tot = C.sum(0)
    V = C.shape[1]
    for s in np.flatnonzero(tested):
        rest = tot - C[s]
        q = (rest + 0.5) / (rest.sum() + 0.5 * V)
        w = q * np.exp(logrho[s])
        qp = w / w.sum()
        n = C[s].sum()
        for v in range(V):
            if C[s, v] or q[v] >= min_share or qp[v] >= min_share:
                out.append({"drug": drug, "split": groups[s],
                            "variant": names[v], "observed": int(C[s, v]),
                            "split_events": int(n),
                            "observed_share": float(C[s, v] / n),
                            "rest_share": float(q[v]),
                            "spectrum_share": float(qp[v])})
    return out


def top_variant(np, s, labels, cats, names, V):
    """The variant most over-represented in group s against the rest:
    name, observed, expected from the rest's shares."""
    lab = np.array(labels)
    cat = np.array(cats)
    mine = np.bincount(cat[lab == s], minlength=V)
    rest = np.bincount(cat[lab != s], minlength=V)
    n = mine.sum()
    exp = n * (rest + 0.5) / (rest.sum() + 0.5 * V)
    v = int(np.argmax(mine - exp))
    return f"{names[v]} {mine[v]}/{exp[v]:.1f}"


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def fmt(v, spec=".3g"):
    return "-" if v is None else format(v, spec)


def report_catalogue(cat, alts, unreached, out):
    print(f"Catalogue, groups {','.join(map(str, sorted(cat.groups)))}: "
          f"graded rows per drug, and what one substitution can reach",
          file=out)
    print(f"  {'drug':<14} {'rows':>5} {'SNV-type':>8} {'unreached':>9} "
          f"{'rules':<14} {'substitutions':>13}  not visible on the tree",
          file=out)
    for d in cat.drugs():
        eff = cat.listed[d]
        rows = sum(eff.values())
        snv = sum(v for k, v in eff.items() if k in SNV_EFFECTS)
        rules = [k.split()[0] for k in eff if k.endswith("rule")]
        other = Counter({k: v for k, v in eff.items()
                         if k not in SNV_EFFECTS and not k.endswith("rule")})
        print(f"  {d:<14} {rows:>5} {snv:>8} {len(unreached[d]):>9} "
              f"{','.join(rules) or '-':<14} {len(alts.get(d, ())):>13}  "
              + ", ".join(f"{k} {v}" for k, v in other.most_common()),
              file=out)
    print("  rows: graded rows (a rule row counts once).  SNV-type: listed "
          "variants one substitution makes.\n  unreached: SNV-type names "
          "no single substitution from the reference makes.  substitutions:"
          "\n  the (site, base) pairs, listed or ruled, that make a "
          "resistance variant.", file=out)


def report_reading(label, rows, drug_rows, args, out):
    print(f"\n=== {label} ===", file=out)
    print(f"  {'drug':<14} {'events':>6} {'variants':>8} {'tested':>6} "
          f"{'llr_total':>9} {'p':>7} {'G_all':>8} {'p':>7}", file=out)
    for d in drug_rows:
        print(f"  {d['drug']:<14} {d['events']:>6} {d['variants_seen']:>8} "
              f"{d['splits_tested']:>6} {fmt(d['llr_total'], '.2f'):>9} "
              f"{fmt_p(d['llr_total_p']):>7} {fmt(d['g_all'], '.1f'):>8} "
              f"{fmt_p(d['g_all_p']):>7}", file=out)
    anyrows = [r for r in rows if r["drug"] == ANY]
    print(f"\n  {ANY}, per split (composition tested at "
          f">= {args.min_events} events):", file=out)
    print(f"  {'split':<12} {'events':>6} {'E_glob':>7} {'factor':>6} "
          f"{'O/E_sp':>6} {'llr':>7} {'llr_exp':>7} {'p':>7} {'q':>7} "
          f"{'G':>7} {'p':>7} {'q':>7}  top variant obs/exp", file=out)
    for r in sorted(anyrows, key=lambda r: -r["events"]):
        print(f"  {r['split']:<12} {r['events']:>6} "
              f"{fmt(r['e_global'], '.1f'):>7} "
              f"{fmt(r['spectrum_factor'], '.2f'):>6} "
              f"{fmt(r['obs_over_spectrum'], '.2f'):>6} "
              f"{fmt(r['llr'], '.2f'):>7} {fmt(r['llr_expected'], '.2f'):>7} "
              f"{fmt_p(r['llr_p']):>7} {fmt_p(r['llr_q']):>7} "
              f"{fmt(r['g'], '.1f'):>7} {fmt_p(r['g_p']):>7} "
              f"{fmt_p(r['g_q']):>7}  {r['top_variant']}", file=out)
    hits = [r for r in rows if r["drug"] != ANY and (
        (r["llr_q"] is not None and r["llr_q"] < args.max_q)
        or (r["g_q"] is not None and r["g_q"] < args.max_q))]
    print(f"\n  Per drug, splits with llr or G q < {args.max_q}: "
          f"{len(hits) or 'none'}", file=out)
    for r in sorted(hits, key=lambda r: (r["drug"], r["llr_q"] or 1)):
        print(f"  {r['drug']:<14} {r['split']:<12} events {r['events']:>4}  "
              f"llr {fmt(r['llr'], '.2f')} (exp {fmt(r['llr_expected'], '.2f')}"
              f", q {fmt_p(r['llr_q'])})  G {fmt(r['g'], '.1f')} "
              f"(q {fmt_p(r['g_q'])})  {r['top_variant']}", file=out)
    f = [r["spectrum_factor"] for r in anyrows
         if r["spectrum_factor"] is not None and r["mutations"] > 0]
    if f:
        print(f"\n  spectrum_factor ({ANY}) across splits: "
              f"{min(f):.2f} to {max(f):.2f}", file=out)


ROW_COLS = ["reading", "drug", "split", "events", "mutations", "e_global",
            "e_spectrum", "spectrum_factor", "obs_over_spectrum", "llr",
            "llr_expected", "llr_p", "llr_q", "g", "g_p", "g_q",
            "top_variant"]
VAR_COLS = ["reading", "drug", "split", "variant", "observed",
            "split_events", "observed_share", "rest_share", "spectrum_share"]
DRUG_COLS = ["reading", "drug", "events", "variants_seen", "substitutions",
             "splits_tested", "llr_total", "llr_total_p", "g_all",
             "g_all_p"]


def write_rows(path, cols, rows):
    def f(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f(r.get(c)) for c in cols) + "\n")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def run_reading(np, tree, events, group, split_of, alts, contexts, masked,
                opp_ctx, sch, ref, internal_only, args, rng):
    cap = args.max_branch_length
    gids = sorted(set(group[i] for i in range(len(group))))
    spectra = group_spectra(np, tree, contexts, masked, group, cap,
                            internal_only, sch.n)
    total = sum(spectra.values()) if spectra else np.zeros(sch.n)
    row_of = {g: k for k, g in enumerate(gids)}
    labels = [split_of.get(g, "(outside)") for g in gids]
    exposure = np.array([spectra[g].sum() if g in spectra else 0.0
                         for g in gids])
    rate_tab = np.stack([rates(np, spectra[g] if g in spectra
                               else np.zeros(sch.n), total, opp_ctx,
                               args.prior) for g in gids])
    rate_of = lambda s, ctx: rate_tab[s, ctx]
    kept = [e for e in events if group[e[0]] >= -1
            and counted(tree, e[0], cap, internal_only)]
    rows, drug_rows, var_rows = [], [], []
    per_drug = defaultdict(list)
    seen = {}
    for e in kept:
        node, p, pb, a, d, g, name = e[:7]
        per_drug[d].append((row_of[group[node]], p, pb, a, f"{g}_{name}"))
        seen.setdefault((node, p, a), (row_of[group[node]], p, pb, a,
                                       f"{g}_{name}"))
    union = {}
    for d in alts:
        union.update(alts[d])
    for d in sorted(per_drug) + [ANY]:
        ev = list(seen.values()) if d == ANY else per_drug[d]
        al = union if d == ANY else alts.get(d, {})
        r, dr, vr = analyse_drug(np, d, ev, list(al), labels, row_of,
                                 exposure, rate_of, sch, ref, args, rng)
        rows.extend(x for x in r if x["split"] != "(outside)")
        drug_rows.append(dr)
        var_rows.extend(x for x in vr if x["split"] != "(outside)")
    return rows, drug_rows, var_rows


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    for opt in ("input_tree", "clade_samples", "catalogue", "annotation",
                "reference"):
        if not getattr(args, opt):
            print(f"error: --{opt} is required (or pass --selftest)",
                  file=sys.stderr)
            return 2
    try:
        import numpy as np
        import bte
    except ImportError as e:
        print(f"error: {e}; run this in the environment you use for "
              f"spectrumSplits.py", file=sys.stderr)
        return 1
    from permutation_null import extract
    out = sys.stdout

    def warn(msg):
        print(f"  WARNING: {msg}", file=out)

    groups = [int(x) for x in args.groups.split(",") if x.strip()]
    cat = read_catalogue(args.catalogue, groups)
    names = set()
    for d in cat.drugs():
        names |= cat.genes_of(d)
    features = load_annotation(args.annotation)
    genes = build_genes(cat, features, names, warn)
    targets = Targets(cat, genes)
    ref = read_fasta(args.reference)
    alts, unreached = targets.alternatives(ref)
    report_catalogue(cat, alts, unreached, out)
    for d in cat.drugs():
        if unreached[d]:
            print(f"  {d} unreached: " + ", ".join(
                f"{g}_{n}" for g, n in unreached[d][:8])
                + (" ..." if len(unreached[d]) > 8 else ""), file=out)

    masked = tn.read_mask(args.mask.split(",") if args.mask else [],
                          len(ref))
    sch = tn.scheme(args.fold)
    opp = np.array(sch.opportunity(ref, masked), dtype=float)
    opp_ctx = np.maximum(np.array([opp[sch.trinuc_of(c)]
                                   for c in range(sch.n)]), 1.0)

    print(f"Reading {args.input_tree} ...", file=sys.stderr)
    mat = bte.MATree(args.input_tree)
    tree = extract(mat.root)
    del mat
    clades = read_clade_samples(args.clade_samples)
    placed = place_splits(tree, clades)
    print(f"\nSplits: {len(clades)} in {args.clade_samples}", file=out)
    roots, split_of = set(), {}
    for s, (node, extra, missing) in sorted(placed.items()):
        if node is None:
            warn(f"{s}: none of its samples is on this tree; left out")
            continue
        if extra or missing:
            warn(f"{s}: placed at {tree.node_ids[node]}, whose clade has "
                 f"{extra} sample(s) not in the split; {missing} of its "
                 f"samples are not on this tree")
        if node in roots:
            warn(f"{s}: same node as {split_of[node]}; left out")
            continue
        roots.add(node)
        split_of[node] = s
    group = groups_of(tree, roots)

    print("Finding resistance events ...", file=sys.stderr)
    forward = defaultdict(set)
    for d in alts:
        for key in alts[d]:
            forward[key].add(d)
    events, reversions = find_events(tree, targets, ref, forward)
    print("Calling contexts ...", file=sys.stderr)
    contexts, _stats = sch.call(tree.parent, tree.children, tree.positions,
                                tree.types, ref, masked)

    by_drug = Counter(e[4] for e in events)
    distinct = len(set((e[0], e[1], e[3]) for e in events))
    print(f"\nResistance events: {distinct} distinct substitutions on "
          f"branches, {len(events)} event-drug pairs", file=out)
    for d, k in sorted(by_drug.items()):
        term = sum(1 for e in events if e[4] == d
                   and not tree.children[e[0]])
        routes = Counter(e[9] or "listed" for e in events if e[4] == d)
        print(f"  {d:<14} {k:>6}  terminal {term / k:.0%}  back-mutations "
              f"{reversions[d]:>4}  routes: "
              + ", ".join(f"{r} {n}" for r, n in routes.most_common(3)),
              file=out)

    rng = np.random.default_rng(args.seed)
    all_rows, all_drugs, all_vars = [], [], []
    for internal_only, label in ((False, "all"), (True, "internal")):
        rows, drug_rows, var_rows = run_reading(
            np, tree, events, group, split_of, alts, contexts, masked,
            opp_ctx, sch, ref, internal_only, args, rng)
        report_reading(f"{label} branches", rows, drug_rows, args, out)
        for r in rows + drug_rows + var_rows:
            r["reading"] = label
        all_rows += rows
        all_drugs += drug_rows
        all_vars += var_rows

    pre = args.prefix
    with open(f"{pre}_events.tsv", "w") as fh:
        fh.write("node\tterminal\tsplit\tpos\tparent_base\talt\tdrug\tgene\t"
                 "variant\teffect\tgroup\troute\n")
        for e in events:
            node = e[0]
            s = split_of.get(group[node], "(outside)")
            fh.write("\t".join(map(str, (
                tree.node_ids[node], int(not tree.children[node]), s)
                + tuple(e[1:9]) + (e[9] or "listed",))) + "\n")
    write_rows(f"{pre}.tsv", ROW_COLS, all_rows)
    write_rows(f"{pre}_drugs.tsv", DRUG_COLS, all_drugs)
    write_rows(f"{pre}_variants.tsv", VAR_COLS, all_vars)
    print(f"Wrote {pre}_events.tsv, {pre}.tsv, {pre}_drugs.tsv, "
          f"{pre}_variants.tsv", file=sys.stderr)
    return 0


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def selftest():
    import numpy as np
    from permutation_null import FlatTree
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + msg)
        ok = ok and cond

    # --- naming, forward and reverse strand --------------------------------
    #         1234567890123456789012345678
    ref = "GGGGGATGTCGTGGTAACCCCCC"     # + gene at 6..17: ATG TCG TGG TAA
    fwd = Gene("fw", "Rv1", 6, 17, 1, False)
    base = lambda p: ref[p - 1]
    nm = lambda g, p, a: [x[0] for x in snv_names(g, p, base, a)]
    check(nm(fwd, 10, "T") == ["c.5C>T", "p.Ser2Leu"], "+ strand missense")
    check(nm(fwd, 14, "A") == ["c.9G>A", "p.Trp3*"], "+ strand stop gained")
    check(nm(fwd, 7, "C") == ["c.2T>C", "p.Met1?"], "+ strand start lost")
    check(nm(fwd, 15, "C") == ["c.10T>C", "p.Ter4Glnext*?"],
          "+ strand stop lost")
    check(nm(fwd, 11, "A") == ["c.6G>A", "p.Ser2="], "+ strand synonymous")
    check(nm(fwd, 3, "T") == ["c.-3G>T"], "+ strand upstream")
    # reverse strand: the same gene read from the other strand
    rc = "".join(COMPLEMENT[b] for b in reversed(ref))   # gene at 7..18
    rev = Gene("rv", "Rv2c", 7, 18, -1, False)
    base_r = lambda p: rc[p - 1]
    nr = lambda p, a: [x[0] for x in snv_names(rev, p, base_r, a)]
    # genomic 14 on rc is coding 5 (C on the coding strand, G on +)
    check(nr(14, "A") == ["c.5C>T", "p.Ser2Leu"], "- strand missense")
    check(nr(21, "A") == ["c.-3G>T"], "- strand upstream")
    rna = Gene("rr", "Rvnr1", 6, 17, 1, True)
    check(nm(rna, 10, "T") == ["n.5C>T"], "rRNA naming")
    # parent background changes the name
    bg = lambda p: "A" if p == 9 else ref[p - 1]        # TCG -> ACG (Thr)
    check([x[0] for x in snv_names(fwd, 10, bg, "T")][1] == "p.Thr2Met",
          "named against the parent's codon")

    # --- catalogue precedence ---------------------------------------------
    rows = [
        {"drug": "D", "gene": "fw", "variant": "fw_p.Ser2Leu",
         "effect": "missense_variant", "grade": 1, "route": None},
        {"drug": "D", "gene": "fw", "variant": "fw_LoF", "effect": "LoF",
         "grade": 2, "route": None},
        {"drug": "D", "gene": "fw", "variant": "fw_p.Trp3*",
         "effect": "stop_gained", "grade": 3, "route": None},
        {"drug": "D", "gene": "fw", "variant": "fw_p.Ser2Pro",
         "effect": "missense_variant", "grade": 4, "route": None},
        {"drug": "D", "gene": "fw", "variant": "fw_c.-3G>T",
         "effect": "upstream_gene_variant", "grade": 2, "route": None},
        {"drug": "R", "gene": "RRDR", "variant": "RRDR", "effect": "RRDR",
         "grade": 1, "route": None},
        # listed only to show that a reversion is never counted
        {"drug": "D", "gene": "fw", "variant": "fw_p.Leu2Ser",
         "effect": "missense_variant", "grade": 2, "route": None},
    ]
    cat = Catalogue(rows, {"fw": ("Rv1", 6, 17)}, [1, 2])
    check(cat.classify("D", "fw", "p.Ser2Leu", "missense_variant", 2)
          == (1, "listed"), "listed group 1")
    check(cat.classify("D", "fw", "p.Ser2Pro", "missense_variant", 2)
          is None, "listed group 4 is not resistance")
    check(cat.classify("D", "fw", "p.Trp3*", "stop_gained", 3) is None,
          "explicit group 3 beats the LoF rule")
    check(cat.classify("D", "fw", "p.Ser2*", "stop_gained", 2)
          == (2, "Rule: LoF (not listed)"), "LoF rule covers unlisted stops")
    check(cat.classify("R", "rpoB", "p.Ser450Leu", "missense_variant", 450)
          == (1, "Rule: RRDR (not listed)"), "RRDR rule inside 426-452")
    check(cat.classify("R", "rpoB", "p.Val170Phe", "missense_variant", 170)
          is None, "RRDR rule stops outside the region")
    targets = Targets(cat, {"fw": fwd})
    alts, unreached = targets.alternatives(ref)
    names_d = sorted(set(n for _g, n in alts["D"].values()))
    check("p.Ser2Leu" in names_d and "c.-3G>T" in names_d
          and "p.Met1?" in names_d and "p.Trp3*" not in names_d,
          f"alternatives {names_d}")

    # --- events and split placement on a small tree ------------------------
    t = FlatTree()

    def add(nid, parent, muts):
        i = len(t.node_ids)
        t.node_ids.append(nid)
        t.parent.append(parent)
        t.children.append([])
        if parent >= 0:
            t.children[parent].append(i)
        t.raw_len.append(len(muts))
        t.positions.append(array("i", [p for p, _ in muts]))
        t.types.append(array("b", [tn.MUT_INDEX[m] for _, m in muts]))
        return i

    r0 = add("root", -1, [])
    a = add("a", r0, [(9, "TA")])               # background Ser -> Thr
    add("s1", a, [(10, "CT")])                  # Thr2Met: not listed
    add("s2", a, [])
    b = add("b", r0, [(10, "CT")])              # Ser2Leu, group 1
    add("s3", b, [(14, "GA")])                  # Trp3* graded 3: none
    add("s4", b, [(8, "GT")])                   # Met1? via the LoF rule
    add("s5", b, [(10, "TC")])                  # Leu2Ser: a reversion
    fwd_alts = {key: {"D"} for key in alts["D"]}
    ev, rev = find_events(t, targets, ref, fwd_alts)
    got = sorted((t.node_ids[e[0]], e[6]) for e in ev)
    check(got == [("b", "p.Ser2Leu"), ("s4", "p.Met1?")],
          f"events named against the parent {got}")
    check(rev == Counter({"D": 1}),
          f"Leu2Ser under Ser2Leu is a back-mutation, not an event "
          f"{dict(rev)}")
    on_leu = lambda p: "T" if p == 10 else ref[p - 1]
    check(restores_reference(fwd, 10, on_leu, "C", base)
          and not restores_reference(fwd, 10, base, "T", base)
          and restores_reference(rna, 10, on_leu, "C", base),
          "restores_reference: codon and nucleotide")
    placed = place_splits(t, {"A": {"s1", "s2"}, "B": {"s3"},
                              "X": {"s1", "s3"}})
    check(placed["A"][:2] == (a, 0) and placed["B"][:2] == (5, 0)
          and placed["X"] == (r0, 3, 0), f"placement {placed}")

    # --- composition statistics -------------------------------------------
    rng = np.random.default_rng(3)
    S, V, n = 4, 30, 300
    logrho = rng.normal(0, 0.5, (S, V))          # each split's spectrum ratios
    tested = np.ones(S, dtype=bool)

    def simulate(shift):
        """Events per split, P(v | s) proportional to exp(shift[s])."""
        labels, cats = [], []
        for s in range(S):
            w = np.exp(shift[s])
            cats += list(rng.choice(V, size=n, p=w / w.sum()))
            labels += [s] * n
        return np.array(labels), np.array(cats)

    def table(lab, cats):
        return np.bincount(lab * V + cats, minlength=S * V).reshape(S, V) \
            .astype(float)

    lab, cats = simulate(logrho)
    res = shuffle_test(np, lab, cats, S, V, logrho, tested, 200, rng)
    rp, rtot = rho_test(np, table(lab, cats), logrho, tested, 400, rng)
    check(rp.max() < 0.01 and rtot < 0.01 and res[0][0][0] > 0,
          f"spectrum-driven composition found (llr p {rp.max():.3g})")
    check(res[1].max() < 0.01, f"and G finds it ({res[1].max():.3g})")
    C = table(lab, cats)
    rest = C[1:].sum(0)
    q0 = (rest + 0.5) / (rest.sum() + 0.5 * V)
    w = q0 * np.exp(logrho[0])
    w /= w.sum()
    draws = rng.multinomial(n, w, size=2000)
    mean_llr = (draws * np.log(w / q0)).sum(1).mean()
    check(abs(res[0][4][0] - mean_llr) < 0.5,
          f"llr_expected ({res[0][4][0]:.2f}) is the mean llr of events "
          f"drawn from q' ({mean_llr:.2f})")
    vr = variant_rows(np, "D", ["a", "b", "c", "d"], C, logrho, tested,
                      [f"v{v}" for v in range(V)], min_share=0.0)
    from_rows = sum(r["observed"] * math.log(r["spectrum_share"]
                                             / r["rest_share"])
                    for r in vr if r["split"] == "a")
    check(abs(from_rows - res[0][0][0]) < 1e-6,
          "the variants table's shares reproduce llr")
    # A composition shifted for some other reason (treatment, lineage),
    # unrelated to the spectrum.  Shuffling split labels calls the spectrum
    # responsible far too often; permuting the ratios across variants does not.
    reps, low_rho, low_lab = 100, 0, 0
    for _ in range(reps):
        lab, cats = simulate(rng.normal(0, 1.0, (S, V)))
        C = table(lab, cats)
        rp, _t = rho_test(np, C, logrho, tested, 99, rng)
        low_rho += rp[0] <= 0.05
        obs = composition_stats(np, C, logrho, tested)[0][0]
        ge = sum(composition_stats(np, table(rng.permutation(lab), cats),
                                   logrho, tested)[0][0] >= obs
                 for _ in range(99))
        low_lab += (ge + 1) / 100 <= 0.05
    check(low_rho / reps <= 0.11,
          f"other-cause shift: ratio permutation p <= 0.05 in "
          f"{low_rho}/{reps}")
    check(low_lab / reps >= 0.2,
          f"(label shuffling would have said so in {low_lab}/{reps})")
    # G and llr agree with direct computation for one table
    C = np.array([[5., 1, 0], [2, 4, 3]])
    lr = np.zeros((2, 3))
    llr, G, _t, Gall, _e = composition_stats(np, C, lr, np.array([1, 1],
                                                                 bool))
    rest = C[1]
    E = np.vstack([C.sum(0) * C[0].sum() / C.sum(),
                   C.sum(0) * C[1].sum() / C.sum()])
    m = C > 0
    G_direct = 2 * (C[m] * np.log(C[m] / E[m])).sum()
    check(abs(G[0] - G_direct) < 1e-9 and abs(Gall - G_direct) < 1e-9
          and abs(llr).max() < 1e-12, "G matches the textbook 2 x V G; "
          "llr is 0 with no spectrum difference")

    # --- rate expectations ------------------------------------------------
    cnt = np.array([[10., 0, 0], [0, 10, 0]])
    rt = np.stack([rates(np, c, cnt.sum(0), np.ones(3), 0.0) for c in cnt])
    check(np.allclose(rt[0], [1, 0, 0]), "rates are shares over opportunity")
    print("all ok" if ok else "SOME TESTS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
