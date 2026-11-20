#!/usr/bin/env python3
"""Genotype-first test: does a candidate gene's mutation shift the spectrum below it?

The problem it addresses.  Reading a mechanism off a spectrumSplits.py split
fails for a structural reason: a deep split's defining branch carries every
lineage-defining mutation, so any repair-gene variant on it co-occurs with the
spectrum by construction.  (On the MTBC tree, the mutT4 R48G hit on the
lineage 2.2.1 branch was exactly this.)  One branch cannot separate a variant
from its lineage.  Independent recurrence can: if a gene matters, then
*every* branch where it was hit should, on average, show the same shift below
it relative to its closest relatives -- and within each such comparison the
lineage is held fixed.

This script starts from the genotype rather than the splits.  It finds every
branch carrying a coding change in a named gene, compares what accumulated
below that branch with what accumulated below its siblings, and pools the
comparisons across all the independent hits on the tree.

POLYTOMIES.  UShER trees are not binary; a parent can have hundreds of
children.  So there is no single "sister clade".  Each hit branch h with
parent P is compared against ALL of P's other children, pooled: every child of
P is a lineage descended from the same ancestral genome, so each is an
equally valid contrast.  The pooled comparator is just their sum.

The null is exchangeability among P's children.  Under "the hit does
nothing", the hit could have landed on any child of P, so the test asks how
unusual the real hit child is among its siblings.  A permutation draws, for
every event independently, a pseudo-hit child from P's children and
recomputes the pooled statistic.  For a binary parent that is exactly a
paired sign-flip test; for a polytomy it is its natural generalisation, and
the size of the polytomy is irrelevant to its validity.

The draw is size-biased.  A mutation is more likely to land on a long branch
than a short one, so the real hit branch is longer than a random sibling on
average.  If draws were uniform, the hit child's subtree would look like it
had more mutations than its siblings even when nothing changed.  So each
child is drawn with probability proportional to (its top-branch SNV count,
not counting the hit itself) + 1 -- which is the null probability that the
hit would have landed there.

What is counted.  For each child c of P: the mutations on c's branch and
everything below it (--top_branch include, the default), or only below it
(--top_branch exclude).  The hit mutation itself is always removed.
Including the top branch is what lets terminal hits count at all -- most
branches are terminal -- at the cost of diluting the hit child with mutations
that may predate the hit on the same branch.  The dilution weakens a real
effect; it does not create a false one, because the siblings are treated
identically.  Use exclude as a cleaner check on internal hits only.

Nested hits.  If the same gene is hit again inside a hit subtree, the inner
subtree is cut out of the outer comparison and tested as its own event, and
its row names the outer one in `nested_under`.  Other hit children of the
same parent are left out of each other's comparators.

Statistics, per test (one gene and one effect set, plus a pooled ALL row):

  diff_<class>   weighted mean, across events, of (hit-child proportion -
                 comparator proportion) for that class.  Weights are
                 n_hit * n_comp / (n_hit + n_comp), the Mantel-Haenszel
                 style weight, so a terminal hit with three mutations moves
                 the estimate far less than an internal hit with three
                 hundred.
  p_<class>      two-sided permutation p, centred on the permutation mean.
  omnibus_p      one statistic over all classes: sum of squared centred diffs
                 scaled by the tree-wide class proportions.  Use it to decide
                 whether to look at the per-class columns at all.
  effective_events
                 (sum of weights)^2 / sum of squared weights.  The weights
                 grow with subtree size, so one deep branch can carry most
                 of a test; this is how many equally weighted events the
                 test is worth.  A handful means the result is about those
                 few branches, whatever the event count.

No rate test.  A mutator should also raise the mutation count below the hit,
but this permutation cannot test that honestly.  The size-biased draw gets
the chance of each child carrying the hit right, but not the extra branch
time a real hit implies, and depth depends on that time.  On the MTBC tree,
hit children came out shallower than the draw expects in both effect sets,
and a depth test failed the synonymous control in 31 of 56 genes.  The
spectrum statistics are proportions, which barely depend on depth, and they
pass that control.

Controls built in.  Synonymous hits in the same genes are tested the same way
by default.  They should show nothing.  If they do, the comparison itself is
biased, and the nonsynonymous result cannot be read either.

One run, three readings.  The list is tested with every hit counted and
(--top_branch both, the default) with internal hits only, where per-genome
error cannot reach; and, with --control, a length-matched control list from
control_genes.py is tested with every hit counted, as the baseline the
list's excess is read against.  The tree is read once.  Each reading starts
its random draws from --seed, so it gives exactly what a run of it alone
would.  The report on stdout then reads the summaries together: is the
synonymous control clean, do the nonsynonymous hits exceed chance and the
control, which genes clear the bar (0.05 over the genes tested), and which
are leads below it.

Outputs, named from --prefix (default hit): PREFIX_summary.tsv and
PREFIX_events.tsv, the same with _internal and _control, and with
--null_draws_genes PREFIX_null_draws.tsv (and _internal), every permutation's
class shifts for those genes, for plot_gene_hits.py --prefix.

What it cannot see.  The tree holds SNVs only, so frameshifts, deletions and
IS6110 insertions -- the commonest ways a gene is actually lost -- are
invisible here.  Promoter changes are not tested.  And the ALL row pools
genes whose mechanisms differ: a gene that raises CT|GA and one that raises
CA|GT cancel there.  ALL is a screen.  The evidence is in the per-gene rows.

Trinucleotide context (--context_genes).  For the genes named, in both
readings, the same events, comparators and controls, counted in
trinucleotide context classes instead of mutation classes, to ask where
inside a class a shift falls.  The summary is PREFIX_summary_context.tsv
(and _internal), and the report ends with each gene's shape p values.  With --fold, the 96
pyrimidine-centred classes (6 classes x 16 flank pairs); without it, the 192
forward-strand classes (12 types x 16), which keep the strands apart as an
unfolded run does.  Needs numpy.  Every summary column of the class-level
reading is still written, computed by summing the contexts of
each class, plus:

  shape_p_<class>  does the class's shift fall on its 16 contexts in
                   proportion to how common each context is in the tree?
                   The statistic is the squared distance between each
                   context's shift and its proportional share of the class
                   shift, scaled by the context's tree-wide proportion, and
                   the p is against the same permutation draws.  A small p
                   says the shift is concentrated in particular contexts;
                   which ones is in PREFIX_context.tsv.  "Proportional" uses the
                   tree-wide context mix, so a comparator whose own mix
                   differs a little from the tree's adds a second-order term.
  shape_omnibus_p  the same, summed over the classes.

PREFIX_context.tsv (and _internal) has one row per test and context: the context's share of
its class tree-wide, its shift, the shift a proportional spread would give,
the excess over that, and a two-sided permutation p for the shift itself.
Contexts that cannot be called -- a masked or unknown flank, or an adjacent
change on the same branch -- are left out of every count, so the class
totals differ slightly from the class-level reading, and the draws differ, so do not
expect its p values to match one digit for digit.

Usage:
    gene_hit_contrasts.py --input_tree tree.pb.gz \\
        --annotation NC_000962.3.gff3.gz --reference NC_000962.3.fna \\
        --genes tb_repair_genes.txt --control control_genes.txt \\
        --fold --max_branch_length 250 --permutations 10000 \\
        --null_draws_genes dnaE2 --context_genes dnaE2 --mask masked.bed \\
        > hit_summary.txt

    gene_hit_contrasts.py --selftest
"""

import argparse
import bisect
import csv
import math
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
### Imported rather than reimplemented so the folding, the parse rules, the
### codon arithmetic and the chi-squared cannot drift from the other tools.
from permutation_null import (MUT_TYPES, chi2_2xk, extract, folded_index,
                              own_spectra)
from annotate_sites import (codon_for, load_annotation, load_fasta,
                            resolve_region_names, smart_open, translate)

DNA_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

### Bacterial start codons (NCBI table 11 allows more, but these three are
### essentially all real MTB starts).  A first codon moving between two of
### them still initiates, so that is synonymous rather than start_lost.
START_CODONS = {"ATG", "GTG", "TTG"}

### Most severe first.  An event (one branch, one gene) takes the most severe
### effect among its codon changes.
SEVERITY = ["nonsense", "start_lost", "stop_lost", "missense", "synonymous",
            "unknown"]

EFFECT_SETS = {
    "nonsynonymous": {"nonsense", "start_lost", "stop_lost", "missense"},
    "lof": {"nonsense", "start_lost"},
    "missense": {"missense"},
    "synonymous": {"synonymous"},
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Test whether coding changes in candidate genes shift the "
                    "mutation spectrum below them, against sibling subtrees.")
    p.add_argument("--input_tree", type=str, default=None,
                   help="Input tree file (protobuf format)")
    p.add_argument("--annotation", type=str, default=None,
                   help="GFF3 or GenBank for the reference (gzip ok)")
    p.add_argument("--reference", type=str, default=None,
                   help="Reference FASTA the tree's coordinates refer to")
    p.add_argument("--seqid", type=str, default=None,
                   help="Sequence to use if the FASTA/annotation has several")
    p.add_argument("--genes", type=str, default=None,
                   help="Gene names / locus tags, one per line, '#' comments "
                        "allowed (tb_repair_genes.txt works as is)")
    p.add_argument("--effects", type=str, default="nonsynonymous,synonymous",
                   help="Comma-separated effect sets to test, from: "
                        + ", ".join(EFFECT_SETS) + ". Default "
                        "nonsynonymous,synonymous (the second is the control)")
    p.add_argument("--fold", action="store_true",
                   help="Score the strand-folded 6-type spectrum instead of 12")
    p.add_argument("--max_branch_length", type=int, default=None,
                   help="Match the search: branches above this contribute no "
                        "mutations (they still count for topology)")
    p.add_argument("--top_branch", choices=("both", "include", "exclude"),
                   default="both",
                   help="include: count each child's own branch (needed for "
                        "terminal hits); exclude: only what lies below it, "
                        "so only internal hits count; both (default): run "
                        "each, the second as the internal-only reading")
    p.add_argument("--control", type=str, default=None,
                   help="A control gene list (control_genes.py --output), "
                        "tested the same way with every hit counted, as the "
                        "baseline the list's excess is read against")
    p.add_argument("--min_focal", type=int, default=1,
                   help="Minimum mutations in the hit child's counted set for "
                        "the event to be used; siblings below this are not "
                        "drawn as pseudo-hits (default 1)")
    p.add_argument("--min_comparator", type=int, default=10,
                   help="Minimum mutations in the pooled siblings for the "
                        "event to be used (default 10)")
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--prefix", type=str, default="hit",
                   help="Prefix of the output TSVs (default hit: "
                        "hit_summary.tsv, hit_events.tsv, ...; see the "
                        "module docstring)")
    p.add_argument("--null_draws_genes", type=str, default=None,
                   help="Comma-separated genes (or ALL) whose every "
                        "permutation's class shifts are written to "
                        "PREFIX_null_draws.tsv (and _internal), for "
                        "plot_gene_hits.py --prefix")
    p.add_argument("--context_genes", type=str, default=None,
                   help="Comma-separated genes from --genes to test again in "
                        "trinucleotide context classes: 96 folded with "
                        "--fold, 192 forward-strand without (needs numpy); "
                        "see the module docstring")
    p.add_argument("--mask", type=str, default=None,
                   help="With --context_genes: comma-separated BED files the "
                        "tree was masked with; a mutation with a masked flank "
                        "has an unknown context and is not counted")
    p.add_argument("--calibrate_context", type=int, default=0, metavar="N",
                   help="Instead of a test on a tree: simulate N small trees "
                        "per scenario and report how often the --context "
                        "tests reach p < 0.05 with no shift, with a whole "
                        "class shifted in proportion, and with one context "
                        "shifted.  Honours --fold, --permutations and --seed; "
                        "needs numpy, not bte")
    p.add_argument("--selftest", action="store_true",
                   help="Check subtree bookkeeping, codon effects and the "
                        "permutation null on synthetic trees, then exit")
    return p.parse_args(argv)


# --------------------------------------------------------------------------
# genes and hits
# --------------------------------------------------------------------------

def read_gene_list(path):
    names = []
    seen = set()
    with smart_open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                name = line.split()[0]
                if name.lower() not in seen:
                    seen.add(name.lower())
                    names.append(name)
    return names


def cds_index_for(features, names):
    """Map each position in a requested gene's CDS to (gene_name, feature).

    Two requested names can resolve to the same CDS (alkA and its synonym
    ada, if a list gives both).  Counting it twice would make
    every hit two events, so each CDS belongs to the first name that claims it.
    """
    matched, missing = resolve_region_names(features, names)
    owner = {}
    aliases = []
    no_cds = []
    for name in names:
        cds = [f for f in matched.get(name, []) if f.ftype == "CDS"]
        if name in matched and not cds:
            no_cds.append(name)
        for f in cds:
            key = (f.start, f.end, f.strand)
            if key in owner:
                if owner[key][0] != name:
                    aliases.append((name, owner[key][0]))
                continue
            owner[key] = (name, f)
    index = defaultdict(list)
    for name, f in owner.values():
        for pos in range(f.start, f.end + 1):
            index[pos].append((name, f))
    genes = []
    for name, _f in owner.values():
        if name not in genes:
            genes.append(name)
    return index, genes, missing, sorted(set(aliases)), no_cds


def base_at(tree, node, pos, ref):
    """The base at pos in node's genome: nearest mutation at or above, else ref."""
    x = node
    while x >= 0:
        for p, t in zip(tree.positions[x], tree.types[x]):
            if p == pos:
                return MUT_TYPES[t][1]
        x = tree.parent[x]
    return ref[pos - 1] if 1 <= pos <= len(ref) else "N"


def codon_positions(f, codon_no):
    """Genomic positions of a codon, in codon order."""
    if f.strand > 0:
        s = f.start + 3 * (codon_no - 1)
        return [s, s + 1, s + 2]
    e = f.end - 3 * (codon_no - 1)
    return [e, e - 1, e - 2]


def codon_string(f, bases):
    if f.strand > 0:
        return "".join(bases)
    return "".join(DNA_COMPLEMENT.get(b, "N") for b in bases)


def classify_codon(codon_no, ref_codon, alt_codon):
    if len(ref_codon) != 3 or "N" in ref_codon or "N" in alt_codon:
        return "unknown", "?"
    if codon_no == 1 and ref_codon in START_CODONS:
        if alt_codon in START_CODONS:
            return "synonymous", f"M1M"
        ra = "M"
        aa = translate(alt_codon) or "X"
        return "start_lost", f"M1{aa}"
    ra = translate(ref_codon) or "X"
    aa = translate(alt_codon) or "X"
    label = f"{ra}{codon_no}{aa}"
    if ra == aa:
        return "synonymous", label
    if aa == "*":
        return "nonsense", label
    if ra == "*":
        return "stop_lost", label
    return "missense", label


def find_hit_events(tree, index, ref):
    """One event per (branch, gene) carrying at least one SNV in that CDS.

    The codon's other bases are taken from the PARENT's genome, reconstructed
    up the tree, not from the reference: a codon already changed on an
    ancestral branch would otherwise be translated wrongly.  All of the
    branch's changes within one codon are applied together.
    """
    by_key = defaultdict(list)
    for i in range(len(tree.node_ids)):
        for k, pos in enumerate(tree.positions[i]):
            hits = index.get(pos)
            if hits:
                t = tree.types[i][k]
                for name, f in hits:
                    by_key[(i, name)].append((pos, t, f, k))

    events = []
    ref_mismatch = 0
    for (node, gene), muts in by_key.items():
        parent = tree.parent[node]
        codons = defaultdict(dict)
        labels = []
        for pos, t, f, _k in muts:
            mt = MUT_TYPES[t]
            labels.append(f"{mt[0]}{pos}{mt[1]}")
            if base_at(tree, parent, pos, ref) != mt[0]:
                ref_mismatch += 1
            codon_no, _cp, _ = codon_for(f, pos, None)
            codons[(id(f), codon_no)][pos] = mt[1]
            codons[(id(f), codon_no)]["_f"] = f
        effects = []
        changes = []
        for (_fid, codon_no), alts in codons.items():
            f = alts.pop("_f")
            gpos = codon_positions(f, codon_no)
            before = [base_at(tree, parent, p, ref) for p in gpos]
            after = [alts.get(p, b) for p, b in zip(gpos, before)]
            eff, lab = classify_codon(codon_no, codon_string(f, before),
                                      codon_string(f, after))
            effects.append(eff)
            changes.append(lab)
        worst = min(effects, key=SEVERITY.index)
        events.append({
            "node": node, "gene": gene, "effect": worst,
            "aa_changes": ",".join(changes), "mutations": ",".join(labels),
            "hit_types": [t for _p, t, _f, _k in muts],
            "hit_k": [k for _p, _t, _f, k in muts],
        })
    return events, ref_mismatch


# --------------------------------------------------------------------------
# subtree bookkeeping
# --------------------------------------------------------------------------

def subtree_sums(tree, own, ncat):
    """Per node, including its own branch: class counts and tips."""
    n = len(tree.node_ids)
    V = [None] * n
    T = [0] * n
    for i in range(n - 1, -1, -1):
        v = list(own[i])
        kids = tree.children[i]
        t = 0 if kids else 1
        for c in kids:
            vc = V[c]
            for k in range(ncat):
                v[k] += vc[k]
            t += T[c]
        V[i], T[i] = v, t
    return V, T


def cut_adjustments(tree, V, T, hits, ncat):
    """What each ancestor of a hit must subtract to exclude the hit's subtree.

    Walks up from each hit and stops at the first hit above it, which is
    where the nested subtree stops being anyone else's business: everything
    above that hit already loses the whole outer subtree.  Returns
    adj[node] = [hits to subtract, tips] and nested[hit] = enclosing hit.
    The vectors are summed only for the nodes child_value() is asked about,
    which with 96 context classes is most of the cost saved.
    """
    adj = {}
    nested = {}
    for h in hits:
        x = tree.parent[h]
        while x >= 0:
            a = adj.get(x)
            if a is None:
                a = adj[x] = [[], 0]
            a[0].append(h)
            a[1] += T[h]
            if x in hits:
                nested[h] = x
                break
            x = tree.parent[x]
    return adj, nested


def child_value(c, V, T, own, adj, top, ncat, remove=None):
    """(vec, tips) for child c as seen from its parent."""
    v = list(V[c])
    t = T[c]
    a = adj.get(c)
    if a is not None:
        for h in a[0]:
            vh = V[h]
            for k in range(ncat):
                v[k] -= vh[k]
        t -= a[1]
    if top == "exclude":
        for k in range(ncat):
            v[k] -= own[c][k]
    elif remove is not None:
        for k in range(ncat):
            v[k] -= remove[k]
    return v, t


# --------------------------------------------------------------------------
# events and the test
# --------------------------------------------------------------------------

class Prepared(object):
    """One event's eligible children, with everything a draw needs."""
    __slots__ = ("wd", "w", "cum", "total", "wd2", "w2")


def prepare_event(tree, h, hits, V, T, own, adj, top, ncat,
                  hit_vec, min_focal, min_comparator):
    """Build the hit child's contrast with its pooled siblings.

    Returns (row_fields, Prepared or None).  Prepared is None when the event
    is not informative.
    """
    P = tree.parent[h]
    if P < 0:
        return None, None
    fv, ft = child_value(h, V, T, own, adj, top, ncat, remove=hit_vec)
    sibs = [c for c in tree.children[P] if c != h and c not in hits]
    vals = [child_value(c, V, T, own, adj, top, ncat) for c in sibs]
    tot_v = list(fv)
    tot_t = ft
    for v, t in vals:
        for k in range(ncat):
            tot_v[k] += v[k]
        tot_t += t
    cv = [tot_v[k] - fv[k] for k in range(ncat)]
    ct = tot_t - ft
    nf, nc = sum(fv), sum(cv)
    row = {
        "parent_id": tree.node_ids[P],
        "terminal": int(not tree.children[h]),
        "parent_children": len(tree.children[P]),
        "siblings_used": len(sibs),
        "focal_tips": ft, "focal_mutations": nf,
        "comparator_tips": ct, "comparator_mutations": nc,
        "chi2": chi2_2xk(fv, cv) if nf and nc else 0.0,
        "focal_vec": fv, "comparator_vec": cv,
    }
    if nf < min_focal or nc < min_comparator:
        row["informative"] = 0
        return row, None
    row["informative"] = 1

    ### Eligible pseudo-hit children: the real one first, then every sibling
    ### that would have passed the same --min_focal had the hit landed there.
    ### Size-bias weight: non-hit SNVs on the top branch, plus one.
    hit_count = sum(hit_vec) if hit_vec is not None else 0
    elig = [(fv, len(tree.types[h]) - hit_count + 1)]
    for c, (v, t) in zip(sibs, vals):
        if sum(v) >= min_focal:
            elig.append((v, len(tree.types[c]) + 1))

    prep = Prepared()
    prep.wd, prep.w, prep.cum = [], [], []
    ### With context classes (96 or 192), also the within-class version: each
    ### class's 16 context shares, weighted by that class's counts.  See
    ### ContextAccumulator for why the shape test needs it.
    blocks = ncat // 16 if ncat % 16 == 0 else 0
    prep.wd2, prep.w2 = [], []
    acc = 0
    for v, wt in elig:
        rest = [tot_v[k] - v[k] for k in range(ncat)]
        if blocks:
            wd2, w2 = [], []
            for b in range(blocks):
                lo, hi = 16 * b, 16 * b + 16
                m1, m2 = sum(v[lo:hi]), sum(rest[lo:hi])
                if m1 > 0 and m2 > 0:
                    w = m1 * m2 / (m1 + m2)
                    wd2.extend(w * (v[k] / m1 - rest[k] / m2)
                               for k in range(lo, hi))
                    w2.append(w)
                else:
                    wd2.extend([0.0] * 16)
                    w2.append(0.0)
            prep.wd2.append(wd2)
            prep.w2.append(w2)
        n1, n2 = sum(v), sum(rest)
        if n1 > 0 and n2 > 0:
            w = n1 * n2 / (n1 + n2)
            prep.wd.append([w * (v[k] / n1 - rest[k] / n2)
                            for k in range(ncat)])
            prep.w.append(w)
        else:
            prep.wd.append([0.0] * ncat)
            prep.w.append(0.0)
        acc += wt
        prep.cum.append(acc)
    prep.total = acc
    return row, prep


def statistic(preps, choice, ncat):
    acc = [0.0] * ncat
    W = 0.0
    for p, j in zip(preps, choice):
        w = p.w[j]
        if w > 0:
            wd = p.wd[j]
            for k in range(ncat):
                acc[k] += wd[k]
            W += w
    return [a / W for a in acc] if W > 0 else [0.0] * ncat


def permutation_test(preps, ncat, pbar, permutations, rng):
    """Observed statistics and permutation p-values; see the module docstring."""
    obs_d = statistic(preps, [0] * len(preps), ncat)
    perm_d = []
    for _ in range(permutations):
        choice = [bisect.bisect_right(p.cum, rng.random() * p.total)
                  for p in preps]
        ### bisect_right on a value in [0, total) never returns len(cum).
        perm_d.append(statistic(preps, choice, ncat))
    if not permutations:
        return obs_d, [float("nan")] * ncat, float("nan"), float("nan"), []
    mean_d = [sum(d[k] for d in perm_d) / permutations for k in range(ncat)]

    def omni(d):
        return sum((d[k] - mean_d[k]) ** 2 / pbar[k]
                   for k in range(ncat) if pbar[k] > 0)

    p_class = []
    for k in range(ncat):
        o = abs(obs_d[k] - mean_d[k])
        ge = sum(1 for d in perm_d if abs(d[k] - mean_d[k]) >= o - 1e-15)
        p_class.append((ge + 1) / (permutations + 1))
    o = omni(obs_d)
    p_omni = (sum(1 for d in perm_d if omni(d) >= o - 1e-15) + 1) / \
        (permutations + 1)
    return obs_d, p_class, o, p_omni, perm_d


def run_test(tree, test_events, V, T, own, ncat, args, pbar, rng,
             accumulator=None):
    """Events for one test (they define the cut set), then the permutation."""
    hits = set(e["node"] for e in test_events)
    adj, nested = cut_adjustments(tree, V, T, hits, ncat)
    cap = args.max_branch_length
    rows = []
    preps = []
    w0 = []
    ### Several genes can be hit on one branch; within one test a branch is
    ### one event, or the ALL row would count it twice.
    seen = set()
    for e in test_events:
        h = e["node"]
        if h in seen:
            continue
        seen.add(h)
        hit_vec = [0] * ncat
        if cap is None or tree.raw_len[h] <= cap:
            for c in e["hit_cats"]:
                hit_vec[c] += 1
        row, prep = prepare_event(tree, h, hits, V, T, own, adj,
                                  args.top_branch, ncat, hit_vec,
                                  args.min_focal, args.min_comparator)
        if row is None:
            continue
        row.update({"node_id": tree.node_ids[h], "gene": e["gene"],
                    "effect": e["effect"], "aa_changes": e["aa_changes"],
                    "mutations": e["mutations"],
                    "nested_under": tree.node_ids[nested[h]]
                    if h in nested else ""})
        rows.append(row)
        if prep is not None:
            w0.append(prep.w[0])
            ### With context classes, holding every event's draws until the end
            ### costs gigabytes on a large list; the accumulator consumes
            ### each event as it is built.
            if accumulator is not None:
                accumulator.add(prep)
            else:
                preps.append(prep)
    result = None
    if w0:
        result = accumulator.finish(pbar) if accumulator is not None else \
            permutation_test(preps, ncat, pbar, args.permutations, rng)
    return rows, len(w0), result, effective_events(w0)


def effective_events(ws):
    """(sum w)^2 / sum w^2 over the real hit children's weights.

    The weights grow with subtree size, so one deep branch can carry most of
    a test.  This says how many equally weighted events the test is worth.
    """
    s2 = sum(w * w for w in ws)
    return sum(ws) ** 2 / s2 if s2 > 0 else 0.0


class ContextAccumulator(object):
    """The permutation test of permutation_test(), for context classes
    (96 folded or 192 forward-strand: any multiple of 16).

    Vectorised with numpy and fed one event at a time.  Each event's draws
    for all permutations are taken at once, and its contribution added to a
    permutations x classes array, so memory does not grow with the number
    of events.  The statistic is the same weighted mean shift: weights depend
    only on mutation totals, so summing a class's 16 contexts gives exactly
    the 6-class statistic for the same draws.
    """

    def __init__(self, ncat, permutations, rng):
        import numpy as np
        self.np = np
        self.ncat = ncat
        self.P = permutations
        self.rng = rng
        self.acc = np.zeros((permutations, ncat))
        self.W = np.zeros(permutations)
        self.obs = np.zeros(ncat)
        self.W0 = 0.0
        ncls = ncat // 16
        self.acc2 = np.zeros((permutations, ncat))
        self.Wc = np.zeros((permutations, ncls))
        self.obs2 = np.zeros(ncat)
        self.Wc0 = np.zeros(ncls)

    def add(self, prep):
        np = self.np
        wd = np.asarray(prep.wd, dtype=float)
        w = np.asarray(prep.w, dtype=float)
        wd2 = np.asarray(prep.wd2, dtype=float)
        w2 = np.asarray(prep.w2, dtype=float)
        if w[0] > 0:
            self.obs += wd[0]
            self.W0 += w[0]
        self.obs2 += wd2[0]
        self.Wc0 += w2[0]
        if self.P:
            cum = np.asarray(prep.cum, dtype=float)
            ### side="right" is bisect_right, as in permutation_test().
            idx = np.searchsorted(cum, self.rng.random(self.P) * prep.total,
                                  side="right")
            self.acc += wd[idx]
            self.W += w[idx]
            self.acc2 += wd2[idx]
            self.Wc += w2[idx]

    def finish(self, pbar):
        np = self.np
        P, ncat = self.P, self.ncat
        ncls = ncat // 16
        obs = self.obs / self.W0 if self.W0 > 0 else np.zeros(ncat)
        D = np.where(self.W[:, None] > 0,
                     self.acc / np.where(self.W > 0, self.W, 1.0)[:, None], 0.0)
        pb = np.asarray(pbar, dtype=float)
        pb6 = pb.reshape(ncls, 16).sum(axis=1)
        obs6 = obs.reshape(ncls, 16).sum(axis=1)
        D6 = D.reshape(P, ncls, 16).sum(axis=2)
        nan = float("nan")
        if not P:
            return (obs6.tolist(), [nan] * ncls, nan, nan, [],
                    {"d": obs.tolist(), "p": [nan] * ncat,
                     "shape_p": [nan] * ncls, "shape_omnibus_p": nan})

        def pval(o, null):
            return (np.sum(null >= o - 1e-15, axis=0) + 1) / (P + 1)

        mean6 = D6.mean(axis=0)
        p6 = pval(np.abs(obs6 - mean6), np.abs(D6 - mean6))
        ok6 = pb6 > 0
        om_obs = np.sum((obs6 - mean6)[ok6] ** 2 / pb6[ok6])
        om_null = np.sum((D6 - mean6)[:, ok6] ** 2 / pb6[ok6], axis=1)
        p_om = pval(om_obs, om_null)

        mean = D.mean(axis=0)
        p96 = pval(np.abs(obs - mean), np.abs(D - mean))

        ### Shape.  Within each class, the weighted mean shift in each
        ### context's SHARE of the class, hit child minus siblings, with
        ### weights from the class's own counts.  Scaled as a chi-squared --
        ### total weight times the sum of squared share shifts over the
        ### tree-wide shares -- it has the same null distribution however
        ### many mutations of the class the hit children carry.  A statistic
        ### built from the overall shifts would not: when a class genuinely
        ### rises, its contexts' overall shifts get noisier than any
        ### no-shift permutation, and a shift spread in proportion across
        ### contexts reads as concentrated.  The self-test checks this.
        q = np.where(np.repeat(pb6, 16) > 0, pb / np.repeat(
            np.where(pb6 > 0, pb6, 1.0), 16), 0.0)
        okc = q > 0

        def shares(acc2, Wc):
            Wr = np.repeat(Wc, 16, axis=-1)
            return np.where(Wr > 0, acc2 / np.where(Wr > 0, Wr, 1.0), 0.0)

        T2_obs = shares(self.obs2, self.Wc0)
        T2 = shares(self.acc2, self.Wc)
        m2 = T2.mean(axis=0)
        p_within = pval(np.abs(T2_obs - m2), np.abs(T2 - m2))

        def shape(t2, Wc):
            e2 = np.where(okc, (t2 - m2) ** 2 / np.where(okc, q, 1.0), 0.0)
            return Wc * e2.reshape(e2.shape[:-1] + (ncls, 16)).sum(axis=-1)

        s_obs = shape(T2_obs, self.Wc0)
        s_null = shape(T2, self.Wc)
        shape_p = pval(s_obs, s_null)
        shape_om_p = pval(s_obs.sum(), s_null.sum(axis=1))
        return (obs6.tolist(), p6.tolist(), float(om_obs), float(p_om),
                D6.tolist(),
                {"d": obs.tolist(), "p": p96.tolist(),
                 "within": T2_obs.tolist(), "within_p": p_within.tolist(),
                 "shape_p": shape_p.tolist(),
                 "shape_omnibus_p": float(shape_om_p),
                 "q": q.tolist()})


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def fmt(x):
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return f"{x:.6g}"
    return str(x)


def calibrate_context(args):
    """False-positive rates and power of --context's tests, by simulation.

    Each tree is a _simulate topology with 25 hit subtrees; every branch
    draws its mutations' contexts from a base mix, and branches below a hit
    from a shifted one.  Three scenarios, N trees each:

      null          shifted = base: every p < 0.05 is a false positive.
      proportional  the C>T class (CT|GA folded, CT unfolded) raised 3x, all
                    its contexts alike: the class test should fire, and the
                    shape test should not, since the class's mix of contexts
                    is unchanged.
      one context   one C>T context raised: the shape test should fire.

    For each, the share of p < 0.05 for: the shape test of every class, the
    shape test of the C>T class, the within-class test of every context, and
    the C>T class test.  These are the rates METHODS.md step 12 quotes.  The trees
    are small, so the shifts are large; the point is calibration under the
    null and the proportional shift, and that power is there at all.
    """
    try:
        import numpy as np
    except ImportError:
        print("error: --calibrate_context needs numpy", file=sys.stderr)
        return 1
    import trinucleotide as tn
    sch = tn.scheme(args.fold)
    n = sch.n
    k = sch.class_labels.index("CT|GA" if args.fold else "CT")
    lo, hi = 16 * k, 16 * k + 16

    class A(object):
        pass
    targs = A()
    targs.max_branch_length = None
    targs.top_branch = "include"
    targs.min_focal = 1
    targs.min_comparator = 5
    targs.permutations = args.permutations
    rng = random.Random(args.seed)
    base = [1.0 + (c % 7) for c in range(n)]
    tot = sum(base)
    pbar = [b / tot for b in base]
    prop = list(base)
    for c in range(lo, hi):
        prop[c] *= 3.0
    one = list(base)
    one[lo + 8] += 60.0
    scenarios = (("null", base), ("proportional 3x C>T", prop),
                 (f"one context +60 ({sch.labels[lo + 8]})", one))
    print(f"{args.calibrate_context} simulated trees per scenario, "
          f"{sch.name} context classes, {args.permutations} permutations, "
          f"seed {args.seed}")
    print(f"share of p < 0.05 for: shape (every class), shape "
          f"({sch.class_labels[k]}), within-class (every context), class "
          f"test ({sch.class_labels[k]})")
    print("scenario	shape_all	shape_CT	within_all	class_CT")
    for name, shifted in scenarios:
        shape_all, shape_k, within, cls = [], [], [], []
        for _rep in range(args.calibrate_context):
            tree, own, evs = _context_tree(rng, base, shifted)
            V, T = subtree_sums(tree, own, n)
            acc = ContextAccumulator(n, args.permutations,
                                     np.random.default_rng(
                                         rng.randrange(1 << 30)))
            res = run_test(tree, evs, V, T, own, n, targs, pbar, rng, acc)[2]
            if res is None:
                continue
            cx = res[5]
            shape_all.extend(cx["shape_p"])
            shape_k.append(cx["shape_p"][k])
            within.extend(cx["within_p"])
            cls.append(res[1][k])

        def frac(ps):
            return sum(1 for p in ps if p < 0.05) / len(ps) if ps else \
                float("nan")
        print(f"{name}\t{frac(shape_all):.3f}\t{frac(shape_k):.3f}\t"
              f"{frac(within):.3f}\t{frac(cls):.3f}")
        sys.stdout.flush()
    return 0


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    if args.calibrate_context:
        return calibrate_context(args)
    for opt in ("input_tree", "annotation", "reference", "genes"):
        if not getattr(args, opt):
            print(f"error: --{opt} is required (or pass --selftest)",
                  file=sys.stderr)
            return 2
    effect_sets = [s.strip() for s in args.effects.split(",") if s.strip()]
    bad = [s for s in effect_sets if s not in EFFECT_SETS]
    if bad:
        print(f"error: unknown effect set(s) {bad}; choose from "
              f"{list(EFFECT_SETS)}", file=sys.stderr)
        return 2

    fold_map, fold_labels = folded_index()
    if args.fold:
        cats = fold_labels
    else:
        cats, fold_map = MUT_TYPES, None
    ctx_genes = [g for g in (args.context_genes or "").split(",") if g]
    if ctx_genes:
        try:
            import numpy  # noqa: F401
        except ImportError:
            print("error: --context_genes needs numpy", file=sys.stderr)
            return 1
    elif args.mask:
        print("note: --mask only matters with --context_genes; ignoring it",
              file=sys.stderr)
    null_genes = set(g for g in (args.null_draws_genes or "").split(",") if g)

    seqs = load_fasta(args.reference)
    if args.seqid:
        ref = seqs.get(args.seqid)
        if ref is None:
            print(f"error: {args.seqid} not in {args.reference}",
                  file=sys.stderr)
            return 2
    elif len(seqs) == 1:
        ref = next(iter(seqs.values()))
    else:
        print(f"error: {args.reference} has {len(seqs)} sequences; pass "
              f"--seqid", file=sys.stderr)
        return 2
    features = load_annotation(args.annotation, args.seqid)
    index, genes = resolve_list(features, args.genes)
    missing_ctx = [g for g in ctx_genes if g not in genes]
    if missing_ctx:
        print(f"error: --context_genes {', '.join(missing_ctx)} not among the "
              f"genes tested from {args.genes}: {', '.join(genes)}",
              file=sys.stderr)
        return 2
    control = resolve_list(features, args.control) if args.control else None

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

    cap = args.max_branch_length if args.max_branch_length is not None \
        else float("inf")
    tops = ["include", "exclude"] if args.top_branch == "both" \
        else [args.top_branch]
    pre = args.prefix
    written = []

    ### The class-level tests: the list (both readings) and its control.
    ncat = len(cats)
    own = own_spectra(tree, tree.types, ncat, fold_map, cap)
    V, T = subtree_sums(tree, own, ncat)
    pbar = proportions(V)
    events = hit_events(tree, index, ref, fold_map, None)
    for top in tops:
        ### run_analysis names the reading: PREFIX_summary.tsv, or
        ### PREFIX_summary_internal.tsv for --top_branch exclude.
        written += run_analysis(
            tree, events, genes, V, T, own, ncat, cats, args, top, pbar,
            pre, null_genes)
    if control is not None:
        c_index, c_genes = control
        c_events = hit_events(tree, c_index, ref, fold_map, None)
        written += run_analysis(tree, c_events, c_genes, V, T, own, ncat,
                                cats, args, "include", pbar,
                                f"{pre}", set(), tag="_control")

    ### The context tests, for the genes named.
    if ctx_genes:
        import trinucleotide as tn
        sch = tn.scheme(args.fold)
        assert sch.class_labels == list(cats), (sch.class_labels, cats)
        masked = tn.read_mask(args.mask.split(",") if args.mask else [],
                              len(ref))
        contexts, cstats = sch.call(tree.parent, tree.children,
                                    tree.positions, tree.types, ref, masked)
        n_snv = sum(cstats[k] for k in ("called", "unknown", "doublet"))
        print(f"contexts: {cstats['called']} of {n_snv} SNVs called, "
              f"{cstats['unknown']} unknown flank, {cstats['doublet']} "
              f"doublet", file=sys.stderr)
        if cstats["ref_mismatch"] > 0.01 * max(n_snv, 1):
            print(f"  WARNING: {cstats['ref_mismatch']} SNVs disagree with the "
                  f"reconstructed parent base; the FASTA is probably not the "
                  f"tree's root genome, and every context is suspect",
                  file=sys.stderr)
        n_ctx = sch.n
        zero = [0] * n_ctx
        own_c = []
        for i, carr in enumerate(contexts):
            if not carr or tree.raw_len[i] > cap:
                own_c.append(zero)
                continue
            vec = [0] * n_ctx
            for c in carr:
                if c >= 0:
                    vec[c] += 1
            own_c.append(vec)
        V_c, T_c = subtree_sums(tree, own_c, n_ctx)
        pbar_c = proportions(V_c)
        ctx_events = [dict(e) for e in events if e["gene"] in ctx_genes]
        for e in ctx_events:
            carr = contexts[e["node"]]
            e["hit_cats"] = [carr[k] for k in e["hit_k"] if carr[k] >= 0]
        for top in tops:
            tag = "" if top == "include" else "_internal"
            written += run_analysis(
                tree, ctx_events, ctx_genes, V_c, T_c, own_c, n_ctx, cats,
                args, top, pbar_c, f"{pre}", set(), tag=f"_context{tag}",
                sch=sch)

    print(f"wrote {', '.join(written)}", file=sys.stderr)

    ### The report: the class-level summaries read together.
    main_rows = read_summary(f"{pre}_summary.tsv") \
        if "include" in tops else read_summary(f"{pre}_summary_internal.tsv")
    internal = read_summary(f"{pre}_summary_internal.tsv") \
        if len(tops) == 2 else None
    ctrl = read_summary(f"{pre}_summary_control.tsv") \
        if control is not None else None
    shim = argparse.Namespace(
        summary=f"{pre}_summary.tsv",
        internal=f"{pre}_summary_internal.tsv" if internal else None,
        control=f"{pre}_summary_control.tsv" if ctrl else None)
    summary_report(main_rows, shim, sys.stdout, internal, ctrl)
    if ctx_genes:
        context_report(pre, tops, cats, sys.stdout)
    return 0


# --------------------------------------------------------------------------
# the report (formerly summarize_gene_hits.py)
# --------------------------------------------------------------------------

NONSYN = "nonsynonymous"
SYN = "synonymous"
FEW_EFFECTIVE = 10


def read_summary(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def classes_of(rows):
    if not rows:
        return []
    return [c[len("diff_"):] for c in rows[0] if c.startswith("diff_")]


def genes(rows, es):
    """Per-gene rows (not ALL) of one effect set that have a test result."""
    return [r for r in rows if r["effect_set"] == es and r["gene"] != "ALL"
            and num(r.get("omnibus_p")) is not None]


def pooled(rows, es):
    for r in rows:
        if r["effect_set"] == es and r["gene"] == "ALL":
            return r
    return None


def permutations_from(rows):
    """Smallest N for which every p is (integer)/(N + 1), or None."""
    ps = [num(r[c]) for r in rows for c in r
          if (c.startswith("p_") or c == "omnibus_p") and num(r[c])]
    if not ps:
        return None
    for n in (100, 1000, 10000, 100000, 1000000):
        if all(abs(p * (n + 1) - round(p * (n + 1))) < 1e-3 * max(1, p * n)
               for p in ps):
            return n
    return None


def binom_tail(k, n, p):
    """P(X >= k), X ~ Binomial(n, p)."""
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i)
               for i in range(k, n + 1))


def fisher_greater(a, b, c, d):
    """One-sided Fisher p that row 1's rate a/(a+b) exceeds c/(c+d)."""
    n1, n2, k = a + b, c + d, a + c

    def h(x):
        return math.comb(n1, x) * math.comb(n2, k - x) / math.comb(n1 + n2, k)
    return sum(h(x) for x in range(a, min(n1, k) + 1))


def counts(rows, es, cls):
    g = genes(rows, es)
    n = len(g)
    o = sum(1 for r in g if num(r["omnibus_p"]) < 0.05)
    tests = [(r, c) for r in g for c in cls if num(r.get(f"p_{c}")) is not None]
    k = sum(1 for r, c in tests if num(r[f"p_{c}"]) < 0.05)
    return n, o, len(tests), k


def best_class(r, cls):
    """(class, shift in points, p) for the class with the smallest p."""
    best = None
    for c in cls:
        p, d = num(r.get(f"p_{c}")), num(r.get(f"diff_{c}"))
        if p is None or d is None:
            continue
        key = (p, -abs(d))
        if best is None or key < best[0]:
            best = (key, c, 100 * d, p)
    return best[1:] if best else (None, None, None)


def fmt_p(p):
    if p is None:
        return "NA"
    return f"{p:.1e}" if p < 0.001 else f"{p:.3g}"


def eff_str(r):
    e = num(r.get("effective_events")) if r else None
    return f"{e:,.0f}" if e is not None else "NA"


def by_gene(rows, es):
    return {r["gene"]: r for r in rows if r["effect_set"] == es}


def count_line(label, n, o, t, k, out):
    eo, ek = 0.05 * n, 0.05 * t
    print(f"  {label:<14} genes at omnibus p < 0.05: {o} of {n} "
          f"(chance {eo:.1f}, P(>= {o}) = {fmt_p(binom_tail(o, n, 0.05))});  "
          f"class tests at p < 0.05: {k} of {t} (chance {ek:.1f}, "
          f"P = {fmt_p(binom_tail(k, t, 0.05))})", file=out)


def summary_report(main, args, out, internal=None, control=None):
    cls = classes_of(main)
    n_perm = permutations_from(main)
    floor = 1 / (n_perm + 1) if n_perm else None

    print(f"GENE HIT CONTRASTS: {args.summary}", file=out)
    for es in (NONSYN, SYN):
        a = pooled(main, es)
        if a:
            print(f"  {es}: {len(genes(main, es))} genes tested, "
                  f"{int(num(a['informative']) or 0):,} informative events"
                  + (f" ({eff_str(a)} effective)"
                     if num(a.get("effective_events")) else ""), file=out)
    if n_perm:
        print(f"  {n_perm:,} permutations (read from the p-values): the "
              f"smallest possible p is {fmt_p(floor)}", file=out)

    # 1. control ------------------------------------------------------------
    print(file=out)
    print("1. IS THE COMPARISON SOUND?  Synonymous hits should look like "
          "chance.", file=out)
    n, o, t, k = counts(main, SYN, cls)
    if n:
        count_line("synonymous", n, o, t, k, out)
        bad = binom_tail(o, n, 0.05) < 0.01 or binom_tail(k, t, 0.05) < 0.01
        print("  -> " + ("FAILS: the synonymous hits show more than chance, "
                         "so the comparison is biased and the nonsynonymous "
                         "rows cannot be read either."
                         if bad else "passes: no more than chance."), file=out)
    else:
        print("  no synonymous rows; run with --effects including "
              "synonymous", file=out)

    # 2. signal -------------------------------------------------------------
    print(file=out)
    print("2. IS THERE SIGNAL?  Nonsynonymous hits against chance"
          + (" and against the control genes." if control else "."),
          file=out)
    n, o, t, k = counts(main, NONSYN, cls)
    count_line("nonsynonymous", n, o, t, k, out)
    if control:
        cn, co, ct, ck = counts(control, NONSYN, cls)
        count_line("control genes", cn, co, ct, ck, out)
        print(f"  list vs control, one-sided Fisher: genes p = "
              f"{fmt_p(fisher_greater(o, n - o, co, cn - co))}, class tests "
              f"p = {fmt_p(fisher_greater(k, t - k, ck, ct - ck))}", file=out)
        for es in (NONSYN, SYN):
            a = pooled(control, es)
            if not a:
                continue
            shifts = ", ".join(f"{c} {100 * num(a[f'diff_{c}']):+.2f}"
                               for c in cls if num(a.get(f"diff_{c}")) is not None)
            print(f"  control baseline, any {es} hit (pooled, omnibus p "
                  f"{fmt_p(num(a['omnibus_p']))}): {shifts} points", file=out)

    # 3. genes that clear the bar ------------------------------------------
    ns = genes(main, NONSYN)
    bar = 0.05 / len(ns) if ns else None
    print(file=out)
    print(f"3. GENES THAT CLEAR THE BAR: omnibus p < 0.05 / {len(ns)} = "
          f"{fmt_p(bar)}", file=out)
    if floor and bar is not None and bar < floor:
        print(f"  NOTE: the bar is below the smallest p {n_perm:,} "
              f"permutations can give ({fmt_p(floor)}), so no gene can clear "
              f"it.  Rerun with more --permutations.", file=out)
    syn = by_gene(main, SYN)
    inn = by_gene(internal, NONSYN) if internal else {}
    clear = sorted((r for r in ns if num(r["omnibus_p"]) < bar),
                   key=lambda r: num(r["omnibus_p"]))
    if not clear:
        print("  none", file=out)
    for r in clear:
        at_floor = floor and num(r["omnibus_p"]) <= floor * 1.0001
        print(f"  {r['gene']}: omnibus p {fmt_p(num(r['omnibus_p']))}"
              + (" (the permutation floor; the true p may be smaller)"
                 if at_floor else ""), file=out)
        shown = sorted((c for c in cls if num(r.get(f"p_{c}")) is not None
                        and num(r[f"p_{c}"]) < 0.05),
                       key=lambda c: num(r[f"p_{c}"]))
        if not shown:
            shown = [best_class(r, cls)[0]]
        ri, rs = inn.get(r["gene"]), syn.get(r["gene"])
        cols = [("all hits", r)] + ([("internal hits", ri)] if ri else []) \
            + ([("synonymous", rs)] if rs else [])
        table = [["  class"] + [name for name, _x in cols]]
        for c in shown:
            line = [f"  {c}"]
            for _name, x in cols:
                d, p = num(x.get(f"diff_{c}")), num(x.get(f"p_{c}"))
                line.append(f"{100 * d:+.2f} pts, p {fmt_p(p)}"
                            if d is not None else "NA")
            table.append(line)
        table.append(["  events"] + [
            f"{int(num(x['informative']) or 0):,} ({eff_str(x)} effective)"
            for _name, x in cols])
        widths = [max(len(t[i]) for t in table) for i in range(len(table[0]))]
        for t in table:
            print("  " + "   ".join(v.ljust(w) for v, w in zip(t, widths))
                  .rstrip(), file=out)
        e = num(rs.get("effective_events")) if rs else None
        if e is not None and e < FEW_EFFECTIVE:
            print(f"    its synonymous control has only {e:.0f} effective "
                  f"events: a weak control; the one across genes (section 1) "
                  f"carries the weight", file=out)
        e = num(r.get("effective_events"))
        if e is not None and e < FEW_EFFECTIVE:
            print(f"    WARNING: only {e:.0f} effective events; the result "
                  f"is about those few branches", file=out)

    # 4. leads --------------------------------------------------------------
    leads = sorted((r for r in ns if bar <= num(r["omnibus_p"]) < 0.05),
                   key=lambda r: num(r["omnibus_p"]))
    print(file=out)
    print(f"4. LEADS: omnibus p < 0.05 but above the bar ({len(leads)}; "
          f"chance alone gives about {0.05 * len(ns):.1f})", file=out)
    for r in leads:
        c, d, p = best_class(r, cls)
        extra = ""
        ri = inn.get(r["gene"])
        if ri and num(ri.get(f"diff_{c}")) is not None:
            extra = (f"; internal {100 * num(ri[f'diff_{c}']):+.2f}, p "
                     f"{fmt_p(num(ri[f'p_{c}']))}")
        print(f"  {r['gene']:<10} omnibus p {fmt_p(num(r['omnibus_p'])):<8} "
              f"{c} {d:+.2f} points (p {fmt_p(p)}), {eff_str(r)} effective"
              f"{extra}", file=out)


def proportions(V):
    tot = sum(V[0])
    return [v / tot if tot else 0.0 for v in V[0]]


def resolve_list(features, path):
    """(CDS index, genes tested) for a gene list, with what did not resolve
    reported on stderr."""
    names = read_gene_list(path)
    index, genes, missing, aliases, no_cds = cds_index_for(features, names)
    print(f"{path}: {len(genes)} of {len(names)} requested genes resolve to "
          f"a CDS", file=sys.stderr)
    if missing:
        print(f"  NOT FOUND in annotation: {', '.join(missing)}",
              file=sys.stderr)
    if no_cds:
        print(f"  found but not protein-coding (not tested): "
              f"{', '.join(no_cds)}", file=sys.stderr)
    for a, b in aliases:
        print(f"  {a} is the same CDS as {b}; counted once, as {b}",
              file=sys.stderr)
    return index, genes


def hit_events(tree, index, ref, fold_map, contexts):
    """Every branch-gene hit event, with each hit's class (or context)."""
    events, ref_mismatch = find_hit_events(tree, index, ref)
    print(f"{len(events)} branch-gene hit events", file=sys.stderr)
    if ref_mismatch:
        n_hit = sum(len(e["hit_types"]) for e in events)
        print(f"  WARNING: {ref_mismatch} of {n_hit} hit mutations have a "
              f"parent base that disagrees with the reference reconstructed "
              f"up the tree.  A few are expected from reversions; many mean "
              f"the FASTA is not the tree's reference.", file=sys.stderr)
    by_effect = defaultdict(int)
    for e in events:
        by_effect[e["effect"]] += 1
    print("  by effect: " + ", ".join(f"{k} {by_effect[k]}" for k in SEVERITY
                                       if by_effect[k]), file=sys.stderr)
    for e in events:
        e["hit_cats"] = [fold_map[t] if fold_map else t
                         for t in e["hit_types"]]
    return events


def run_analysis(tree, events, genes, V, T, own, ncat, cats, args, top, pbar,
                 prefix, null_genes, tag=None, sch=None):
    """One reading: every effect set, every gene and ALL, under one
    --top_branch, written to PREFIX_summary{tag}.tsv and PREFIX_events{tag}.tsv
    (tag "" or "_internal" by top, unless given), plus the null draws for
    null_genes and, with sch, the per-context table.  The random streams
    start from --seed, so each reading reproduces on its own.  Returns the
    paths written."""
    import copy
    a = copy.copy(args)
    a.top_branch = top
    if tag is None:
        tag = "" if top == "include" else "_internal"
    context = sch is not None
    nsum = len(cats)
    effect_sets = [x.strip() for x in args.effects.split(",") if x.strip()]

    def summary_vec(v):
        if not context:
            return v
        return [sum(v[16 * k:16 * k + 16]) for k in range(nsum)]

    rng = random.Random(args.seed)
    ev_cols = (["effect_set", "gene", "node_id", "parent_id", "effect",
                "aa_changes", "mutations", "terminal", "nested_under",
                "parent_children", "siblings_used", "informative",
                "focal_tips", "focal_mutations", "comparator_tips",
                "comparator_mutations", "chi2"]
               + [f"focal_{c}" for c in cats]
               + [f"comparator_{c}" for c in cats])
    sum_cols = (["effect_set", "gene", "events", "informative",
                 "informative_terminal", "effective_events"]
                + [f"diff_{c}" for c in cats] + [f"p_{c}" for c in cats]
                + ["omnibus", "omnibus_p"])
    paths = [f"{prefix}_summary{tag}.tsv", f"{prefix}_events{tag}.tsv"]
    ctx_out = null_out = None
    if context:
        import numpy
        sum_cols += [f"shape_p_{c}" for c in cats] + ["shape_omnibus_p"]
        np_rng = numpy.random.default_rng(args.seed)
        ### The context table is named for the reading, not the summary:
        ### PREFIX_context.tsv and PREFIX_context_internal.tsv.
        ctx_path = f"{prefix}{tag}.tsv"
        paths.append(ctx_path)
        ctx_out = open(ctx_path, "w")
        ctx_out.write("effect_set\tgene\tclass\tcontext\tclass_share\t"
                      "diff\tproportional_diff\texcess\tp\t"
                      "within_class_diff\twithin_class_p\n")
    if null_genes:
        paths.append(f"{prefix}_null_draws{tag}.tsv")
        null_out = open(paths[-1], "w")
        null_out.write("\t".join(["effect_set", "gene", "permutation"]
                                 + [f"diff_{c}" for c in cats]) + "\n")
    print(f"testing {prefix}{tag or ''} ({top} top branch) ...",
          file=sys.stderr)
    with open(paths[1], "w") as ev_out, open(paths[0], "w") as sum_out:
        ev_out.write("\t".join(ev_cols) + "\n")
        sum_out.write("\t".join(sum_cols) + "\n")
        for es in effect_sets:
            in_set = [e for e in events if e["effect"] in EFFECT_SETS[es]]
            for gene in genes + ["ALL"]:
                test_events = in_set if gene == "ALL" else \
                    [e for e in in_set if e["gene"] == gene]
                accumulator = ContextAccumulator(
                    ncat, args.permutations, np_rng) if context else None
                rows, n_inf, result, n_eff = run_test(
                    tree, test_events, V, T, own, ncat, a, pbar, rng,
                    accumulator)
                if gene != "ALL":
                    for r in rows:
                        fv = summary_vec(r["focal_vec"])
                        cv = summary_vec(r["comparator_vec"])
                        nf, nc = sum(fv), sum(cv)
                        out = [es] + [r[c] for c in ev_cols[1:17]]
                        out += [fv[k] / nf if nf else float("nan")
                                for k in range(nsum)]
                        out += [cv[k] / nc if nc else float("nan")
                                for k in range(nsum)]
                        ev_out.write("\t".join(fmt(x) for x in out) + "\n")
                n_inf_term = sum(1 for r in rows
                                 if r["informative"] and r["terminal"])
                out = [es, gene, len(rows), n_inf, n_inf_term,
                       round(n_eff, 1)]
                if result is None:
                    out += [float("nan")] * (2 * nsum + 2)
                    if context:
                        out += [float("nan")] * (nsum + 1)
                else:
                    d, pc, om, pom, perm_d = result[:5]
                    out += d + pc + [om, pom]
                    if context:
                        cx = result[5]
                        out += cx["shape_p"] + [cx["shape_omnibus_p"]]
                        for c in range(ncat):
                            k = c // 16
                            prop = d[k] * cx["q"][c]
                            ctx_out.write("\t".join(
                                [es, gene, cats[k], sch.labels[c]]
                                + [fmt(x) for x in (
                                    cx["q"][c], cx["d"][c], prop,
                                    cx["d"][c] - prop, cx["p"][c],
                                    cx["within"][c], cx["within_p"][c])])
                                + "\n")
                    if null_out and gene in null_genes:
                        for i, pd in enumerate(perm_d, 1):
                            null_out.write("\t".join(
                                [es, gene, str(i)] + [fmt(x) for x in pd])
                                + "\n")
                sum_out.write("\t".join(fmt(x) for x in out) + "\n")
                if result is not None and gene == "ALL":
                    print(f"  {es} ALL: {n_inf} informative events, omnibus "
                          f"p {fmt(result[3])}", file=sys.stderr)
    for fh in (ctx_out, null_out):
        if fh:
            fh.close()
    return paths


def context_report(pre, tops, cats, out):
    """The context tests' shape p values, per gene, per reading."""
    print("=" * 72, file=out)
    print("TRINUCLEOTIDE CONTEXT (--context_genes)", file=out)
    print("=" * 72, file=out)
    print("  shape p: does a class's shift fall on its 16 contexts in "
          "proportion to how", file=out)
    print("  common each is (large p), or in particular ones (small p)?  "
          "Read the omnibus", file=out)
    print(f"  first; which contexts is in {pre}_context.tsv.", file=out)
    for top in tops:
        tag = "" if top == "include" else "_internal"
        rows = read_summary(f"{pre}_summary_context{tag}.tsv")
        print(f"\n  {'all hits' if top == 'include' else 'internal hits'} "
              f"({pre}_summary_context{tag}.tsv):", file=out)
        print("    " + "\t".join(["effect", "gene", "events", "effective",
                                  "omnibus_p", "shape_omnibus_p"]
                                 + [f"shape_{c}" for c in cats]), file=out)
        for r in rows:
            print("    " + "\t".join(
                [r["effect_set"], r["gene"], r["informative"],
                 r["effective_events"], fmt_p(num(r["omnibus_p"])),
                 fmt_p(num(r["shape_omnibus_p"]))]
                + [fmt_p(num(r[f"shape_p_{c}"])) for c in cats]), file=out)
    print(file=out)


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

class _N(object):
    def __init__(self, i, muts, kids=()):
        self.id = i
        self.mutations = list(muts)
        self.children = list(kids)


def _brute_child(tree, c, hits, own, top, ncat):
    """Counts and tips for c by explicit traversal, skipping hit subtrees."""
    v = [0] * ncat
    tips = 0
    stack = [c]
    while stack:
        x = stack.pop()
        if x != c and x in hits:
            continue
        mine = own[x] if (x != c or top == "include") else [0] * ncat
        for k in range(ncat):
            v[k] += mine[k]
        kids = tree.children[x]
        if not kids:
            tips += 1
        stack.extend(kids)
    return v, tips


def _simulate(rng, n_parents, effect, ncat):
    """Star-of-polytomies tree with one hit child per parent.

    Each parent has 2-40 children, each a small subtree.  Every child draws
    its mutations from the same base spectrum, except the hit child's
    subtree, which has `effect` added to class 0's weight.
    """
    base = [1.0] * 12
    counter = [0]

    def muts(n, weights):
        out = []
        for _ in range(n):
            t = rng.choices(range(12), weights=weights)[0]
            counter[0] += 1
            mt = MUT_TYPES[t]
            out.append(f"{mt[0]}{counter[0]}{mt[1]}")
        return out

    def subtree(prefix, weights, top_len):
        kids = []
        for j in range(rng.randint(0, 3)):
            kids.append(_N(f"{prefix}.{j}", muts(rng.randint(0, 6), weights)))
        return _N(prefix, muts(top_len, weights), kids)

    shifted = list(base)
    shifted[0] += effect
    parents = []
    hit_ids = []
    for i in range(n_parents):
        kids = []
        n_kids = rng.randint(2, 40)
        h = rng.randrange(n_kids)
        for j in range(n_kids):
            if j == h:
                node = subtree(f"p{i}.c{j}", shifted, rng.randint(1, 6))
                hit_ids.append(node.id)
            else:
                node = subtree(f"p{i}.c{j}", base, rng.randint(0, 5))
            kids.append(node)
        parents.append(_N(f"p{i}", muts(2, base), kids))
    return _N("root", [], parents), hit_ids


def selftest():
    fold_map, labels = folded_index()

    # ---- codon effects, both strands ----------------------------------
    from annotate_sites import Feature
    ### Forward gene 1..9: ATG GCT TAA.  Reverse gene 11..19 whose codons,
    ### read on the minus strand from 19 down, are ATG TGG TAA.
    ref = "ATGGCTTAA" + "N" + "TTACCACAT"
    fwd = Feature("s", "CDS", 1, 9, 1, "fwd", "F1", None)
    rev = Feature("s", "CDS", 11, 19, -1, "rev", "R1", None)
    index = defaultdict(list)
    for f in (fwd, rev):
        for p in range(f.start, f.end + 1):
            index[p].append((f.gene, f))
    root = _N("root", [], [
        _N("a", ["C5A"]),            # GCT -> GAT  A2D missense
        _N("b", ["T6C"]),            # GCT -> GCC  synonymous
        _N("c", ["A1G"]),            # ATG -> GTG  start kept
        _N("d", ["G3C"]),            # ATG -> ATC  start_lost
        _N("e", ["T7C"]),            # TAA -> CAA  stop_lost
        _N("f", ["G4T"], [           # GCT -> TCT  A2S on an ancestor ...
            _N("g", ["C5A"])]),      # ... so TCT -> TAT: S2Y, not A2D
        _N("h", ["C15T"]),           # rev codon 2 TGG: pos 15 is codon pos 2,
                                     # G->A on minus: TAG nonsense W2*
        _N("i", ["C5A", "T6G"]),     # two in one codon: GCT -> GAG  A2E
    ])
    tree = extract(root)
    events, mism = find_hit_events(tree, index, ref)
    got = {tree.node_ids[e["node"]]: (e["effect"], e["aa_changes"])
           for e in events}
    want = {"a": ("missense", "A2D"), "b": ("synonymous", "A2A"),
            "c": ("synonymous", "M1M"), "d": ("start_lost", "M1I"),
            "e": ("stop_lost", "*3Q"), "f": ("missense", "A2S"),
            "g": ("missense", "S2Y"), "h": ("nonsense", "W2*"),
            "i": ("missense", "A2E")}
    for k, v in want.items():
        assert got.get(k) == v, (k, got.get(k), v)
    assert mism == 0, mism
    print("ok: codon effects on both strands, with ancestral codon state and "
          "multi-hit codons")

    # ---- subtree bookkeeping against brute force ---------------------
    ### P has children: h (hit, with a nested hit h2 below), s1 (internal),
    ### s2, s3 (tips), x (another hit at the same parent).
    h2 = _N("h2", ["A10G"], [_N("h2a", ["C11T"]), _N("h2b", ["C12A"])])
    h = _N("h", ["C1T", "G2A"], [h2, _N("ha", ["T3C"]),
                                 _N("hb", ["A4G"], [_N("hb1", ["C13T"]),
                                                    _N("hb2", [])])])
    s1 = _N("s1", ["C5T"], [_N("s1a", ["C6T", "A7G"]), _N("s1b", ["G8T"])])
    P = _N("P", ["A20T"], [h, s1, _N("s2", ["C9A"]), _N("s3", []),
                           _N("x", ["T14C"], [_N("xa", ["C15T"])])])
    root = _N("root", ["G21A"], [P, _N("other", ["C22T"])])
    tree = extract(root)
    idx = {nid: i for i, nid in enumerate(tree.node_ids)}
    for ncat, fm in ((12, None), (6, fold_map)):
        own = own_spectra(tree, tree.types, ncat, fm, float("inf"))
        V, T = subtree_sums(tree, own, ncat)
        hits = {idx["h"], idx["h2"], idx["x"]}
        adj, nested = cut_adjustments(tree, V, T, hits, ncat)
        assert tree.node_ids[nested[idx["h2"]]] == "h"
        assert idx["h"] not in nested and idx["x"] not in nested
        for top in ("include", "exclude"):
            for name in ("h", "s1", "s2", "s3", "x", "h2", "P"):
                c = idx[name]
                got = child_value(c, V, T, own, adj, top, ncat)
                want = _brute_child(tree, c, hits, own, top, ncat)
                assert tuple(got[0]) == tuple(want[0]) and \
                    got[1] == want[1], (ncat, top, name, got, want)
    print("ok: cut subtree counts and tips match brute force, "
          "including a nested hit, with and without the top branch")

    ### The event's comparator leaves out x, the other hit child of P.
    own = own_spectra(tree, tree.types, 12, None, float("inf"))
    V, T = subtree_sums(tree, own, 12)
    hits = {idx["h"], idx["x"]}
    adj, _ = cut_adjustments(tree, V, T, hits, 12)
    row, prep = prepare_event(tree, idx["h"], hits, V, T, own, adj,
                              "include", 12, [0] * 12, 1, 1)
    ### siblings s1 (4 muts), s2 (1), s3 (0): comparator = 5 mutations.
    assert row["siblings_used"] == 3 and row["comparator_mutations"] == 5, row
    ### eligible: h, s1, s2 (s3 has none); weights 2+1, 1+1, 1+1.
    assert len(prep.cum) == 3 and prep.cum == [3, 5, 7], prep.cum
    print("ok: comparator pools all non-hit siblings; draws are size-biased")

    # ---- the null is calibrated, and the test has power --------------
    class A(object):
        pass
    targs = A()
    targs.max_branch_length = None
    targs.top_branch = "include"
    targs.min_focal = 1
    targs.min_comparator = 5
    targs.permutations = 200
    rng = random.Random(7)
    pvals = []
    for rep in range(40):
        root, hit_ids = _simulate(rng, 25, 0.0, 12)
        tree = extract(root)
        idx = {nid: i for i, nid in enumerate(tree.node_ids)}
        own = own_spectra(tree, tree.types, 12, None, float("inf"))
        V, T = subtree_sums(tree, own, 12)
        evs = [{"node": idx[x], "gene": "g", "effect": "missense",
                "aa_changes": "", "mutations": "", "hit_types": [],
                "hit_cats": []}
               for x in hit_ids]
        pbar = [1 / 12.0] * 12
        _rows, n_inf, res, _n_eff = run_test(tree, evs, V, T, own, 12,
                                     targs, pbar, rng)
        pvals.append((res[1][0], res[3]))
    frac_class = sum(1 for p, _ in pvals if p < 0.05) / len(pvals)
    frac_omni = sum(1 for _, p in pvals if p < 0.05) / len(pvals)
    assert frac_class <= 0.2 and frac_omni <= 0.2, (frac_class, frac_omni)
    print(f"ok: under the null, p < 0.05 in {frac_class:.0%} (class 0) and "
          f"{frac_omni:.0%} (omnibus) of 40 simulated trees")

    hits_p = []
    for rep in range(10):
        root, hit_ids = _simulate(rng, 25, 3.0, 12)
        tree = extract(root)
        idx = {nid: i for i, nid in enumerate(tree.node_ids)}
        own = own_spectra(tree, tree.types, 12, None, float("inf"))
        V, T = subtree_sums(tree, own, 12)
        evs = [{"node": idx[x], "gene": "g", "effect": "missense",
                "aa_changes": "", "mutations": "", "hit_types": [],
                "hit_cats": []}
               for x in hit_ids]
        _rows, n_inf, res, _n_eff = run_test(tree, evs, V, T, own, 12,
                                     targs, [1 / 12.0] * 12, rng)
        assert res[0][0] > 0, res[0]
        hits_p.append(res[1][0])
    power = sum(1 for p in hits_p if p < 0.05) / len(hits_p)
    assert power >= 0.8, hits_p
    print(f"ok: a 4x shift in one class on 25 hit subtrees is detected in "
          f"{power:.0%} of 10 simulated trees, in the right direction")

    if _selftest_summary():
        return 1
    if _selftest_run():
        return 1

    try:
        import numpy
    except ImportError:
        print("skipped: --context checks need numpy")
        print("all self-tests passed")
        return 0
    _selftest_context(numpy, targs)
    print("all self-tests passed")
    return 0


def _context_tree(rng, base, shifted):
    """A _simulate topology with context-class vectors (len(base) of them)
    drawn per branch.

    Branches inside a hit subtree draw contexts from `shifted`, all others
    from `base`.
    """
    root, hit_ids = _simulate(rng, 25, 0.0, 12)
    tree = extract(root)
    idx = {nid: i for i, nid in enumerate(tree.node_ids)}
    hit_set = set(idx[x] for x in hit_ids)
    inside = [False] * len(tree.node_ids)
    own = []
    for i in range(len(tree.node_ids)):
        p = tree.parent[i]
        inside[i] = i in hit_set or (p >= 0 and inside[p])
        vec = [0] * len(base)
        for c in rng.choices(range(len(base)),
                             weights=shifted if inside[i] else base,
                             k=len(tree.types[i])):
            vec[c] += 1
        own.append(vec)
    evs = [{"node": idx[x], "gene": "g", "effect": "missense",
            "aa_changes": "", "mutations": "", "hit_types": [],
            "hit_cats": []} for x in hit_ids]
    return tree, own, evs


def _selftest_context(np, targs):
    rng = random.Random(11)
    base = [1.0 + (c % 7) for c in range(96)]
    tot = sum(base)
    pbar = [b / tot for b in base]

    def run(shifted, permutations):
        tree, own, evs = _context_tree(rng, base, shifted)
        V, T = subtree_sums(tree, own, 96)
        acc = ContextAccumulator(96, permutations, np.random.default_rng(
            rng.randrange(1 << 30)))
        _r, _n, res, _e = run_test(tree, evs, V, T, own, 96, targs, pbar,
                                   rng, acc)
        return tree, own, evs, res

    ### Summing contexts gives exactly the 6-class statistic for the same
    ### events: the weights depend only on mutation totals.
    shifted = list(base)
    shifted[40] += 2.0
    tree, own, evs, res = run(shifted, 50)
    own6 = [[sum(v[16 * k:16 * k + 16]) for k in range(6)] for v in own]
    V6, T6 = subtree_sums(tree, own6, 6)
    targs0 = type(targs)()
    targs0.__dict__.update(targs.__dict__)
    targs0.permutations = 0
    pbar6 = [sum(pbar[16 * k:16 * k + 16]) for k in range(6)]
    _r, _n, res6, _e = run_test(tree, evs, V6, T6, own6, 6, targs0, pbar6, rng)
    assert all(abs(a - b) < 1e-12 for a, b in zip(res[0], res6[0])), \
        (res[0], res6[0])
    assert abs(sum(res[5]["d"][32:48]) - res[0][2]) < 1e-12
    print("ok: --context class shifts are exactly the 6-class shifts")

    ### The same with 192 forward-strand classes, against the 12 types.
    base192 = [1.0 + (c % 5) for c in range(192)]
    t192 = sum(base192)
    pbar192 = [b / t192 for b in base192]
    sh192 = list(base192)
    sh192[5 * 16 + 3] += 2.0
    tree, own, evs = _context_tree(rng, base192, sh192)
    V, T = subtree_sums(tree, own, 192)
    acc = ContextAccumulator(192, 50, np.random.default_rng(3))
    _r, _n, res, _e = run_test(tree, evs, V, T, own, 192, targs, pbar192,
                               rng, acc)
    own12 = [[sum(v[16 * k:16 * k + 16]) for k in range(12)] for v in own]
    V12, T12 = subtree_sums(tree, own12, 12)
    pbar12 = [sum(pbar192[16 * k:16 * k + 16]) for k in range(12)]
    _r, _n, res12, _e = run_test(tree, evs, V12, T12, own12, 12, targs0,
                                 pbar12, rng)
    assert all(abs(a - b) < 1e-12 for a, b in zip(res[0], res12[0]))
    assert len(res[5]["shape_p"]) == 12
    print("ok: with 192 forward-strand classes, the class shifts are exactly "
          "the 12-type shifts, with a shape p per type")

    ### Shape test, on C>T (contexts 32..47, class 2).  No effect; the whole
    ### class raised in proportion to its contexts, which should move the
    ### class but not the shape; and one context raised.
    def ps(shifted, n):
        out = []
        for _ in range(n):
            res = run(shifted, 200)[3]
            out.append((res[5]["shape_p"][2], res[1][2]))
        return out

    prop_shift = list(base)
    for c in range(32, 48):
        prop_shift[c] *= 3.0
    one = list(base)
    one[40] += 60.0
    null, prop, conc = ps(base, 30), ps(prop_shift, 30), ps(one, 10)

    def frac(rows, j):
        return sum(1 for r in rows if r[j] < 0.05) / len(rows)

    assert frac(null, 0) <= 0.2, null
    assert frac(prop, 0) <= 0.2 and frac(prop, 1) >= 0.8, prop
    assert frac(conc, 0) >= 0.8, conc
    print(f"ok: C>T shape p < 0.05 in {frac(null, 0):.0%} of null trees and "
          f"{frac(prop, 0):.0%} of trees with a proportional C>T shift (whose "
          f"class test fires in {frac(prop, 1):.0%}); a one-context shift is "
          f"found in {frac(conc, 0):.0%}")


def _selftest_run():
    """run_analysis on a simulated tree: both readings write their files,
    null draws go only to the genes named, and a reading's results do not
    depend on what ran before it."""
    import tempfile
    rng = random.Random(5)
    root, hit_ids = _simulate(rng, 25, 3.0, 12)
    tree = extract(root)
    idx = {nid: i for i, nid in enumerate(tree.node_ids)}
    own = own_spectra(tree, tree.types, 12, None, float("inf"))
    V, T = subtree_sums(tree, own, 12)
    events = [{"node": idx[x], "gene": "g" if i % 2 else "h",
               "effect": "missense", "aa_changes": "", "mutations": "",
               "hit_types": [], "hit_cats": [], "hit_k": []}
              for i, x in enumerate(hit_ids)]
    args = parse_args(["--permutations", "50", "--effects", "nonsynonymous"])
    d = tempfile.mkdtemp()
    pre = os.path.join(d, "t")
    first = run_analysis(tree, events, ["g", "h"], V, T, own, 12, MUT_TYPES,
                         args, "include", proportions(V), pre, {"g"})
    internal = run_analysis(tree, events, ["g", "h"], V, T, own, 12,
                            MUT_TYPES, args, "exclude", proportions(V), pre,
                            set())
    with open(first[0]) as fh:
        a = fh.read()
    again = run_analysis(tree, events, ["g", "h"], V, T, own, 12, MUT_TYPES,
                         args, "include", proportions(V), pre + "2", {"g"})
    with open(again[0]) as fh:
        b = fh.read()
    ok = (first[0].endswith("t_summary.tsv")
          and internal[0].endswith("t_summary_internal.tsv")
          and first[2].endswith("t_null_draws.tsv") and a == b
          and all(os.path.getsize(x) > 0 for x in first + internal))
    with open(first[2]) as fh:
        genes_drawn = {line.split("\t")[1] for line in list(fh)[1:]}
    ok = ok and genes_drawn == {"g"}
    print(("ok" if ok else "FAIL") + ": both readings write their files, "
          "null draws only for the genes named, and a reading reproduces "
          "whatever ran before it")
    return 0 if ok else 1


def _selftest_summary():
    import io
    import random
    ok = True

    def check(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    rng = random.Random(3)
    cls = ["A", "B"]
    n_perm = 1000

    def p_of(k):
        return f"{(k + 1) / (n_perm + 1):.6g}"

    def row(es, gene, omni_k, a_k, a_diff, eff=50.0, inf=100):
        return {"effect_set": es, "gene": gene, "events": str(inf),
                "informative": str(inf), "informative_terminal": "0",
                "effective_events": str(eff),
                "diff_A": str(a_diff), "diff_B": "0.0",
                "p_A": p_of(a_k), "p_B": p_of(500),
                "omnibus": "1", "omnibus_p": p_of(omni_k)}

    rows = []
    for i in range(40):
        g = f"g{i}"
        if g == "g0":
            rows.append(row(NONSYN, g, 0, 0, 0.01))       # clears the bar
        elif g == "g1":
            rows.append(row(NONSYN, g, 20, 20, 0.02))     # a lead
        else:
            k = rng.randrange(100, 1000)
            rows.append(row(NONSYN, g, k, k, 0.0))
        rows.append(row(SYN, g, rng.randrange(100, 1000), 600, 0.0,
                        eff=2.0 if g == "g0" else 40.0))
    rows.append(row(NONSYN, "ALL", 500, 500, 0.001, eff=900, inf=4000))
    rows.append(row(SYN, "ALL", 500, 500, 0.0, eff=800, inf=3000))
    check(permutations_from(rows) == 1000, "permutations read from p-values")
    args = argparse.Namespace(summary="test.tsv", internal=None,
                              control=None)
    buf = io.StringIO()
    summary_report(rows, args, buf)
    text = buf.getvalue()
    check("passes: no more than chance" in text, "synonymous control passes")
    check("g0: omnibus p" in text and "permutation floor" in text,
          "g0 clears the bar, at the floor")
    check("a weak control" in text, "g0's synonymous control flagged weak")
    check("\n  g1 " in text and "LEADS" in text, "g1 listed as a lead")
    check("bar is below the smallest p" not in text,
          "0.05 / 40 is within reach of 1,000 permutations")
    more = rows + [row(NONSYN, f"h{i}", 600, 600, 0.0) for i in range(20)]
    buf = io.StringIO()
    summary_report(more, args, buf)
    text = buf.getvalue()
    check("bar is below the smallest p" in text and "g0: omnibus" not in text,
          "0.05 / 60 is not, and nothing clears it")
    check(abs(binom_tail(0, 10, 0.05) - 1.0) < 1e-12
          and abs(fisher_greater(5, 5, 0, 10) - 0.01625) < 1e-4,
          "binomial and Fisher tails")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1



if __name__ == "__main__":
    sys.exit(main())
