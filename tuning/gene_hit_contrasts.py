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

What it cannot see.  The tree holds SNVs only, so frameshifts, deletions and
IS6110 insertions -- the commonest ways a gene is actually lost -- are
invisible here.  Promoter changes are not tested.  And the ALL row pools
genes whose mechanisms differ: a gene that raises CT|GA and one that raises
CA|GT cancel there.  ALL is a screen.  The evidence is in the per-gene rows.

Usage:
    gene_hit_contrasts.py --input_tree tree.pb.gz \\
        --annotation NC_000962.3.gff3.gz --reference NC_000962.3.fna \\
        --genes tb_repair_genes.txt --fold --max_branch_length 250 \\
        --output_events hit_events.tsv --output_summary hit_summary.tsv

    gene_hit_contrasts.py --selftest
"""

import argparse
import bisect
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


def parse_args():
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
    p.add_argument("--top_branch", choices=("include", "exclude"),
                   default="include",
                   help="Count each child's own branch (default; needed for "
                        "terminal hits) or only what lies below it")
    p.add_argument("--min_focal", type=int, default=1,
                   help="Minimum mutations in the hit child's counted set for "
                        "the event to be used; siblings below this are not "
                        "drawn as pseudo-hits (default 1)")
    p.add_argument("--min_comparator", type=int, default=10,
                   help="Minimum mutations in the pooled siblings for the "
                        "event to be used (default 10)")
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--output_events", type=str, default="hit_events.tsv")
    p.add_argument("--output_summary", type=str, default="hit_summary.tsv")
    p.add_argument("--null_draws", type=str, default=None,
                   help="Write every permutation's class shifts for the "
                        "--null_draws_genes to this TSV, e.g. for "
                        "plot_gene_hits.py")
    p.add_argument("--null_draws_genes", type=str, default=None,
                   help="Comma-separated genes (or ALL) whose permutation "
                        "draws --null_draws writes")
    p.add_argument("--selftest", action="store_true",
                   help="Check subtree bookkeeping, codon effects and the "
                        "permutation null on synthetic trees, then exit")
    return p.parse_args()


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
                    by_key[(i, name)].append((pos, t, f))

    events = []
    ref_mismatch = 0
    for (node, gene), muts in by_key.items():
        parent = tree.parent[node]
        codons = defaultdict(dict)
        labels = []
        for pos, t, f in muts:
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
            "hit_types": [t for _p, t, _f in muts],
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
    adj[node] = [vec, tips] and nested[hit] = enclosing hit.
    """
    adj = {}
    nested = {}
    for h in hits:
        x = tree.parent[h]
        while x >= 0:
            a = adj.get(x)
            if a is None:
                a = adj[x] = [[0] * ncat, 0]
            vh = V[h]
            for k in range(ncat):
                a[0][k] += vh[k]
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
        for k in range(ncat):
            v[k] -= a[0][k]
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
    __slots__ = ("wd", "w", "cum", "total")


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
    acc = 0
    for v, wt in elig:
        rest = [tot_v[k] - v[k] for k in range(ncat)]
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


def run_test(tree, test_events, V, T, own, ncat, fold_map, args,
             pbar, rng):
    """Events for one test (they define the cut set), then the permutation."""
    hits = set(e["node"] for e in test_events)
    adj, nested = cut_adjustments(tree, V, T, hits, ncat)
    cap = args.max_branch_length
    rows = []
    preps = []
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
            for t in e["hit_types"]:
                hit_vec[fold_map[t] if fold_map else t] += 1
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
            preps.append(prep)
    result = None
    if preps:
        result = permutation_test(preps, ncat, pbar, args.permutations, rng)
    return rows, len(preps), result, effective_events(preps)


def effective_events(preps):
    """(sum w)^2 / sum w^2 over the real hit children's weights.

    The weights grow with subtree size, so one deep branch can carry most of
    a test.  This says how many equally weighted events the test is worth.
    """
    ws = [p.w[0] for p in preps]
    s2 = sum(w * w for w in ws)
    return sum(ws) ** 2 / s2 if s2 > 0 else 0.0


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def fmt(x):
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return f"{x:.6g}"
    return str(x)


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
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
        ncat, cats = 6, fold_labels
    else:
        ncat, cats, fold_map = 12, MUT_TYPES, None

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
    names = read_gene_list(args.genes)
    index, genes, missing, aliases, no_cds = cds_index_for(features, names)
    print(f"{len(genes)} of {len(names)} requested genes resolve to a CDS",
          file=sys.stderr)
    if missing:
        print(f"  NOT FOUND in annotation: {', '.join(missing)}",
              file=sys.stderr)
    if no_cds:
        print(f"  found but not protein-coding (not tested): "
              f"{', '.join(no_cds)}", file=sys.stderr)
    for a, b in aliases:
        print(f"  {a} is the same CDS as {b}; counted once, as {b}",
              file=sys.stderr)

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
    own = own_spectra(tree, tree.types, ncat, fold_map, cap)
    V, T = subtree_sums(tree, own, ncat)
    root_tot = sum(V[0])
    pbar = [V[0][k] / root_tot if root_tot else 0.0 for k in range(ncat)]

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

    null_genes = set(args.null_draws_genes.split(",")) \
        if args.null_draws_genes else set()
    if args.null_draws and not null_genes:
        print("error: --null_draws needs --null_draws_genes", file=sys.stderr)
        return 2
    null_out = open(args.null_draws, "w") if args.null_draws else None
    if null_out:
        null_out.write("\t".join(["effect_set", "gene", "permutation"]
                                 + [f"diff_{c}" for c in cats]) + "\n")
    with open(args.output_events, "w") as ev_out, \
            open(args.output_summary, "w") as sum_out:
        ev_out.write("\t".join(ev_cols) + "\n")
        sum_out.write("\t".join(sum_cols) + "\n")
        for es in effect_sets:
            in_set = [e for e in events if e["effect"] in EFFECT_SETS[es]]
            for gene in genes + ["ALL"]:
                test_events = in_set if gene == "ALL" else \
                    [e for e in in_set if e["gene"] == gene]
                rows, n_inf, result, n_eff = run_test(
                    tree, test_events, V, T, own, ncat, fold_map,
                    args, pbar, rng)
                if gene != "ALL":
                    for r in rows:
                        fv, cv = r["focal_vec"], r["comparator_vec"]
                        nf, nc = sum(fv), sum(cv)
                        out = [es] + [r[c] for c in ev_cols[1:17]]
                        out += [fv[k] / nf if nf else float("nan")
                                for k in range(ncat)]
                        out += [cv[k] / nc if nc else float("nan")
                                for k in range(ncat)]
                        ev_out.write("\t".join(fmt(x) for x in out) + "\n")
                n_inf_term = sum(1 for r in rows
                                 if r["informative"] and r["terminal"])
                out = [es, gene, len(rows), n_inf, n_inf_term,
                       round(n_eff, 1)]
                if result is None:
                    out += [float("nan")] * (2 * ncat + 2)
                else:
                    d, pc, om, pom, perm_d = result
                    out += d + pc + [om, pom]
                    if null_out and gene in null_genes:
                        for i, pd in enumerate(perm_d, 1):
                            null_out.write("\t".join(
                                [es, gene, str(i)] + [fmt(x) for x in pd])
                                + "\n")
                sum_out.write("\t".join(fmt(x) for x in out) + "\n")
                if result is not None and gene == "ALL":
                    print(f"{es} ALL: {n_inf} informative events, omnibus "
                          f"p {fmt(result[3])}",
                          file=sys.stderr)
    if null_out:
        null_out.close()
        print(f"wrote {args.null_draws}", file=sys.stderr)
    print(f"wrote {args.output_events} and {args.output_summary}",
          file=sys.stderr)
    return 0


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
                "aa_changes": "", "mutations": "", "hit_types": []}
               for x in hit_ids]
        pbar = [1 / 12.0] * 12
        _rows, n_inf, res, _n_eff = run_test(tree, evs, V, T, own, 12,
                                     None, targs, pbar, rng)
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
                "aa_changes": "", "mutations": "", "hit_types": []}
               for x in hit_ids]
        _rows, n_inf, res, _n_eff = run_test(tree, evs, V, T, own, 12,
                                     None, targs, [1 / 12.0] * 12, rng)
        assert res[0][0] > 0, res[0]
        hits_p.append(res[1][0])
    power = sum(1 for p in hits_p if p < 0.05) / len(hits_p)
    assert power >= 0.8, hits_p
    print(f"ok: a 4x shift in one class on 25 hit subtrees is detected in "
          f"{power:.0%} of 10 simulated trees, in the right direction")

    print("all self-tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
