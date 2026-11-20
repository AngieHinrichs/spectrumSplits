#!/usr/bin/env python3
"""Trinucleotide context for every mutation on a tree, from the parent's genome.

A shared helper, imported by profile_tree.py and gene_hit_contrasts.py, with
its own --selftest.  Nothing here is organism-specific: the reference and the
mask are inputs.

Classes.  The 96 pyrimidine-centred classes, in the order the mutational
signature literature uses (COSMIC SBS): the six strand-folded classes C>A,
C>G, C>T, T>A, T>C, T>G -- the same order as permutation_null.folded_index(),
so class k of the 6-type spectrum is contexts 16k .. 16k+15 -- and within each
class, 5' base then 3' base, each A, C, G, T.  A mutation whose reference base
is a purine is reverse-complemented, trinucleotide and all, so G>A in xGy is
counted as C>T in comp(y) C comp(x).

Folding to 96 has the same blind spot as --fold: a process tied to the leading
or transcribed strand cancels between the two replication arms.  Oriented 192
classes would be the context-level counterpart of spectrumSplits.py
--orient (parked); this module does not do that yet.

Where the flanking bases come from.  Not the reference: the PARENT's genome,
i.e. the reference with every mutation on the path from the root applied.  One
depth-first pass keeps a dictionary of the positions changed along the
current root-to-node path, applying a node's mutations on the way down and
undoing them on the way back up, so the whole tree costs one step per
mutation.  On a shallow tree this rarely changes a call, and the report says
how often it did.  The reference must be the ROOT's genome: on a rerooted tree
that is the rerooted reference, not the original assembly.  Every mutation's
own reference character is checked against the reconstructed parent base, and
the count of disagreements is reported; more than a handful means the FASTA
is not the tree's root.

What is not called (context index < 0):
  UNKNOWN (-1)  a flank is N in the reference, masked, or off the end of the
                sequence.  A masked position's state in the real genomes is
                unknown, so a context through it is unknown too.  Ends are not
                wrapped even for a circular genome; that loses at most two
                sites.
  DOUBLET (-2)  an adjacent position changed on the same branch.  CC>TT-type
                double substitutions are a signature of their own, and calling
                each half in a context that contains the other blurs both.
                Every mutation in a run of adjacent same-branch changes is
                marked.

Opportunity.  For comparisons within one tree the genome is shared and
opportunity mostly cancels.  For anything external -- reference signatures
built on another genome's composition, or a statement that one context is
"hot" -- counts must be divided by how often each trinucleotide occurs.
opportunity() counts the 32 pyrimidine-centred trinucleotides in the
reference, over centres whose trinucleotide is fully known and unmasked.
"""

import os
import sys
from array import array
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from permutation_null import MUT_TYPES

BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

### (pyrimidine ref, alt) -> class index, matching folded_index()'s order.
CLASS_ORDER = [("C", "A"), ("C", "G"), ("C", "T"),
               ("T", "A"), ("T", "C"), ("T", "G")]
CLASS_INDEX = {ra: i for i, ra in enumerate(CLASS_ORDER)}
CLASS_LABELS = ["CA|GT", "CG|GC", "CT|GA", "TA|AT", "TC|AG", "TG|AC"]

UNKNOWN = -1
DOUBLET = -2
N_CONTEXTS = 96


def context_label(idx):
    """COSMIC-style label, e.g. A[C>T]G."""
    cls, rest = divmod(idx, 16)
    left, right = divmod(rest, 4)
    ref, alt = CLASS_ORDER[cls]
    return f"{BASES[left]}[{ref}>{alt}]{BASES[right]}"


CONTEXT_LABELS = [context_label(i) for i in range(N_CONTEXTS)]


def context_index(left, ref, right, alt):
    """0..95 for a fully known context, else UNKNOWN."""
    if left not in COMPLEMENT or right not in COMPLEMENT \
            or ref not in COMPLEMENT or alt not in COMPLEMENT or ref == alt:
        return UNKNOWN
    if ref in "GA":
        left, ref, right, alt = (COMPLEMENT[right], COMPLEMENT[ref],
                                 COMPLEMENT[left], COMPLEMENT[alt])
    return (CLASS_INDEX[(ref, alt)] * 16 + BASE_INDEX[left] * 4
            + BASE_INDEX[right])


def class_of(idx):
    """The 6-type folded class a context belongs to."""
    return idx // 16


def trinuc_of_context(idx):
    """Index 0..31 of the pyrimidine-centred trinucleotide a context sits in.

    0..15 are xCy, 16..31 xTy, each by 5' then 3' base.
    """
    cls, rest = divmod(idx, 16)
    return (0 if cls < 3 else 16) + rest


TRINUC_LABELS = [f"{BASES[l]}{c}{BASES[r]}" for c in "CT"
                 for l in range(4) for r in range(4)]


def read_mask(paths, length):
    """bytearray, 1 at each masked 1-based position, from BED files."""
    masked = bytearray(length + 2)
    if not paths:
        return masked
    from annotate_sites import read_bed_intervals
    for s, e in read_bed_intervals(paths):
        s = max(s, 1)
        e = min(e, length)
        if s <= e:
            masked[s:e + 1] = b"\x01" * (e - s + 1)
    return masked


def call_contexts(parent, children, positions, types, ref, masked=None):
    """Context index for every parsed SNV on every branch.

    parent, children, positions, types: as in permutation_null.FlatTree
    (types index MUT_TYPES; parents precede children).  ref: the root genome
    as an upper-case string.  masked: optional bytearray from read_mask().

    Returns (contexts, stats): contexts[i] is an array('h') parallel to
    positions[i]; stats counts called, unknown, doublet, ref_mismatch (the
    mutation's own reference base disagrees with the reconstructed parent
    base) and flank_changed (a called context whose flank differs from the
    reference, i.e. one that reference flanks would have got wrong).
    """
    n = len(parent)
    L = len(ref)
    contexts = [None] * n
    stats = Counter()
    state = {}
    undo = {}

    def flank(p):
        if p < 1 or p > L or (masked is not None and masked[p]):
            return "N"
        b = state.get(p)
        return ref[p - 1] if b is None else b

    roots = [i for i in range(n) if parent[i] < 0]
    stack = [(r, False) for r in reversed(roots)]
    while stack:
        i, leaving = stack.pop()
        if leaving:
            for p, old in reversed(undo.pop(i)):
                if old is None:
                    del state[p]
                else:
                    state[p] = old
            continue
        pos_arr = positions[i]
        typ_arr = types[i]
        here = set(pos_arr)
        out = array("h")
        for p, t in zip(pos_arr, typ_arr):
            mt = MUT_TYPES[t]
            pb = state.get(p)
            if pb is None:
                pb = ref[p - 1] if 1 <= p <= L else "N"
            if pb != mt[0]:
                stats["ref_mismatch"] += 1
            if (p - 1) in here or (p + 1) in here:
                out.append(DOUBLET)
                stats["doublet"] += 1
                continue
            left, right = flank(p - 1), flank(p + 1)
            idx = context_index(left, mt[0], right, mt[1])
            out.append(idx)
            if idx < 0:
                stats["unknown"] += 1
                continue
            stats["called"] += 1
            if (p >= 2 and left != ref[p - 2]) or (p < L and right != ref[p]):
                stats["flank_changed"] += 1
        contexts[i] = out
        changes = []
        for p, t in zip(pos_arr, typ_arr):
            changes.append((p, state.get(p)))
            state[p] = MUT_TYPES[t][1]
        undo[i] = changes
        stack.append((i, True))
        for c in reversed(children[i]):
            stack.append((c, False))
    return contexts, stats


def opportunity(ref, masked=None):
    """Occurrences of each of the 32 pyrimidine-centred trinucleotides.

    A purine-centred trinucleotide counts as its reverse complement.  Only
    centres whose three bases are A/C/G/T and unmasked are counted.
    """
    counts = [0] * 32
    L = len(ref)
    for p in range(2, L):            # 1-based centre; needs p-1 and p+1
        if masked is not None and (masked[p - 1] or masked[p]
                                   or masked[p + 1]):
            continue
        left, c, right = ref[p - 2], ref[p - 1], ref[p]
        if left not in COMPLEMENT or c not in COMPLEMENT \
                or right not in COMPLEMENT:
            continue
        if c in "GA":
            left, c, right = COMPLEMENT[right], COMPLEMENT[c], COMPLEMENT[left]
        counts[(0 if c == "C" else 16) + BASE_INDEX[left] * 4
               + BASE_INDEX[right]] += 1
    return counts


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def _brute(tree, ref, masked):
    """Context per mutation by walking up the tree for every flank."""
    def base_at(node, p):
        if p < 1 or p > len(ref) or (masked is not None and masked[p]):
            return "N"
        x = node
        while x >= 0:
            for q, t in zip(tree.positions[x], tree.types[x]):
                if q == p:
                    return MUT_TYPES[t][1]
            x = tree.parent[x]
        return ref[p - 1]

    out = []
    for i in range(len(tree.node_ids)):
        here = set(tree.positions[i])
        row = []
        for p, t in zip(tree.positions[i], tree.types[i]):
            mt = MUT_TYPES[t]
            if (p - 1) in here or (p + 1) in here:
                row.append(DOUBLET)
                continue
            par = tree.parent[i]
            row.append(context_index(base_at(par, p - 1), mt[0],
                                     base_at(par, p + 1), mt[1]))
        out.append(row)
    return out


def selftest():
    import random
    from permutation_null import extract

    class N(object):
        def __init__(self, i, muts, kids=()):
            self.id = i
            self.mutations = list(muts)
            self.children = list(kids)

    # ---- labels and folding -----------------------------------------
    assert CONTEXT_LABELS[0] == "A[C>A]A" and CONTEXT_LABELS[95] == "T[T>G]T"
    assert CONTEXT_LABELS[2 * 16 + 0 * 4 + 2] == "A[C>T]G"
    ### A[C>T]G on the other strand is C[G>A]T.
    assert context_index("A", "C", "G", "T") == \
        context_index("C", "G", "T", "A") == 2 * 16 + 2
    ### Every one of the 192 strand-specific types lands on exactly the
    ### 96 classes, two to each.
    hits = Counter()
    for l in BASES:
        for r in BASES:
            for ref in BASES:
                for alt in BASES:
                    if ref != alt:
                        hits[context_index(l, ref, r, alt)] += 1
    assert len(hits) == 96 and set(hits.values()) == {2}, hits
    from permutation_null import folded_index
    _fold_map, labels = folded_index()
    assert labels == CLASS_LABELS, labels
    for t, mt in enumerate(MUT_TYPES):
        idx = context_index("A", mt[0], "A", mt[1])
        assert CLASS_LABELS[class_of(idx)] == labels[_fold_map[t]], mt
    assert context_index("N", "C", "A", "T") == UNKNOWN
    print("ok: 96 classes in COSMIC order; the two strands fold together and "
          "agree with folded_index()")

    # ---- parent-genome flanks on a hand-built tree -------------------
    ref = "ACGTACGTAC"
    ### positions:  1A 2C 3G 4T 5A 6C 7G 8T 9A 10C
    root = N("root", [], [
        N("a", ["T4C"], [                 # a changes 4 to C ...
            N("a1", ["A5G"]),             # ... so a1's left flank is C: C[A>G]C
            N("a2", ["C4T"], [            # reversion at 4 ...
                N("a2x", ["A5T"])])]),    # ... flank back to T: T[A>T]C
        N("b", ["A5G"]),                  # sibling of a: flank is ref T
        N("c", ["C6T", "G7A"]),           # adjacent on one branch: doublets
        N("d", ["C2A"]),                  # A[C>A]G
        N("e", ["A1G"]),                  # off the end: unknown
    ])
    tree = extract(root)
    ctx, stats = call_contexts(tree.parent, tree.children, tree.positions,
                               tree.types, ref)
    idx = {nid: i for i, nid in enumerate(tree.node_ids)}
    want = {
        "a1": context_index("C", "A", "C", "G"),
        "a2x": context_index("T", "A", "C", "T"),
        "b": context_index("T", "A", "C", "G"),
        "d": context_index("A", "C", "G", "A"),
        "e": UNKNOWN,
    }
    for k, v in want.items():
        assert ctx[idx[k]][0] == v, (k, ctx[idx[k]][0], v)
    assert list(ctx[idx["c"]]) == [DOUBLET, DOUBLET]
    assert stats["ref_mismatch"] == 0, stats
    ### a1's left flank is changed from the reference; a2x's has reverted.
    assert stats["flank_changed"] == 1, stats
    print("ok: flanks come from the parent's genome, are undone between "
          "siblings, follow reversions, and doublets are marked")

    ### A mask on a flank makes the context unknown.
    masked = read_mask([], len(ref))
    masked[3] = 1
    ctx2, _ = call_contexts(tree.parent, tree.children, tree.positions,
                            tree.types, ref, masked)
    assert ctx2[idx["d"]][0] == UNKNOWN
    print("ok: a masked flank makes the context unknown")

    ### A wrong reference shows up as ref_mismatch.
    _c, bad = call_contexts(tree.parent, tree.children, tree.positions,
                            tree.types, "GGGGGGGGGG")
    assert bad["ref_mismatch"] > 0
    print("ok: a reference that is not the root's genome is detected")

    # ---- random trees against brute force ---------------------------
    rng = random.Random(3)
    for rep in range(30):
        L = 60
        ref = "".join(rng.choice("ACGTN" if rep % 5 == 0 else "ACGT")
                      for _ in range(L))
        masked = read_mask([], L)
        if rep % 3 == 0:
            for p in rng.sample(range(1, L + 1), 6):
                masked[p] = 1
        genome = {}
        counter = [0]

        def grow(depth, genome):
            muts = []
            used = set()
            for _ in range(rng.randint(0, 4)):
                p = rng.randint(1, L)
                if p in used:
                    continue
                cur = genome.get(p, ref[p - 1])
                if cur not in COMPLEMENT:
                    continue
                alt = rng.choice([b for b in BASES if b != cur])
                muts.append(f"{cur}{p}{alt}")
                used.add(p)
            g2 = dict(genome)
            for m in muts:
                g2[int(m[1:-1])] = m[-1]
            counter[0] += 1
            kids = [grow(depth + 1, g2) for _ in range(
                rng.randint(0, 3) if depth < 4 else 0)]
            return N(f"n{counter[0]}", muts, kids)

        tree = extract(grow(0, genome))
        got, st = call_contexts(tree.parent, tree.children, tree.positions,
                                tree.types, ref, masked)
        want = _brute(tree, ref, masked)
        assert [list(a) for a in got] == want, rep
        assert st["ref_mismatch"] == 0, (rep, st)
    print("ok: 30 random trees with reversions, masks and Ns match a "
          "brute-force walk up the tree")

    # ---- opportunity -------------------------------------------------
    opp = opportunity("ACGTT")
    ### centres 2..4: ACG -> xCy ACG; CGT -> G centre -> ACG; GTT -> GTT.
    want = [0] * 32
    want[BASE_INDEX["A"] * 4 + BASE_INDEX["G"]] = 2
    want[16 + BASE_INDEX["G"] * 4 + BASE_INDEX["T"]] = 1
    assert opp == want, opp
    m = read_mask([], 5)
    m[1] = 1
    assert sum(opportunity("ACGTT", m)) == 2
    print("ok: opportunity folds purine centres and skips masked or "
          "unknown trinucleotides")
    print("all self-tests passed")
    return 0


if __name__ == "__main__":
    if "--selftest" in sys.argv[1:]:
        sys.exit(selftest())
    print(__doc__)
    print("This module is imported by profile_tree.py and "
          "gene_hit_contrasts.py; run it with --selftest to check it.")
    sys.exit(0)
