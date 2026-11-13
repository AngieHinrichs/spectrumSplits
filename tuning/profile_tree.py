#!/usr/bin/env python3
"""Profile a mutation-annotated tree to choose spectrumSplits.py parameters.

This is step 1 of parameter tuning: measure the scale of the tree before
guessing thresholds.  spectrumSplits.py's chi-squared statistic grows roughly
in proportion to (size of spectral difference) x (effective mutation count),
so the thresholds that work on one tree do not transfer to a tree with a
different mutation budget, different clade sizes, or a different spectrum
shape.

Reports:
  * tree size (nodes, leaves, internal nodes)
  * total mutations, and the per-branch mutation count distribution
    (informs --max_branch_length)
  * the subtree mutation count distribution, and how many nodes remain
    eligible as split candidates at a range of --min_mutations values
  * the global 12-type spectrum, the strand-folded 6-type spectrum, and the
    minimum --min_mutations needed to keep expected chi-squared cell counts
    at or above 5 (informs --min_mutations)
  * the mutations:tips ratio distribution (informs qc/prune_mutation_sample_ratio.py)

Nothing is modified; the tree is read only.

Usage:
    profile_tree.py --input_tree tree.pb.gz
    profile_tree.py --input_tree tree.pb.gz --node_table nodes.tsv
"""

import argparse
import math
import sys
from collections import Counter, defaultdict

# Must match the ordering used in spectrumSplits.py find_splits().
MUT_TYPES = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
MUT_INDEX = {m: i for i, m in enumerate(MUT_TYPES)}

# Strand-complementary pairs, for the folded 6-type spectrum.
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
FOLDED_TYPES = ["CA|GT", "CG|GC", "CT|GA", "TA|AT", "TC|AG", "TG|AC"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Profile a mutation-annotated tree to choose spectrumSplits.py parameters."
    )
    parser.add_argument("--input_tree", type=str, required=True,
                        help="Input tree file (protobuf format)")
    parser.add_argument("--node_table", type=str, default=None,
                        help="Optional TSV of per-node statistics (one row per node)")
    parser.add_argument("--min_mutations_grid", type=str,
                        default="100,250,500,1000,2500,5000,10000",
                        help="Comma-separated --min_mutations values to evaluate")
    parser.add_argument("--max_branch_length_grid", type=str,
                        default="10,25,50,100,250,500,1000",
                        help="Comma-separated --max_branch_length values to evaluate")
    parser.add_argument("--top_branches", type=int, default=20,
                        help="Number of longest branches to list")
    return parser.parse_args()


def parse_mutation(mutation):
    """Return (position, 'XY') for a SNV string like 'C241T', else None.

    Mutations involving ambiguous or non-ACGT characters are rejected so they
    do not pollute the spectrum; spectrumSplits.py silently buckets them into
    extra mutation types, which the 12-type chi-squared test then ignores.
    """
    if len(mutation) < 3:
        return None
    ref, alt = mutation[0].upper(), mutation[-1].upper()
    if ref not in COMPLEMENT or alt not in COMPLEMENT or ref == alt:
        return None
    try:
        pos = int(mutation[1:-1])
    except ValueError:
        return None
    return pos, ref + alt


def index_tree(root):
    """Flatten the tree into parallel arrays.

    Returns (nodes, parent, children) where every parent index is strictly
    less than its children's indices, so iterating indices in reverse visits
    children before parents (a valid post-order for accumulation).
    """
    nodes = []
    parent = []
    children = []
    stack = [(root, -1)]
    while stack:
        node, p = stack.pop()
        i = len(nodes)
        nodes.append(node)
        parent.append(p)
        children.append([])
        if p >= 0:
            children[p].append(i)
        for child in node.children:
            stack.append((child, i))
    return nodes, parent, children


def quantiles(sorted_values, probs):
    """Linear-interpolated quantiles of an already-sorted list."""
    if not sorted_values:
        return {p: float("nan") for p in probs}
    out = {}
    n = len(sorted_values)
    for p in probs:
        if n == 1:
            out[p] = float(sorted_values[0])
            continue
        h = (n - 1) * p
        lo = int(math.floor(h))
        hi = min(lo + 1, n - 1)
        out[p] = float(sorted_values[lo]) + (h - lo) * float(sorted_values[hi] - sorted_values[lo])
    return out


def describe(name, values, out):
    """Print a quantile summary of an unsorted list of numbers."""
    if not values:
        print(f"{name}: (none)", file=out)
        return
    sv = sorted(values)
    probs = [0.0, 0.5, 0.9, 0.99, 0.999, 1.0]
    q = quantiles(sv, probs)
    total = sum(sv)
    print(f"{name}", file=out)
    print(f"  n={len(sv)}  sum={total}  mean={total / len(sv):.3f}", file=out)
    print(f"  min={q[0.0]:.0f}  median={q[0.5]:.0f}  p90={q[0.9]:.0f}  "
          f"p99={q[0.99]:.0f}  p99.9={q[0.999]:.0f}  max={q[1.0]:.0f}", file=out)


def fold(mut_type):
    """Map a 12-type label onto its strand-folded pair label."""
    ref, alt = mut_type[0], mut_type[1]
    a = ref + alt
    b = COMPLEMENT[ref] + COMPLEMENT[alt]
    # Canonical label is the pair with a C or T reference first.
    if ref in ("C", "T"):
        return a + "|" + b
    return b + "|" + a


def main():
    args = parse_args()
    out = sys.stdout

    try:
        import bte
    except ImportError:
        print("error: could not import bte; run this in the environment you use "
              "for spectrumSplits.py", file=sys.stderr)
        return 1

    print(f"Reading {args.input_tree} ...", file=sys.stderr)
    tree = bte.MATree(args.input_tree)
    nodes, parent, children = index_tree(tree.root)
    n_nodes = len(nodes)
    print(f"Indexed {n_nodes} nodes", file=sys.stderr)

    # ---- per-branch pass -------------------------------------------------
    own_count = [0] * n_nodes          # usable SNVs on this branch
    own_spectrum = [None] * n_nodes    # list of 12 counts, lazily created
    raw_count = [0] * n_nodes          # every mutation string on this branch
    global_spectrum = [0] * 12
    site_counts = Counter()
    skipped = Counter()
    n_leaves = 0

    for i, node in enumerate(nodes):
        muts = node.mutations
        raw_count[i] = len(muts)
        if not children[i]:
            n_leaves += 1
        spec = None
        for mutation in muts:
            parsed = parse_mutation(mutation)
            if parsed is None:
                skipped[mutation[:1] + "..." + mutation[-1:] if len(mutation) > 2 else mutation] += 1
                continue
            pos, mut_type = parsed
            idx = MUT_INDEX[mut_type]
            if spec is None:
                spec = [0] * 12
            spec[idx] += 1
            own_count[i] += 1
            global_spectrum[idx] += 1
            site_counts[pos] += 1
        own_spectrum[i] = spec

    total_mutations = sum(global_spectrum)
    total_skipped = sum(skipped.values())

    # ---- subtree accumulation (post-order over reversed indices) ---------
    subtree_mutations = list(own_count)
    subtree_tips = [1 if not children[i] else 0 for i in range(n_nodes)]
    for i in range(n_nodes - 1, 0, -1):
        p = parent[i]
        subtree_mutations[p] += subtree_mutations[i]
        subtree_tips[p] += subtree_tips[i]

    # ---- report ----------------------------------------------------------
    print("=" * 72, file=out)
    print("TREE SIZE", file=out)
    print("=" * 72, file=out)
    print(f"  nodes:    {n_nodes}", file=out)
    print(f"  leaves:   {n_leaves}", file=out)
    print(f"  internal: {n_nodes - n_leaves}", file=out)
    print(f"  usable SNVs (ACGT->ACGT): {total_mutations}", file=out)
    print(f"  skipped mutation strings: {total_skipped}", file=out)
    if total_skipped:
        print("  WARNING: every mutation on an UShER tree should be an unambiguous", file=out)
        print("  [ACGT]->[ACGT] substitution, because UShER imputes a base for each", file=out)
        print("  ambiguous input base. A non-zero count here means the tree is not", file=out)
        print("  what spectrumSplits.py expects; investigate before going further.", file=out)
    if skipped:
        print("  most common skipped patterns:", file=out)
        for pat, c in skipped.most_common(5):
            print(f"    {pat}\t{c}", file=out)
    print(f"  distinct mutated sites: {len(site_counts)}", file=out)
    print(file=out)

    print("=" * 72, file=out)
    print("GLOBAL SPECTRUM", file=out)
    print("=" * 72, file=out)
    print("  12-type (the categories spectrumSplits.py tests):", file=out)
    min_p_12 = 1.0
    for i, mut_type in enumerate(MUT_TYPES):
        p = global_spectrum[i] / total_mutations if total_mutations else 0.0
        min_p_12 = min(min_p_12, p) if p > 0 else min_p_12
        print(f"    {mut_type}\t{global_spectrum[i]:>12}\t{p:7.4f}", file=out)

    folded = defaultdict(int)
    for i, mut_type in enumerate(MUT_TYPES):
        folded[fold(mut_type)] += global_spectrum[i]
    print("  6-type, strand-folded (the categories to consider for TB):", file=out)
    min_p_6 = 1.0
    for label in FOLDED_TYPES:
        c = folded.get(label, 0)
        p = c / total_mutations if total_mutations else 0.0
        if p > 0:
            min_p_6 = min(min_p_6, p)
        print(f"    {label}\t{c:>12}\t{p:7.4f}", file=out)
    print(file=out)
    print("  Expected-cell-count floor on --min_mutations:", file=out)
    print("    The chi-squared approximation degrades when an expected cell count", file=out)
    print("    falls below ~5.  The binding cell is the rarest mutation type on", file=out)
    print("    the smaller side of the split, so the smaller side needs at least", file=out)
    print("    5 / p_rarest mutations.", file=out)
    if min_p_12 < 1.0:
        print(f"    12-type: rarest class p={min_p_12:.5f} -> "
              f"--min_mutations >= {math.ceil(5 / min_p_12)}", file=out)
    if min_p_6 < 1.0:
        print(f"    6-type:  rarest class p={min_p_6:.5f} -> "
              f"--min_mutations >= {math.ceil(5 / min_p_6)}", file=out)
    print(file=out)

    print("=" * 72, file=out)
    print("PER-BRANCH MUTATION COUNTS  (informs --max_branch_length)", file=out)
    print("=" * 72, file=out)
    describe("  branch mutation counts, all nodes", own_count, out)
    internal = [own_count[i] for i in range(n_nodes) if children[i]]
    leafward = [own_count[i] for i in range(n_nodes) if not children[i]]
    describe("  branch mutation counts, internal nodes", internal, out)
    describe("  branch mutation counts, terminal branches", leafward, out)
    print(file=out)
    print("  Effect of capping --max_branch_length "
          "(branches above the cap contribute nothing):", file=out)
    print("    cap\tbranches_excluded\tmutations_excluded\tpct_of_total", file=out)
    for cap_s in args.max_branch_length_grid.split(","):
        cap_s = cap_s.strip()
        if not cap_s:
            continue
        cap = int(cap_s)
        n_ex = 0
        m_ex = 0
        for i in range(n_nodes):
            if own_count[i] > cap:
                n_ex += 1
                m_ex += own_count[i]
        pct = 100.0 * m_ex / total_mutations if total_mutations else 0.0
        print(f"    {cap}\t{n_ex}\t{m_ex}\t{pct:.2f}%", file=out)
    print(file=out)
    order = sorted(range(n_nodes), key=lambda i: own_count[i], reverse=True)
    print(f"  {args.top_branches} longest branches:", file=out)
    print("    node_id\tbranch_mutations\tsubtree_tips", file=out)
    for i in order[:args.top_branches]:
        print(f"    {nodes[i].id}\t{own_count[i]}\t{subtree_tips[i]}", file=out)
    print(file=out)

    print("=" * 72, file=out)
    print("SPLIT CANDIDATE AVAILABILITY  (informs --min_mutations)", file=out)
    print("=" * 72, file=out)
    describe("  subtree mutation counts, all nodes", subtree_mutations, out)
    print(file=out)
    print("  A node can be tested as a split of the whole tree only if both", file=out)
    print("  sides clear --min_mutations (spectrumSplits.py lines 160-163).", file=out)
    print("  Counts below are for the first pass, before any splitting; later", file=out)
    print("  passes work on smaller subtrees and so have fewer candidates.", file=out)
    print("    min_mutations\teligible_nodes\tpct_of_nodes\tsmallest_eligible_subtree_tips", file=out)
    for thr_s in args.min_mutations_grid.split(","):
        thr_s = thr_s.strip()
        if not thr_s:
            continue
        thr = int(thr_s)
        eligible = 0
        min_tips = None
        for i in range(n_nodes):
            below = subtree_mutations[i]
            above = total_mutations - below
            if below >= thr and above >= thr:
                eligible += 1
                if min_tips is None or subtree_tips[i] < min_tips:
                    min_tips = subtree_tips[i]
        pct = 100.0 * eligible / n_nodes if n_nodes else 0.0
        print(f"    {thr}\t{eligible}\t{pct:.3f}%\t"
              f"{min_tips if min_tips is not None else '-'}", file=out)
    print(file=out)

    print("=" * 72, file=out)
    print("MUTATIONS:TIPS RATIO  (informs qc/prune_mutation_sample_ratio.py)", file=out)
    print("=" * 72, file=out)
    ratios = [subtree_mutations[i] / subtree_tips[i]
              for i in range(n_nodes) if subtree_tips[i] > 0]
    describe("  subtree mutations / subtree tips", ratios, out)
    if ratios:
        mean_r = sum(ratios) / len(ratios)
        var_r = sum((r - mean_r) ** 2 for r in ratios) / len(ratios)
        sd_r = math.sqrt(var_r)
        auto_thr = mean_r + 2 * sd_r
        n_over = sum(1 for r in ratios if r >= auto_thr)
        print(f"  prune_mutation_sample_ratio.py default threshold "
              f"(mean + 2*SD) = {auto_thr:.3f}", file=out)
        print(f"  nodes at or above that threshold: {n_over}", file=out)
        print("  Inspect what this prunes before accepting it: on a sparsely", file=out)
        print("  sampled, diverged tree a high ratio is often real divergence", file=out)
        print("  rather than an artifact.", file=out)
    print(file=out)

    print("=" * 72, file=out)
    print("SITE RECURRENCE  (independence violations that inflate chi-squared)", file=out)
    print("=" * 72, file=out)
    site_vals = list(site_counts.values())
    describe("  mutations per site", site_vals, out)
    if total_mutations:
        top_sites = site_counts.most_common(20)
        share = 100.0 * sum(c for _, c in top_sites) / total_mutations
        print(f"  top 20 sites hold {share:.2f}% of all mutations", file=out)
        print("    position\tmutations\tpct_of_total", file=out)
        for pos, c in top_sites:
            print(f"    {pos}\t{c}\t{100.0 * c / total_mutations:.3f}%", file=out)
    print("  The chi-squared test assumes mutations are independent draws.", file=out)
    print("  Sites mutating repeatedly within one clade break that assumption", file=out)
    print("  and are the main reason --min_chi must exceed any nominal", file=out)
    print("  significance cutoff.  qc/mask_site_splits.py targets these.", file=out)
    print(file=out)

    if args.node_table:
        with open(args.node_table, "w") as fh:
            header = (["node_id", "is_leaf", "branch_mutations", "raw_branch_mutations",
                       "subtree_mutations", "subtree_tips", "mutations_per_tip"]
                      + MUT_TYPES)
            fh.write("\t".join(header) + "\n")
            for i in range(n_nodes):
                spec = own_spectrum[i] or [0] * 12
                ratio = (subtree_mutations[i] / subtree_tips[i]) if subtree_tips[i] else 0.0
                row = ([nodes[i].id,
                        "1" if not children[i] else "0",
                        own_count[i], raw_count[i],
                        subtree_mutations[i], subtree_tips[i],
                        f"{ratio:.4f}"]
                       + spec)
                fh.write("\t".join(str(x) for x in row) + "\n")
        print(f"Per-node table written to {args.node_table}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
