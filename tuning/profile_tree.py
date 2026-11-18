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
  * with --reference or --terminus: strand asymmetry by replication arm, by
    gene strand (with --annotation), and its spread across annotated clades
    (informs --fold; see the STRAND ASYMMETRY section's notes)

Nothing is modified; the tree is read only.

Usage:
    profile_tree.py --input_tree tree.pb.gz
    profile_tree.py --input_tree tree.pb.gz --node_table nodes.tsv
    profile_tree.py --input_tree tree.pb.gz --reference ref.fa --annotation genes.gff3
"""

import argparse
import gzip
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
    parser.add_argument("--reference", type=str, default=None,
                        help="Reference FASTA in the tree's coordinates, optionally gzipped.  Enables the "
                        "STRAND ASYMMETRY section; the origin and terminus of "
                        "replication are estimated from cumulative GC skew unless "
                        "--origin/--terminus are given")
    parser.add_argument("--origin", type=int, default=None,
                        help="Replication origin position (1-based).  Default: the "
                        "minimum of cumulative G-C skew with --reference, else 1")
    parser.add_argument("--terminus", type=int, default=None,
                        help="Replication terminus position (1-based).  Default: the "
                        "maximum of cumulative G-C skew; required for the STRAND "
                        "ASYMMETRY section when --reference is not given")
    parser.add_argument("--annotation", type=str, default=None,
                        help="GFF3 annotation, optionally gzipped; adds a gene-strand breakdown to the "
                        "STRAND ASYMMETRY section (gene and pseudogene features)")
    parser.add_argument("--strand_min_clade", type=int, default=25000,
                        help="Smallest clade, in mutations, listed in the per-clade "
                        "strand asymmetry table")
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


def open_text(path):
    """Open a plain or gzipped text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path):
    """Concatenated, upper-cased sequence of every record in a FASTA file."""
    parts = []
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def gc_skew_extremes(seq):
    """(argmin, argmax) of cumulative G-C count, as 1-based positions.

    In most bacteria the leading strand is G-rich, so walking the forward
    strand the cumulative G-C count falls to its minimum at the origin and
    rises to its maximum at the terminus.
    """
    c = 0
    lo, lo_pos, hi, hi_pos = 0, 1, 0, 1
    for i, b in enumerate(seq, 1):
        if b == "G":
            c += 1
        elif b == "C":
            c -= 1
        if c < lo:
            lo, lo_pos = c, i
        if c > hi:
            hi, hi_pos = c, i
    return lo_pos, hi_pos


def forward_is_leading(pos, origin, terminus):
    """True if the forward strand is the leading strand at pos.

    That is the arm replicated travelling forward from the origin to the
    terminus, wrapping around the end of a circular genome when the origin
    lies past the terminus.
    """
    if origin <= terminus:
        return origin <= pos < terminus
    return pos >= origin or pos < terminus


def read_gene_strands(path):
    """bytearray indexed by 1-based position: 0 none, 1 +, 2 -, 3 both."""
    spans = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene") or f[6] not in "+-":
                continue
            spans.append((int(f[3]), int(f[4]), 1 if f[6] == "+" else 2))
    strands = bytearray((max(e for _, e, _ in spans) if spans else 0) + 2)
    for start, end, bit in spans:
        for p in range(start, end + 1):
            strands[p] |= bit
    return strands


def complement_type(mut_type):
    return COMPLEMENT[mut_type[0]] + COMPLEMENT[mut_type[1]]


def pair_rows(spectrum):
    """(label, a, b, log2(a/b), SE) for each folded pair of a 12-count list."""
    rows = []
    for label in FOLDED_TYPES:
        x, y = label.split("|")
        a, b = spectrum[MUT_INDEX[x]], spectrum[MUT_INDEX[y]]
        if a and b:
            lr = math.log2(a / b)
            se = math.sqrt(1 / a + 1 / b) / math.log(2)
        else:
            lr = se = float("nan")
        rows.append((label, a, b, lr, se))
    return rows


def print_pair_table(title, spectrum, out):
    print(f"  {title}  (n={sum(spectrum)})", file=out)
    print("    pair\tfirst\tsecond\tlog2(first/second)\tSE", file=out)
    for label, a, b, lr, se in pair_rows(spectrum):
        print(f"    {label}\t{a}\t{b}\t{lr:+.3f}\t{se:.3f}", file=out)


def homogeneity_chi2(pairs):
    """Pearson chi-squared for a 2 x K table given as [(a_k, b_k), ...]."""
    ta = sum(a for a, _ in pairs)
    tb = sum(b for _, b in pairs)
    n = ta + tb
    if not ta or not tb:
        return 0.0
    chi = 0.0
    for a, b in pairs:
        m = a + b
        ea, eb = m * ta / n, m * tb / n
        chi += (a - ea) ** 2 / ea + (b - eb) ** 2 / eb
    return chi


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

    # ---- strand-asymmetry setup -------------------------------------------
    do_strand = args.reference is not None or args.terminus is not None
    origin, terminus, skew_note = args.origin, args.terminus, ""
    gene_strands = None
    if do_strand:
        if args.reference:
            seq = read_fasta(args.reference)
            skew_min, skew_max = gc_skew_extremes(seq)
            if origin is None:
                origin = skew_min
                skew_note += " origin from GC skew;"
            if terminus is None:
                terminus = skew_max
                skew_note += " terminus from GC skew;"
            del seq
        if origin is None:
            origin = 1
        if args.annotation:
            gene_strands = read_gene_strands(args.annotation)
    arm_spectrum = {"leading": [0] * 12, "lagging": [0] * 12}   # forward-strand counts
    gene_spectrum = {k: [0] * 12 for k in ("+ gene", "- gene", "both strands", "intergenic")}
    clade_leading = defaultdict(lambda: [0] * 12)   # leading-strand orientation
    clade_of = [""] * n_nodes
    gene_label = {0: "intergenic", 1: "+ gene", 2: "- gene", 3: "both strands"}

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
        if do_strand:
            clade_of[i] = next((a for a in node.annotations if a), "") \
                or (clade_of[parent[i]] if parent[i] >= 0 else "")
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
            if do_strand:
                leading = forward_is_leading(pos, origin, terminus)
                arm_spectrum["leading" if leading else "lagging"][idx] += 1
                oriented = idx if leading else MUT_INDEX[complement_type(mut_type)]
                clade_leading[clade_of[i]][oriented] += 1
                if gene_strands is not None:
                    st = gene_strands[pos] if pos < len(gene_strands) else 0
                    gene_spectrum[gene_label[st]][idx] += 1
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

    if do_strand:
        print("=" * 72, file=out)
        print("STRAND ASYMMETRY  (informs --fold)", file=out)
        print("=" * 72, file=out)
        print(f"  origin {origin}, terminus {terminus}"
              + (f"  ({skew_note.strip(' ;')})" if skew_note else ""), file=out)
        print("  Each pair's log2 ratio compares a mutation type with its", file=out)
        print("  complement, counted on the forward strand.  0 means symmetric.", file=out)
        print(file=out)
        print_pair_table("forward strand is LEADING (origin -> terminus)",
                         arm_spectrum["leading"], out)
        print_pair_table("forward strand is LAGGING",
                         arm_spectrum["lagging"], out)
        print(file=out)
        if gene_strands is not None:
            for key in ("+ gene", "- gene", "intergenic"):
                print_pair_table(f"forward strand, {key}", gene_spectrum[key], out)
            print(f"  ({sum(gene_spectrum['both strands'])} mutations in genes on both "
                  "strands left out)", file=out)
            print(file=out)

        clades = [(c, s) for c, s in clade_leading.items()
                  if sum(s) >= args.strand_min_clade]
        clades.sort(key=lambda cs: -sum(cs[1]))
        print("  Per clade, leading-strand orientation (lagging-arm mutations", file=out)
        print("  complemented), log2(first/second) for each pair.  Mutations are", file=out)
        print("  assigned to the nearest annotated ancestor, so a parent clade's", file=out)
        print("  row excludes its annotated subclades; (none) is mutations above", file=out)
        print(f"  every annotation.  Clades with >= {args.strand_min_clade} mutations:", file=out)
        print("    clade\tmutations\t" + "\t".join(FOLDED_TYPES), file=out)
        for c, s in clades:
            vals = [f"{lr:+.3f}" for _, _, _, lr, _ in pair_rows(s)]
            print(f"    {c or '(none)'}\t{sum(s)}\t" + "\t".join(vals), file=out)
        if len(clades) >= 2:
            print("  Heterogeneity across those clades (2 x K chi-squared on each", file=out)
            print("  pair's leading-oriented counts; chi2/df near 1 means the clades", file=out)
            print("  agree within sampling noise; spread = max - min log2 ratio, which", file=out)
            print("  includes that noise for the smaller clades):", file=out)
            print("    pair\tchi2\tdf\tchi2/df\tspread\tlowest\thighest", file=out)
            for j, label in enumerate(FOLDED_TYPES):
                x, y = label.split("|")
                counts = [(s[MUT_INDEX[x]], s[MUT_INDEX[y]]) for _, s in clades]
                lrs = [(pair_rows(s)[j][3], c or "(none)") for c, s in clades]
                lrs = [t for t in lrs if not math.isnan(t[0])]
                chi = homogeneity_chi2(counts)
                lo, hi = min(lrs), max(lrs)
                df = len(counts) - 1
                print(f"    {label}\t{chi:.1f}\t{df}\t{chi / df:.1f}\t{hi[0] - lo[0]:.3f}\t"
                      f"{lo[1]} {lo[0]:+.3f}\t{hi[1]} {hi[0]:+.3f}", file=out)
        print(file=out)
        print("  Reading this:", file=out)
        print("  * A pair whose ratio flips sign between the LEADING and LAGGING arms", file=out)
        print("    (or between + and - genes) marks a strand-asymmetric process", file=out)
        print("    that cancels in the GLOBAL SPECTRUM.  Near-equal global pairs do", file=out)
        print("    not show that no such process exists.", file=out)
        print("  * spectrumSplits.py counts forward-strand types, so the same", file=out)
        print("    cancellation happens inside every clade whose mutations span both", file=out)
        print("    arms.  A clade-to-clade difference in this asymmetry is invisible", file=out)
        print("    with or without --fold; unfolding adds degrees of freedom that", file=out)
        print("    carry mostly noise and arm-imbalance artifacts.", file=out)
        print("  * Heterogeneity chi-squared is inflated by recurrent sites like any", file=out)
        print("    statistic on this tree; read the spread as the effect size.", file=out)
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
