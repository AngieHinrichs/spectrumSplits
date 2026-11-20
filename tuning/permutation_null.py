#!/usr/bin/env python3
"""Calibrate spectrumSplits.py's --min_chi against a permutation null.

This is step 7 of parameter tuning.  spectrumSplits.py scores every node in a
subtree and keeps the single best chi-squared value, so the largest chi-squared
it finds is a maximum over hundreds of thousands of correlated tests.  A
nominal significance cutoff is therefore meaningless: --min_chi has to clear
whatever the search turns up when there is no clade-specific spectrum at all.

This script measures that.  It destroys any association between mutation type
and position in the tree, then runs the same search and the same statistic
spectrumSplits.py uses, and reports the chi-squared values it still finds.
Set --min_chi comfortably above them.

The printed summary's recommendation rests on a --quantile (default 0.95) of
the replicates' maxima, not on the largest one.  The largest of n draws grows
with n and is set by a single replicate, so a threshold built on it moves
with the replicate count and jumps between runs.  A quantile stays put as
replicates are added and only gets more precise.  The largest of 20
replicates sits near the 95th percentile, so the default keeps thresholds
comparable with 20-replicate runs; it needs about 100 replicates to resolve.

Null models (--mode):

  types  Shuffle mutation-type labels across all mutations in the tree, keeping
         the topology and the number of mutations per branch fixed.  This is
         the strict independence null: it removes clade-specific spectra and
         also removes the type clumping caused by recurrent sites, so it tends
         to give the more permissive (lower) threshold.

  sites  Give every site a single mutation type (its most common one), then
         shuffle that site-to-type assignment across sites.  Which branches
         carry mutations at which site is untouched, so a site that mutates
         500 times inside one clade still dumps 500 counts of one type into
         that clade - the type is just randomized.  This preserves the
         independence violation that recurrent sites cause, and so gives the
         more conservative (higher) threshold.  This is usually the one to
         calibrate against, especially for a tree with convergent selection.

The statistic matches spectrumSplits.py: a 2 x 12 contingency table of the
candidate subtree's spectrum against the rest of the current split, scored
with Pearson's chi-squared (no continuity correction, as scipy applies none
for tables larger than 2 x 2).  Columns whose total is zero are skipped, where
scipy would raise; this only arises when a mutation type is entirely absent.
Pass --fold to score the strand-folded 6-type spectrum instead.

One invocation, the whole of step 7.  By default (--mode both) it runs both
nulls with the same --seed, so their replicates are paired: `types` over the
whole tree, and `sites` following the greedy cascade --max_splits (default
25) deep with --min_chi 1, so it also measures the null for smaller
subtrees.  The unpermuted tree is scored through the same code path as a
reference.  The report then has four sections:

  1. PERMUTATION NULL SUMMARY.  Each null's replicate maxima side by side,
     how far `sites` sits above `types` (what residual recurrent-site
     clumping costs), and the recommended --min_chi: --margin (default 2)
     times the --quantile of the `sites` maxima.  Always calibrate against
     `sites`.

  2. NULL CEILING BY SUBTREE SIZE.  The cascade's null splits in
     factor-of-2 bins of the size of the table each was scored on, with the
     median, --null_quantile (default 0.99) and maximum chi-squared per bin,
     the lowest --min_chi that clears the whole-tree null, and that
     threshold's margin over each bin.  Chi-squared scales with the counts in
     the table, so a threshold calibrated on the whole tree is too strict for
     deeper splits; this shows by how much.

  3. OBSERVED SPLITS AGAINST A SIZE-MATCHED NULL.  Each observed cascade
     split against the --null_quantile of null splits within --null_window
     fold of its size (split_confidence.py's ceiling).  Splits the
     recommended --min_chi (or --flag_min_chi) rejects though they clear
     their null by --margin, and splits it accepts that don't clear it, are
     flagged.  The observed cascade accepts every split at --min_chi 1, so
     deep rows approximate what a production run would test; for the
     production run's own splits use split_confidence.py (step 10).

  4. HOTSPOT CHECK.  Under the sites null each site keeps its branches, so
     a heavily recurrent site inside one clade dumps clumped counts there in
     every replicate, and wins big when it draws a rare type.  A clade that
     is the whole-tree maximum in more than --hotspot_fraction of the
     replicates, or holds the --top_maxima largest maxima, is flagged, with
     its most recurrent sites and an annotate_sites.py command for them.  A
     node and a child holding nearly all of it (--same_clade) count as one
     clade.  If a hotspot is real, mask it and rerun: it is setting the
     threshold.

Outputs, named from --prefix (default null): PREFIX_types.tsv and
PREFIX_sites_cascade.tsv (PREFIX_sites.tsv with --max_splits 1), every split
scored, per replicate; step 10's split_confidence.py --null reads the
cascade.  The report goes to stdout.

Usage:
    permutation_null.py --input_tree tree.pb.gz --fold --min_mutations 500 \
        --max_branch_length 250 --replicates 100 --seed 1 --nthreads 8 \
        --annotation genes.gff3.gz --reference ref.fa > null_report.txt

    permutation_null.py --selftest
"""

import argparse
import collections
import math
import os
import random
import statistics
import sys
from array import array
from collections import Counter, defaultdict

# Must match the ordering used in spectrumSplits.py find_splits().
MUT_TYPES = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
MUT_INDEX = {m: i for i, m in enumerate(MUT_TYPES)}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def folded_label(mut_type):
    """Label for the strand-complementary pair a mutation class belongs to."""
    ref, alt = mut_type[0], mut_type[1]
    reverse = COMPLEMENT[ref] + COMPLEMENT[alt]
    if ref in ("C", "T"):
        return mut_type + "|" + reverse
    return reverse + "|" + mut_type


def folded_index():
    """Map each of the 12 type indices onto one of 6 strand-folded indices.

    Labels are sorted so the category order matches profile_tree.py and
    spectrumSplits.py; the statistic is order-independent but the reported
    order should not differ between the tools.
    """
    labels = sorted({folded_label(m) for m in MUT_TYPES})
    mapping = [labels.index(folded_label(m)) for m in MUT_TYPES]
    return mapping, labels


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Calibrate spectrumSplits.py --min_chi against permutation "
                    "nulls, by subtree size, and check them for a hotspot.")
    parser.add_argument("--input_tree", type=str, default=None,
                        help="Input tree file (protobuf format)")
    parser.add_argument("--prefix", type=str, default="null",
                        help="Prefix of the per-split TSVs (default null: "
                             "null_types.tsv, null_sites_cascade.tsv)")
    parser.add_argument("--mode", type=str, default="both",
                        choices=["both", "sites", "types"],
                        help="Null model(s) to run (default both; see the "
                             "module docstring)")
    parser.add_argument("--replicates", type=int, default=100,
                        help="Number of permutation replicates (default 100)")
    parser.add_argument("--quantile", type=float, default=0.95,
                        help="Quantile of the replicate maxima the recommended "
                             "--min_chi rests on (default 0.95)")
    parser.add_argument("--margin", type=float, default=2.0,
                        help="The recommended --min_chi is this times the "
                             "--quantile; also the margin over its "
                             "size-matched null for a split below it to be "
                             "flagged as rejected (default 2)")
    parser.add_argument("--min_mutations", type=int, default=500,
                        help="Match the --min_mutations you intend to use for the real run")
    parser.add_argument("--max_branch_length", type=int, default=100000,
                        help="Match the --max_branch_length you intend to use")
    parser.add_argument("--min_chi", type=float, default=1.0,
                        help="Threshold the null searches stop at (default 1, "
                             "so the cascade runs its full depth; a higher "
                             "value cuts the null off at the depth being "
                             "measured)")
    parser.add_argument("--max_splits", type=int, default=25,
                        help="Splits the sites null follows per replicate "
                             "(default 25).  1 scores only the whole-tree "
                             "pass; the types null always scores only that")
    parser.add_argument("--fold", action="store_true",
                        help="Score the strand-folded 6-type spectrum instead of 12 types")
    parser.add_argument("--seed", type=int, default=1,
                        help="Base random seed; replicate r uses seed + r in "
                             "both nulls")
    parser.add_argument("--nthreads", type=int, default=1,
                        help="Run replicates concurrently (requires fork; falls back to serial)")
    g = parser.add_argument_group("size-matched null and hotspot check "
                                  "(from the sites cascade)")
    g.add_argument("--flag_min_chi", type=float, default=None,
                   help="The --min_chi to flag observed splits against "
                        "(default: the recommended one)")
    g.add_argument("--null_window", type=float, default=1.5,
                   help="Fold range of subtree size accepted as size-matched, "
                        "as in split_confidence.py (default 1.5)")
    g.add_argument("--null_quantile", type=float, default=0.99,
                   help="Quantile of the size-matched null taken as the "
                        "ceiling, as in split_confidence.py (default 0.99)")
    g.add_argument("--min_nulls", type=int, default=10,
                   help="Mark a ceiling drawn from fewer null splits than "
                        "this (default 10)")
    g.add_argument("--hotspot_fraction", type=float, default=0.5,
                   help="Flag a clade that is the whole-tree null maximum in "
                        "more than this fraction of replicates (default 0.5)")
    g.add_argument("--top_maxima", type=int, default=2,
                   help="Flag a clade that holds this many of the largest "
                        "whole-tree null maxima (default 2)")
    g.add_argument("--same_clade", type=float, default=0.95,
                   help="A nested node whose mutations_below is at least this "
                        "fraction of its ancestor's counts as the same clade "
                        "(default 0.95)")
    g.add_argument("--hotspot_sites", type=int, default=10,
                   help="Recurrent sites to list per flagged clade (default 10)")
    g.add_argument("--hotspot_min", type=int, default=20,
                   help="Sites with at least this many mutations in a flagged "
                        "clade go into the printed annotate_sites.py command "
                        "(default 20)")
    g.add_argument("--annotation", default=None,
                   help="Used only in the printed annotate_sites.py command")
    g.add_argument("--reference", default=None,
                   help="Used only in the printed annotate_sites.py command")
    parser.add_argument("--selftest", action="store_true",
                        help="Run the self-tests and exit")
    return parser.parse_args(argv)


# --------------------------------------------------------------------------
# chi-squared
# --------------------------------------------------------------------------

def chi2_2xk(below, above):
    """Pearson chi-squared for a 2 x k contingency table.

    Matches scipy.stats.chi2_contingency([below, above]) for k > 2, which
    applies no continuity correction.  Columns with a zero total are skipped.
    Returns 0.0 if either row total is zero.
    """
    rb = 0
    ra = 0
    for v in below:
        rb += v
    for v in above:
        ra += v
    n = rb + ra
    if rb <= 0 or ra <= 0 or n <= 0:
        return 0.0
    chi2 = 0.0
    for j in range(len(below)):
        col = below[j] + above[j]
        if col == 0:
            continue
        e_b = rb * col / n
        e_a = ra * col / n
        d_b = below[j] - e_b
        d_a = above[j] - e_a
        chi2 += d_b * d_b / e_b + d_a * d_a / e_a
    return chi2


def selftest():
    """Verify chi2_2xk against a naive implementation and against scipy."""
    rng = random.Random(12345)
    worst = 0.0
    for trial in range(200):
        k = rng.choice([6, 12])
        below = [rng.randint(0, 500) for _ in range(k)]
        above = [rng.randint(0, 5000) for _ in range(k)]
        mine = chi2_2xk(below, above)

        # naive reference
        rows = [below, above]
        row_tot = [sum(r) for r in rows]
        col_tot = [below[j] + above[j] for j in range(k)]
        n = sum(row_tot)
        ref = 0.0
        for i in range(2):
            for j in range(k):
                if col_tot[j] == 0:
                    continue
                e = row_tot[i] * col_tot[j] / n
                ref += (rows[i][j] - e) ** 2 / e
        if abs(mine - ref) > 1e-9 * max(1.0, ref):
            print(f"FAIL naive: trial {trial} mine={mine} ref={ref}", file=sys.stderr)
            return 1
        worst = max(worst, abs(mine - ref))
    print(f"ok: matches naive reference on 200 random tables (max abs diff {worst:.3e})")

    ### The recommendation rests on the quantile, not the single largest
    ### replicate, and warns when there are too few replicates to resolve it.
    import io
    import types as _types
    a = _types.SimpleNamespace(quantile=0.95, margin=2.0, min_mutations=10,
                               max_branch_length=250, max_splits=25,
                               min_chi=1.0)
    maxima = sorted([10.0 + i / 10 for i in range(99)] + [100.0])
    types_max = sorted(v * 0.9 for v in maxima)
    buf = io.StringIO()
    rec = summarize({"types": types_max, "sites": maxima}, 1000.0, a, 6, buf)
    text = buf.getvalue()
    q = quantile(maxima, 0.95)
    if (rec != math.ceil(2 * q)
            or f"Recommended (--margin 2): --min_chi {rec}" not in text
            or "--min_chi 200" in text or "Raise --replicates" in text
            or f"{1000 / q:.1f}x" not in text
            or "sites p95 is 1.11x types p95" not in text):
        print("FAIL: summary threshold\n" + text, file=sys.stderr)
        return 1
    buf = io.StringIO()
    summarize({"sites": maxima[:20]}, None, a, 6, buf)
    if "Raise --replicates" not in buf.getvalue():
        print("FAIL: no warning for 20 replicates\n" + buf.getvalue(),
              file=sys.stderr)
        return 1
    print("ok: the recommended --min_chi rests on the sites quantile, sites "
          "is set against types, and 20 replicates draw a warning")
    _selftest_size()

    try:
        from scipy.stats import chi2_contingency
    except ImportError:
        print("note: scipy not importable here; skipped the scipy comparison. "
              "Re-run --selftest in your spectrumSplits environment to confirm.")
        print("all self-tests passed")
        return 0
    worst = 0.0
    for trial in range(200):
        k = rng.choice([6, 12])
        below = [rng.randint(1, 500) for _ in range(k)]
        above = [rng.randint(1, 5000) for _ in range(k)]
        mine = chi2_2xk(below, above)
        ref = chi2_contingency([below, above])[0]
        if abs(mine - ref) > 1e-8 * max(1.0, ref):
            print(f"FAIL scipy: trial {trial} mine={mine} ref={ref}", file=sys.stderr)
            return 1
        worst = max(worst, abs(mine - ref))
    print(f"ok: matches scipy.stats.chi2_contingency on 200 random tables "
          f"(max abs diff {worst:.3e})")
    print("all self-tests passed")
    return 0


# --------------------------------------------------------------------------
# tree extraction
# --------------------------------------------------------------------------

class FlatTree(object):
    """Topology plus per-branch mutations, independent of bte.

    Parent indices are strictly less than child indices, so iterating indices
    in reverse visits children before parents.
    """

    def __init__(self):
        self.node_ids = []
        self.parent = array("i")
        self.children = []
        self.raw_len = array("i")      # len(node.mutations), as spectrumSplits sees it
        self.positions = []            # array('i') per node, parsed SNVs only
        self.types = []                # array('b') per node, aligned with positions


def parse_mutation(mutation):
    if len(mutation) < 3:
        return None
    ref, alt = mutation[0].upper(), mutation[-1].upper()
    if ref not in COMPLEMENT or alt not in COMPLEMENT or ref == alt:
        return None
    try:
        pos = int(mutation[1:-1])
    except ValueError:
        return None
    return pos, MUT_INDEX[ref + alt]


def extract(root):
    tree = FlatTree()
    stack = [(root, -1)]
    while stack:
        node, p = stack.pop()
        i = len(tree.node_ids)
        tree.node_ids.append(node.id)
        tree.parent.append(p)
        tree.children.append([])
        if p >= 0:
            tree.children[p].append(i)
        muts = node.mutations
        tree.raw_len.append(len(muts))
        pos_arr = array("i")
        type_arr = array("b")
        for mutation in muts:
            parsed = parse_mutation(mutation)
            if parsed is None:
                continue
            pos_arr.append(parsed[0])
            type_arr.append(parsed[1])
        tree.positions.append(pos_arr)
        tree.types.append(type_arr)
        for child in node.children:
            stack.append((child, i))
    return tree


# --------------------------------------------------------------------------
# permutation
# --------------------------------------------------------------------------

def permute_types(tree, rng):
    """Shuffle type labels across all mutations in the tree."""
    flat = array("b")
    for t in tree.types:
        flat.extend(t)
    shuffled = list(flat)
    rng.shuffle(shuffled)
    out = []
    k = 0
    for t in tree.types:
        n = len(t)
        out.append(array("b", shuffled[k:k + n]))
        k += n
    return out


def site_modal_types(tree):
    """Most common mutation type at each site."""
    per_site = defaultdict(Counter)
    for pos_arr, type_arr in zip(tree.positions, tree.types):
        for pos, t in zip(pos_arr, type_arr):
            per_site[pos][t] += 1
    return {pos: counter.most_common(1)[0][0] for pos, counter in per_site.items()}


def permute_sites(tree, rng, modal):
    """Shuffle the site-to-type assignment, keeping site positions in place."""
    sites = list(modal.keys())
    assigned = [modal[p] for p in sites]
    rng.shuffle(assigned)
    new_map = dict(zip(sites, assigned))
    out = []
    for pos_arr in tree.positions:
        out.append(array("b", [new_map[p] for p in pos_arr]))
    return out


# --------------------------------------------------------------------------
# the search, mirroring spectrumSplits.py find_splits()
# --------------------------------------------------------------------------

def own_spectra(tree, types, ncat, fold_map, max_branch_length):
    """Per-branch spectrum vectors, honouring --max_branch_length."""
    n = len(tree.node_ids)
    out = [None] * n
    zero = [0] * ncat
    for i in range(n):
        t = types[i]
        if not t or tree.raw_len[i] > max_branch_length:
            out[i] = zero
            continue
        vec = [0] * ncat
        if fold_map is None:
            for x in t:
                vec[x] += 1
        else:
            for x in t:
                vec[fold_map[x]] += 1
        out[i] = vec
    return out


def scan_subtree(tree, own, ncat, split_root, stop_set, min_mutations):
    """Score every candidate node in one subtree; return (max_chi2, node, below, above).

    Mirrors the inner loop of find_splits(): descend from split_root, skipping
    any child that is an accepted split, accumulate subtree spectra, and score
    each node's subtree against the remainder of the split root's spectrum.
    """
    # in-scope nodes, parents before children
    scope = []
    stack = [split_root]
    while stack:
        v = stack.pop()
        scope.append(v)
        for c in tree.children[v]:
            if c in stop_set:
                continue
            stack.append(c)

    local = {v: k for k, v in enumerate(scope)}
    spec = [list(own[v]) for v in scope]
    for k in range(len(scope) - 1, 0, -1):
        tgt = spec[local[tree.parent[scope[k]]]]
        src = spec[k]
        for j in range(ncat):
            tgt[j] += src[j]

    root_spec = spec[0]
    root_total = 0
    for v in root_spec:
        root_total += v

    max_chi = 0.0
    max_node = None
    max_below = 0
    for k in range(len(scope)):
        below = spec[k]
        n_below = 0
        for v in below:
            n_below += v
        if n_below < min_mutations:
            continue
        if root_total - n_below < min_mutations:
            continue
        above = [root_spec[j] - below[j] for j in range(ncat)]
        chi = chi2_2xk(below, above)
        if chi > max_chi:
            max_chi = chi
            max_node = scope[k]
            max_below = n_below
    return max_chi, max_node, max_below, root_total - max_below


def run_search(tree, own, ncat, min_chi, min_mutations, max_splits, log, label):
    """Greedy split search; returns the list of accepted splits as records."""
    accepted = {0}
    finalized = set()
    records = []
    iteration = 0
    while len(accepted) > len(finalized) and len(records) < max_splits:
        iteration += 1
        new_splits = set()
        for split_root in sorted(accepted):
            if split_root in finalized:
                continue
            chi, node, n_below, n_above = scan_subtree(
                tree, own, ncat, split_root, accepted, min_mutations)
            if node is not None and chi > min_chi and node not in accepted:
                new_splits.add(node)
                records.append({
                    "label": label,
                    "iteration": iteration,
                    "split_root_id": tree.node_ids[split_root],
                    "node_id": tree.node_ids[node],
                    "chi2": chi,
                    "mutations_below": n_below,
                    "mutations_above": n_above,
                })
                print(f"  [{label}] iter {iteration}: split at {tree.node_ids[node]} "
                      f"chi2={chi:.1f} below={n_below} above={n_above}", file=log)
                if len(records) >= max_splits:
                    break
            else:
                finalized.add(split_root)
                if node is None:
                    print(f"  [{label}] iter {iteration}: no eligible candidate under "
                          f"{tree.node_ids[split_root]} at --min_mutations "
                          f"{min_mutations}", file=log)
        accepted |= new_splits
        if not new_splits:
            break
    return records


# --------------------------------------------------------------------------
# replicate driver
# --------------------------------------------------------------------------

def one_replicate(tree, args, ncat, fold_map, modal, replicate, mode,
                  max_splits):
    """Score one replicate of one null; replicate 0 is the unpermuted tree."""
    if replicate == 0:
        types = tree.types
        label = "observed"
    else:
        rng = random.Random(args.seed + replicate)
        if mode == "types":
            types = permute_types(tree, rng)
        else:
            types = permute_sites(tree, rng, modal)
        label = f"perm{replicate}"
    own = own_spectra(tree, types, ncat, fold_map, args.max_branch_length)
    print(f"Scoring {mode} {label} ...", file=sys.stderr)
    return run_search(tree, own, ncat, args.min_chi, args.min_mutations,
                      max_splits, sys.stderr, label)


def quantile(sorted_values, p):
    if not sorted_values:
        return float("nan")
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    h = (n - 1) * p
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    return sorted_values[lo] + (h - lo) * (sorted_values[hi] - sorted_values[lo])


def summarize(maxima, obs, args, ncat, out):
    """The PERMUTATION NULL SUMMARY: each null's replicate maxima, sites
    against types, the recommended --min_chi from the --quantile of the
    sites maxima (types if that is all there is), and the observed maximum
    against it.  maxima maps mode -> sorted replicate maxima.  Returns the
    recommended --min_chi, or None."""
    q = args.quantile
    ql = f"p{100 * q:g}"
    print("=" * 72, file=out)
    print("PERMUTATION NULL SUMMARY", file=out)
    print("=" * 72, file=out)
    print(f"  categories {ncat}, --min_mutations {args.min_mutations}, "
          f"--max_branch_length {args.max_branch_length}", file=out)
    print(f"  types: whole tree only; sites: the cascade, --max_splits "
          f"{args.max_splits}, --min_chi {args.min_chi:g}", file=out)
    print(file=out)
    print("  Maximum chi-squared per null replicate (each replicate's maximum "
          "is", file=out)
    print("  in practice its whole-tree split):", file=out)
    print(f"    null	replicates	min	median	{ql}	max", file=out)
    for mode in ("types", "sites"):
        m = maxima.get(mode)
        if m is None:
            continue
        if not m:
            print(f"    {mode}	0	-	-	-	-", file=out)
            continue
        print(f"    {mode}	{len(m)}	{m[0]:.1f}	{quantile(m, 0.5):.1f}	"
              f"{quantile(m, q):.1f}	{m[-1]:.1f}", file=out)
    basis = "sites" if maxima.get("sites") else "types"
    m = maxima.get(basis) or []
    if maxima.get("sites") and maxima.get("types"):
        s_q = quantile(maxima["sites"], q)
        t_q = quantile(maxima["types"], q)
        above = sum(1 for v in maxima["sites"] if v > maxima["types"][-1])
        print(f"  sites {ql} is {s_q / t_q:.2f}x types {ql}; "
              f"{above} of {len(maxima['sites'])} sites replicates "
              f"{'exceeds' if above == 1 else 'exceed'} the largest types "
              f"one.", file=out)
        print("  sites keeps the clumping recurrent sites cause, so it comes "
              "out at or above", file=out)
        print("  types; the gap is what residual recurrent-site clumping "
              "costs.  Calibrate", file=out)
        print("  against sites.", file=out)
    print(file=out)
    if not m:
        print("  No null replicate produced a scoreable split.  Either "
              "--min_mutations", file=out)
        print("  excludes every candidate, or --min_chi is above the null "
              "range.", file=out)
        print(file=out)
        return None
    floor = quantile(m, q)
    rec = math.ceil(args.margin * floor)
    print(f"  In {100 * (1 - q):g}% of {basis} replicates the null search "
          f"reached chi-squared {floor:.1f}", file=out)
    print("  or more with no real clade-specific spectrum present.  "
          "--min_chi must exceed it:", file=out)
    print(f"    bare floor         --min_chi {math.ceil(floor)}", file=out)
    print(f"    with 1.5x margin   --min_chi {math.ceil(1.5 * floor)}",
          file=out)
    print(f"    with 2x margin     --min_chi {math.ceil(2.0 * floor)}",
          file=out)
    print(f"  Recommended (--margin {args.margin:g}): --min_chi {rec}",
          file=out)
    need = math.ceil(5 / (1 - q)) if q < 1 else 0
    if len(m) < need:
        print(f"  With {len(m)} replicates the {ql} rests on the top one or "
              f"two; about {need}", file=out)
        print("  resolve it.  Raise --replicates.", file=out)
    if obs is not None:
        print(file=out)
        print(f"  Observed (unpermuted) maximum chi-squared: {obs:.1f}",
              file=out)
        ratio = obs / floor if floor > 0 else float("inf")
        print(f"  Ratio to the {basis} {ql}: {ratio:.1f}x (to its largest "
              f"replicate: "
              f"{obs / m[-1] if m[-1] > 0 else float('inf'):.1f}x)", file=out)
        if ratio < 2:
            print("  That is a thin margin: the strongest real split is close "
                  "to", file=out)
            print("  what chance alone produces.  Treat any split here as",
                  file=out)
            print("  provisional and check it against the bootstrap.",
                  file=out)
    print(file=out)
    return rec



# --------------------------------------------------------------------------
# size-matched null and hotspot check (sections 2-4)
# --------------------------------------------------------------------------

class TreeView(object):
    """The node interface the hotspot check walks (id, parent, children,
    mutations, get_node, root), over a FlatTree, so it needs no second read
    of the tree.  Mutations are rebuilt from the parsed SNVs, which is all
    the check counts."""

    class Node(object):
        __slots__ = ("t", "i")

        def __init__(self, t, i):
            self.t, self.i = t, i

        @property
        def id(self):
            return self.t.flat.node_ids[self.i]

        @property
        def parent(self):
            p = self.t.flat.parent[self.i]
            return TreeView.Node(self.t, p) if p >= 0 else None

        @property
        def children(self):
            return [TreeView.Node(self.t, c) for c in self.t.flat.children[self.i]]

        @property
        def mutations(self):
            f = self.t.flat
            return [f"{MUT_TYPES[ty][0]}{p}{MUT_TYPES[ty][1]}"
                    for p, ty in zip(f.positions[self.i], f.types[self.i])]

    def __init__(self, flat):
        self.flat = flat
        self.index = {nid: i for i, nid in enumerate(flat.node_ids)}
        self.root = TreeView.Node(self, 0)

    def get_node(self, nid):
        return TreeView.Node(self, self.index[nid])


def rows_from_records(records):
    """run_search records as the rows the size report reads."""
    return [{"replicate": r["label"], "iteration": r["iteration"],
             "root": r["split_root_id"], "node": r["node_id"],
             "chi2": r["chi2"], "below": r["mutations_below"],
             "size": r["mutations_below"] + r["mutations_above"]}
            for r in records]


def parse_cascade(lines):
    """All rows as dicts, observed and null alike."""
    rows = []
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if i == 0 or len(f) < 7:
            continue
        rows.append({
            "replicate": f[0], "iteration": int(f[1]), "root": f[2],
            "node": f[3], "chi2": float(f[4]), "below": int(f[5]),
            "size": int(f[5]) + int(f[6]),
        })
    return rows


def size_bins(null_rows):
    """Group null rows into factor-of-2 bins, counted down from the largest."""
    if not null_rows:
        return []
    top = max(r["size"] for r in null_rows)
    bins = collections.defaultdict(list)
    for r in null_rows:
        k = int(math.floor(math.log2(top / r["size"]))) if r["size"] > 0 else 99
        bins[k].append(r)
    return [bins[k] for k in sorted(bins)]


def lowest_min_chi(whole):
    """The smallest integer --min_chi above every whole-tree null maximum."""
    return math.floor(max(r["chi2"] for r in whole)) + 1 if whole else None


def fmt_size(n):
    if n >= 1e6:
        return f"{n / 1e6:.2f}M"
    if n >= 1e3:
        return f"{n / 1e3:.0f}k"
    return str(n)


def same_clades(tree, first, same):
    """node -> the clade it counts as: itself, or the nearest ancestor among
    the whole-tree maxima, if its mutations_below is at least `same` times
    that ancestor's (a node and a child holding nearly all of it)."""
    below = {r["node"]: r["below"] for r in first}
    clade = {}
    for node in sorted(below, key=lambda n: -below[n]):
        clade[node] = node
        if tree is None:
            continue
        p = tree.get_node(node).parent
        while p is not None:
            if p.id in below:
                if below[node] >= same * below[p.id]:
                    clade[node] = clade[p.id]
                break
            p = p.parent
    return clade


def hotspot(null_rows, fraction, top_maxima, clade=None):
    """The whole-tree maxima, ranked, and the clades flagged, with why."""
    first = [r for r in null_rows if r["iteration"] == 1]
    clade = clade or {}

    def key(r):
        return clade.get(r["node"], r["node"])

    counts = collections.Counter(key(r) for r in first)
    n = len({r["replicate"] for r in first})
    ranked = sorted(first, key=lambda r: -r["chi2"])
    flagged = {}
    for c, k in counts.most_common():
        if n and k / n > fraction:
            flagged[c] = [f"it is the maximum in {k} of {n} replicates"]
    if top_maxima > 1 and len(ranked) >= top_maxima:
        lead = key(ranked[0])
        if all(key(r) == lead for r in ranked[:top_maxima]):
            rest = [r["chi2"] for r in ranked if key(r) != lead]
            why = (f"it holds the {top_maxima} largest maxima ("
                   + ", ".join(f"{r['chi2']:.1f}" for r in ranked[:top_maxima])
                   + ")" + (f"; no other clade's exceeds {max(rest):.1f}"
                            if rest else ""))
            flagged.setdefault(lead, []).append(why)
    return {"counts": counts, "n": n, "ranked": ranked, "flagged": flagged,
            "key": key}


def subtree(tree, root=None):
    """Every node at and below root (the whole tree if None)."""
    stack = [tree.get_node(root) if root else tree.root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(node.children)


def site_counts(nodes):
    counts = collections.Counter()
    for node in nodes:
        for m in node.mutations:
            if len(m) > 2 and m[0] in "ACGT" and m[-1] in "ACGT":
                counts[int(m[1:-1])] += 1
    return counts


def clade_sites(tree, node, top, whole):
    """(clade's mutations, [(position, in clade, tree-wide)]) for its most
    recurrent sites.  The clade includes the node's own branch, as
    mutations_below does."""
    inside = site_counts(subtree(tree, node))
    return sum(inside.values()), [(pos, k, whole[pos])
                                  for pos, k in inside.most_common(top)]


def size_report(rows, args, out, tree=None, flag_chi=None):
    """Sections 2-4 from the sites cascade rows; flag_chi is the --min_chi the
    observed splits are flagged against."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from split_confidence import ceiling
    null_rows = [r for r in rows if r["replicate"] != "observed"]
    observed = [r for r in rows if r["replicate"] == "observed"]
    nulls = [(r["size"], r["chi2"]) for r in null_rows]
    n_reps = len({r["replicate"] for r in null_rows})

    print("=" * 72, file=out)
    print("NULL CEILING BY SUBTREE SIZE", file=out)
    print("=" * 72, file=out)
    print(f"  {len(null_rows)} null splits from {n_reps} replicates, in factor-of-2 "
          "bins of the", file=out)
    print("  size of the table each was scored on (mutations_below + "
          "mutations_above).", file=out)
    print("  The whole-tree pass (iteration 1) is listed first; its maxima "
          "are the", file=out)
    print("  summary's sites row.", file=out)
    ql = f"p{100 * args.null_quantile:g}"
    whole = [r for r in null_rows if r["iteration"] == 1]
    deeper = [r for r in null_rows if r["iteration"] > 1]
    lowest = lowest_min_chi(whole)
    mcol = f"\tmargin at {lowest}" if lowest else ""
    print(f"    size range\tnull splits\tmedian chi2\t{ql} chi2\tmax chi2\t"
          f"max at{mcol}", file=out)
    groups = ([("whole tree", whole)] if whole else []) + [(None, b) for b in size_bins(deeper)]
    for label, b in groups:
        lo = min(r["size"] for r in b)
        hi = max(r["size"] for r in b)
        top = max(b, key=lambda r: r["chi2"])
        rng = fmt_size(lo) if fmt_size(lo) == fmt_size(hi) else f"{fmt_size(lo)}-{fmt_size(hi)}"
        if label:
            rng = f"{rng} ({label})"
        q = quantile(sorted(r["chi2"] for r in b), args.null_quantile)
        m = f"\t{lowest / q:.1f}x" if lowest and q > 0 else ""
        print(f"    {rng}\t{len(b)}\t"
              f"{statistics.median(r['chi2'] for r in b):.1f}\t{q:.1f}\t"
              f"{top['chi2']:.1f}\t{top['replicate']} {top['node']}{m}", file=out)
    if lowest:
        top = max(r["chi2"] for r in whole)
        print(f"  Lowest --min_chi that clears the whole-tree null: {lowest} "
              f"(its largest", file=out)
        print(f"  replicate, {top:.1f}, rounded up), so no null replicate's "
              "first split passes.", file=out)
        print(f"  The last column is that threshold's margin over each bin's "
              f"{ql}: where it", file=out)
        print(f"  is {args.margin:g}x or more, a split just at the threshold "
              "clears its size-matched", file=out)
        print("  null by --margin.  Lower in the table means smaller subtrees; "
              "near the root", file=out)
        print("  the margin is thin, but real splits there score far above "
              "the threshold,", file=out)
        print("  and each split's own margin is checked against its size "
              "anyway.", file=out)
    print(file=out)

    if observed:
        print("=" * 72, file=out)
        print("OBSERVED SPLITS AGAINST A SIZE-MATCHED NULL", file=out)
        print("=" * 72, file=out)
        print(f"  ceiling = {ql} of null chi2 from subtrees within "
              f"{args.null_window}x of the split's", file=out)
        print("  size (split_confidence.py's rule); margin = chi2 / ceiling.", file=out)
        if flag_chi is not None:
            print(f"  flag: 'rejected' if chi2 <= --min_chi {flag_chi:g} but "
                  f"margin >= {args.margin:g};", file=out)
            print("        'within null' if chi2 > --min_chi but margin < 1.", file=out)
        print(f"  ('few' marks a ceiling from fewer than {args.min_nulls} null splits.)",
              file=out)
        print("    iter\tsplit_root\tnode\tsize\tchi2\tceiling\tnulls\tmargin\tflag",
              file=out)
        n_rejected = n_within = n_fail = 0
        for r in observed:
            ceil, n = ceiling(nulls, r["size"], args.null_window,
                              args.null_quantile)
            flags = []
            if ceil is None:
                margin_s, ceil_s = "-", "-"
                margin = None
            else:
                margin = r["chi2"] / ceil
                margin_s, ceil_s = f"{margin:.1f}", f"{ceil:.1f}"
                if margin < 1:
                    n_fail += 1
            if n < args.min_nulls:
                flags.append("few")
            if flag_chi is not None and margin is not None:
                if r["chi2"] <= flag_chi and margin >= args.margin:
                    flags.append("rejected")
                    n_rejected += 1
                elif r["chi2"] > flag_chi and margin < 1:
                    flags.append("within null")
                    n_within += 1
            print(f"    {r['iteration']}\t{r['root']}\t{r['node']}\t{r['size']}\t"
                  f"{r['chi2']:.1f}\t{ceil_s}\t{n}\t{margin_s}\t{','.join(flags)}",
                  file=out)
        print(f"  {n_fail} of {len(observed)} observed splits fall below their "
              "size-matched ceiling.", file=out)
        if flag_chi is not None:
            print(f"  At --min_chi {flag_chi:g}: {n_rejected} rejected splits clear "
                  f"their null by >= {args.margin:g}x;", file=out)
            print(f"  {n_within} accepted splits do not clear their null.", file=out)
        print(file=out)

    print("=" * 72, file=out)
    print("HOTSPOT CHECK", file=out)
    print("=" * 72, file=out)
    first = [r for r in null_rows if r["iteration"] == 1]
    clade = same_clades(tree, first, args.same_clade)
    h = hotspot(null_rows, args.hotspot_fraction, args.top_maxima, clade)
    key, n = h["key"], h["n"]
    merged = tree is not None and any(c != node for node, c in clade.items())
    print(f"  The largest whole-tree null maxima (iteration 1), of {n} "
          "replicates:", file=out)
    print("    chi2\treplicate\tnode\tmutations_below\tclade", file=out)
    for r in h["ranked"][:5]:
        print(f"    {r['chi2']:.1f}\t{r['replicate']}\t{r['node']}\t"
              f"{r['below']}\t{key(r)}", file=out)
    print("  Replicates in which each clade is the maximum:", file=out)
    print("    clade\treplicates", file=out)
    for c, k in h["counts"].most_common(5):
        print(f"    {c}\t{k}", file=out)
    if merged:
        print(f"  (Nested nodes whose mutations_below agree within "
              f"{1 - args.same_clade:.0%} count as one clade.)", file=out)
    print(f"  Flags: a clade that is the maximum in more than "
          f"{args.hotspot_fraction:.0%} of replicates", file=out)
    print(f"  (--hotspot_fraction), or that holds the {args.top_maxima} "
          "largest maxima (--top_maxima).", file=out)
    if not h["flagged"]:
        print("  Verdict: no hotspot.", file=out)
    whole = site_counts(subtree(tree)) if tree is not None and h["flagged"] \
        else None
    for c, why in h["flagged"].items():
        print(f"  Verdict: HOTSPOT at {c}: " + "; ".join(why) + ".", file=out)
        print("  A recurrent site inside that clade is likely setting the "
              "threshold.", file=out)
        if tree is None:
            continue
        total, sites = clade_sites(tree, c, args.hotspot_sites, whole)
        print(f"  {c}'s most recurrent sites ({total} mutations in the "
              "clade):", file=out)
        print("    position\tin clade\ttree-wide\tshare of the clade", file=out)
        for pos, k, w in sites:
            print(f"    {pos}\t{k}\t{w}\t{100 * k / total:.2f}%", file=out)
        hot = [str(pos) for pos, k, w in sites if k >= args.hotspot_min]
        if hot:
            print(f"  Sites with >= {args.hotspot_min} mutations in the clade: "
                  "annotate them, and mask them", file=out)
            print("  if they look like an artifact (step 7):", file=out)
            cmd = "annotate_sites.py"
            if args.annotation:
                cmd += f" --annotation {args.annotation}"
            if args.reference:
                cmd += f" --reference {args.reference}"
            print(f"    {cmd} --positions {','.join(hot)}", file=out)
    print(file=out)


def _selftest_size():
    import io
    lines = ["replicate\titeration\tsplit_root_id\tnode_id\tchi2\t"
             "mutations_below\tmutations_above\n"]
    # Whole tree 1000: observed splits big; null maxima, node_7 wins 3 of 4.
    lines.append("observed\t1\tnode_1\tnode_5\t500\t300\t700\n")
    lines.append("observed\t2\tnode_1\tnode_9\t30\t100\t500\n")      # size 600
    lines.append("observed\t2\tnode_5\tnode_6\t8\t100\t200\n")       # size 300
    for rep, node, chi in (("perm1", "node_7", 40), ("perm2", "node_7", 35),
                           ("perm3", "node_7", 38), ("perm4", "node_8", 33)):
        lines.append(f"{rep}\t1\tnode_1\t{node}\t{chi}\t200\t800\n")
        lines.append(f"{rep}\t2\tnode_1\tnode_x\t{chi / 3:.1f}\t150\t550\n")
        lines.append(f"{rep}\t2\tnode_7\tnode_y\t{chi / 4:.1f}\t50\t250\n")
    rows = parse_cascade(lines)
    assert len(rows) == 15
    null_rows = [r for r in rows if r["replicate"] != "observed"]
    bins = size_bins([r for r in null_rows if r["iteration"] > 1])
    assert [len(b) for b in bins] == [4, 4], [len(b) for b in bins]
    h = hotspot(null_rows, 0.5, 2)
    assert h["n"] == 4 and list(h["flagged"]) == ["node_7"], h["flagged"]
    assert len(h["flagged"]["node_7"]) == 2, h["flagged"]   # both signs
    h = hotspot(null_rows, 0.75, 5)
    assert h["flagged"] == {}, h["flagged"]
    args = argparse.Namespace(null_window=1.5, null_quantile=1.0, min_chi=50,
                              margin=2, min_nulls=5,
                              hotspot_fraction=0.5, top_maxima=2,
                              same_clade=0.95, hotspot_sites=3, hotspot_min=3,
                              annotation="ann.gff3", reference="ref.fa")
    buf = io.StringIO()
    size_report(rows, args, buf, flag_chi=50)
    text = buf.getvalue()
    # node_9 at 600: the 1.5x window (400-900) holds only the size-700 nulls,
    # whose max is 40/3 -> margin 2.25, below --min_chi 50, so rejected.
    line9 = next(l for l in text.splitlines() if "\tnode_9\t" in l)
    assert line9.split("\t")[5] == "13.3" and "rejected" in line9, line9
    # node_6 at 300: ceiling 10.0 (40/4), chi 8 -> fails its null.
    line6 = next(l for l in text.splitlines() if "\tnode_6\t" in l)
    assert line6.split("\t")[7] == "0.8", line6
    assert "HOTSPOT at node_7: it is the maximum in 3 of 4" in text, text
    assert "1k (whole tree)\t4\t36.5\t40.0\t40.0\tperm1 node_7\t1.0x" in text, text
    # 41 clears the largest whole-tree null, 40.0; it is 41/13.3 = 3.1x the
    # 700-mutation bin's ceiling, and 4.1x the 300-mutation bin's.
    assert "Lowest --min_chi that clears the whole-tree null: 41" in text, text
    assert "700\t4\t12.2\t13.3\t13.3\tperm1 node_x\t3.1x" in text, text
    assert "300\t4\t9.2\t10.0\t10.0\tperm1 node_y\t4.1x" in text, text

    # A clade that rarely wins but wins big, under two IDs: node_20 and its
    # child node_21, which holds nearly all of it.  It is the maximum in 2 of
    # 6 replicates, under the frequency limit, but holds the 2 largest.
    # Two sites inside it carry most of its mutations.
    class Node:
        def __init__(self, nid, muts, children=()):
            self.id, self.mutations = nid, list(muts)
            self.children, self.parent = list(children), None
            for c in self.children:
                c.parent = self

    class Tree:
        def __init__(self, root):
            self.root, self.nodes, stack = root, {}, [root]
            while stack:
                n = stack.pop()
                self.nodes[n.id] = n
                stack.extend(n.children)

        def get_node(self, nid):
            return self.nodes[nid]

    leaves = [Node(f"leaf{i}", ["A50G", "T52G"] if i % 2 else ["G50A"])
              for i in range(8)]
    n21 = Node("node_21", ["C9T"], leaves)
    n20 = Node("node_20", ["C7T"], [n21])
    tree = Tree(Node("node_1", [], [n20, Node("node_30", ["A50G", "G3T"])]
                     + [Node(f"node_{i}", ["G3T"]) for i in (31, 32, 33)]))
    lines = lines[:1] + [f"perm{i}\t1\tnode_1\t{node}\t{chi}\t{below}\t900\n"
                         for i, (node, chi, below) in enumerate(
                             [("node_20", 45, 100), ("node_21", 31, 98),
                              ("node_30", 24, 300), ("node_31", 22, 200),
                              ("node_32", 21, 200), ("node_33", 20, 200)])]
    rows = parse_cascade(lines)
    h = hotspot(rows, 0.5, 2)
    assert h["flagged"] == {}, "without the tree, two IDs hide the clade"
    first = [r for r in rows if r["iteration"] == 1]
    clade = same_clades(tree, first, 0.95)
    assert clade["node_21"] == "node_20" and clade["node_30"] == "node_30"
    h = hotspot(rows, 0.5, 2, clade)
    assert list(h["flagged"]) == ["node_20"], h["flagged"]
    assert "no other clade's exceeds 24.0" in h["flagged"]["node_20"][0]
    total, sites = clade_sites(tree, "node_20", 3,
                               site_counts(subtree(tree)))
    assert total == 14 and sites[0] == (50, 8, 9) and sites[1] == (52, 4, 4), \
        (total, sites)
    buf = io.StringIO()
    size_report(rows, args, buf, tree, flag_chi=50)
    text = buf.getvalue()
    assert "HOTSPOT at node_20: it holds the 2 largest maxima (45.0, 31.0)" \
        in text, text
    assert ("annotate_sites.py --annotation ann.gff3 --reference ref.fa "
            "--positions 50,52") in text, text
    print("ok: size bins, size-matched ceilings, --min_chi flags and the hotspot "
          "check behave on a constructed cascade")
    print("ok: a clade under two IDs that holds the largest maxima is flagged, "
          "with its recurrent sites")
    ### The same clade logic on the extracted flat tree, through TreeView.
    class N(object):
        def __init__(self, nid, muts, kids=()):
            self.id, self.mutations, self.children = nid, list(muts), list(kids)
    flat = extract(N("node_1", [], [
        N("node_20", ["C7T"], [N("node_21", ["C9T"], [
            N(f"leaf{i}", ["A50G", "T52G"] if i % 2 else ["G50A"])
            for i in range(8)])]),
        N("node_30", ["A50G", "G3T"])]
        + [N(f"node_{i}", ["G3T"]) for i in (31, 32, 33)]))
    view = TreeView(flat)
    assert same_clades(view, first, 0.95) == clade
    assert clade_sites(view, "node_20", 3, site_counts(subtree(view))) == \
        (total, sites)
    print("ok: the hotspot check reads the extracted tree as it read bte's")



def main(argv=None):
    args = parse_args(argv)

    if args.selftest:
        return selftest()

    ### Not enforced by argparse so that --selftest runs without a tree.
    if not args.input_tree:
        print("error: --input_tree is required (or pass --selftest)",
              file=sys.stderr)
        return 2

    fold_map, fold_labels = folded_index()
    if args.fold:
        ncat = 6
        cats = fold_labels
    else:
        ncat = 12
        fold_map = None
        cats = MUT_TYPES

    try:
        import bte
    except ImportError:
        print("error: could not import bte; run this in the environment you use "
              "for spectrumSplits.py", file=sys.stderr)
        return 1

    print(f"Reading {args.input_tree} ...", file=sys.stderr)
    mat = bte.MATree(args.input_tree)
    tree = extract(mat.root)
    del mat
    n_nodes = len(tree.node_ids)
    n_muts = sum(len(t) for t in tree.types)
    print(f"Extracted {n_nodes} nodes, {n_muts} usable SNVs, "
          f"scoring {ncat} categories ({', '.join(cats)})", file=sys.stderr)

    modes = ["types", "sites"] if args.mode == "both" else [args.mode]
    modal = site_modal_types(tree) if "sites" in modes else None
    if modal is not None:
        print(f"Site-level null over {len(modal)} distinct sites", file=sys.stderr)

    ### Replicate 0 is the unpermuted tree, scored through the same code.
    replicates = list(range(0, args.replicates + 1))
    use_fork = args.nthreads > 1
    if use_fork:
        try:
            import multiprocessing
            ctx = multiprocessing.get_context("fork")
        except (ImportError, ValueError):
            print("warning: fork start method unavailable; running serially. "
                  "Run several instances with different --seed to parallelize.",
                  file=sys.stderr)
            use_fork = False
    pool = ctx.Pool(processes=args.nthreads) if use_fork else None

    records, paths = {}, []
    try:
        for mode in modes:
            depth = 1 if mode == "types" else args.max_splits
            results = []
            if pool is not None:
                jobs = [pool.apply_async(one_replicate,
                                         (tree, args, ncat, fold_map, modal,
                                          r, mode, depth))
                        for r in replicates]
                for job in jobs:
                    results.extend(job.get())
            else:
                for r in replicates:
                    results.extend(one_replicate(tree, args, ncat, fold_map,
                                                 modal, r, mode, depth))
            records[mode] = results
            name = "types" if mode == "types" else \
                ("sites_cascade" if depth > 1 else "sites")
            paths.append(f"{args.prefix}_{name}.tsv")
            with open(paths[-1], "w") as fh:
                fh.write("\t".join(["replicate", "iteration", "split_root_id",
                                    "node_id", "chi2", "mutations_below",
                                    "mutations_above"]) + "\n")
                for rec in results:
                    fh.write("\t".join([
                        rec["label"], str(rec["iteration"]),
                        rec["split_root_id"], rec["node_id"],
                        f"{rec['chi2']:.4f}", str(rec["mutations_below"]),
                        str(rec["mutations_above"]),
                    ]) + "\n")
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    print(f"\nPer-split records written to {', '.join(paths)}",
          file=sys.stderr)

    # ---- the report ------------------------------------------------------
    out = sys.stdout
    maxima, obs = {}, None
    for mode, results in records.items():
        per_label = defaultdict(list)
        for rec in results:
            per_label[rec["label"]].append(rec["chi2"])
        maxima[mode] = sorted(max(v) for k, v in per_label.items()
                              if k != "observed" and v)
        if per_label.get("observed"):
            obs = max(per_label["observed"])
    rec_chi = summarize(maxima, obs, args, ncat, out)
    if "sites" in records:
        flag = args.flag_min_chi if args.flag_min_chi is not None else rec_chi
        if args.max_splits == 1:
            print("  (--max_splits 1: no cascade, so the size bins below hold "
                  "only the whole tree.)", file=out)
        size_report(rows_from_records(records["sites"]), args, out,
                    TreeView(tree), flag)
    print("  Next: run spectrumSplits.py with the recommended --min_chi "
          "(step 8), unless the", file=out)
    print("  hotspot check flagged a clade: then mask its sites and rerun "
          "this first.", file=out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
