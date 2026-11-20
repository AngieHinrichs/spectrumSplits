#!/usr/bin/env python3
"""Does trinucleotide context separate a split beyond its class shift?

spectrumSplits.py scores a spectrum of mutation classes: 6 strand-folded
classes with --fold, 12 forward-strand types without.  Two clades with the
same class proportions can still differ inside a class -- more C>T at CpG and
less elsewhere, say -- and a split found on a class shift may or may not
carry a context shift too.  This script asks that of every accepted split,
after the fact, without touching the search.

Match the search: with --fold, the 96 folded context classes (6 classes x 16
flank pairs); without it, the 192 forward-strand ones (12 types x 16), which
keep apart the two strands an unfolded search keeps apart.

For each split, its residual row (its mutations minus those of splits nested
in it, as spectra.tsv reports) is set against its comparator's residual row as
a 2 x 96 (or 2 x 192) table of contexts.  The likelihood-ratio statistic G of
that table splits exactly in two:

  G(contexts) = G(between)  +  G(within)

  G(between)  the 2 x 6 (or 2 x 12) table of classes: the shift
              spectrumSplits scored.
  G(within)   the sum over the classes of a 2 x 16 table of that class's
              contexts.  It is zero when the two rows have the same context
              shares inside every class, whatever their class proportions.

G(within) is the test.  Its asymptotic chi-squared p value is not used:
mutations are not independent draws (recurrent sites, branches whose
mutations share a cause), so the null comes from shuffling whole branches
between the split and its comparator, each branch keeping all its mutations,
and recomputing.  p = (1 + shuffles at or above the observed) / (1 +
shuffles).  The same shuffles give a p value per class.  "excess" is the
observed G(within) over the mean of the shuffles -- about 1 under the null,
and comparable between splits of different size, which G itself is not.
q is Benjamini-Hochberg over the splits tested.

The comparator.  By default, the split's clade parent: the nearest split root
above it, i.e. the residual row it was carved out of in the final
decomposition.  --against search uses the subtree the search found it in
(from --splits_log), which is what split_drivers.py compares against.

Contexts come from the parent's genome (trinucleotide.call_contexts), so
--reference must be the ROOT's genome -- the rerooted reference on a rerooted
tree.  Mutations with an unknown or masked flank, and both halves of an
adjacent same-branch pair, are left out of every table, including G(between),
so between + within = G(contexts) exactly.  --mask should be the masks the tree was
built with.

Two runs, one invocation.  Everything above is done twice: on all branches,
and on internal branches only, which keeps out per-genome error (as in step
9's branch_class_spectra.py).  Contexts are called once; each run has its own
shuffles and controls, from the same --seed.

ctrl_q.  ctrl_p is one test per split, so it is corrected for the splits
tested like the shuffle's p: ctrl_q is its Benjamini-Hochberg value within
the run.  Splits share controls, so their ctrl_p are positively dependent,
under which BH still holds.  ctrl_p can't fall below 1/(controls in the
split's size window + 1), and a split with few controls in its window may be
unable to reach --control_alpha whatever its excess.  Each run names those,
and says which fix applies: more --controls, or, where every clade of that
size on the tree has already been drawn, a wider --control_window.

The summary.  A split whose context shift is worth reading passes both nulls
-- q below --max_q and ctrl_q at or below --control_alpha -- on BOTH runs.
The report ends with those splits, sorted by all-branch excess (--summary_top
N keeps the N with the lowest all-branch ctrl_p, if you want fewer), each
with its clades (from --composition, split_composition.py's output), excess,
ctrl_p, ctrl_q and q on both runs, the --summary_contexts contexts carrying most of its all-branch
G(within) with their shares of class, and each context's top-site share.  A
listed context where one site carries --one_site or more of the split's count
is flagged: a hotspot standing in for a context.

Outputs, named from --prefix (default split_contexts):
  PREFIX.tsv, PREFIX_internal.tsv   one row per split, per run
  PREFIX_96.tsv, PREFIX_internal_96.tsv
                                    one row per split and context, per run
                                    (_192 without --fold)
  PREFIX_table.tsv                  the summary table
The report goes to stdout.

What it cannot tell you: whether a search scoring contexts would find splits
this one did not.  That needs the search to score contexts.

Usage:
    split_contexts.py --input_tree tree.pruned.masked.pb.gz \\
        --spectra spectra.tsv --fold --max_branch_length 250 \\
        --reference root_genome.fa --mask selection_loci.bed,alignment_loci.bed \\
        --composition split_composition.tsv --summary_top 12 \\
        > split_contexts.txt

    split_contexts.py --selftest
"""

import argparse
import csv
import gzip
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from permutation_null import extract
from branch_class_spectra import residual_groups
import trinucleotide as tn

PER = 16


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Test each split for a trinucleotide-context shift "
                    "within its mutation classes.")
    p.add_argument("--input_tree", help="The tree the search ran on (protobuf)")
    p.add_argument("--spectra",
                   help="spectrumSplits.py --output_spectrum TSV; split roots "
                        "and residual totals are read from it")
    p.add_argument("--fold", action="store_true",
                   help="Use the 96 folded context classes, for a search run "
                        "with --fold; without it, the 192 forward-strand "
                        "classes of an unfolded search")
    p.add_argument("--max_branch_length", type=int, default=None,
                   help="Match the search; branches above it contribute "
                        "nothing")
    p.add_argument("--reference",
                   help="FASTA of the tree's ROOT genome (on a rerooted tree, "
                        "the rerooted reference)")
    p.add_argument("--mask", default=None,
                   help="Comma-separated BED files of masked regions; a "
                        "context with a masked flank is not called")
    p.add_argument("--against", choices=["clade", "search"], default="clade",
                   help="Comparator: the nearest enclosing split (clade, "
                        "default) or the subtree the search found it in "
                        "(search, needs --splits_log)")
    p.add_argument("--splits_log", default=None,
                   help="The search's stderr, for --against search")
    p.add_argument("--permutations", type=int, default=1000,
                   help="Branch shuffles per split (default 1000)")
    p.add_argument("--controls", type=int, default=30,
                   help="Random non-split clades drawn per split, within "
                        "--control_window of its size, as a reference for "
                        "how much any clade differs (default 30; 0 skips)")
    p.add_argument("--control_window", type=float, default=2.0,
                   help="Fold range of size a control may differ from a "
                        "split by (default 2)")
    p.add_argument("--control_permutations", type=int, default=200,
                   help="Branch shuffles per control clade (default 200)")
    p.add_argument("--control_alpha", type=float, default=0.05,
                   help="ctrl_q (BH-corrected ctrl_p) a split must be at or "
                        "below to pass the random-clade null, for the "
                        "detailed listing and the summary (default 0.05)")
    p.add_argument("--max_q", type=float, default=0.05,
                   help="Branch-shuffle q a split must be below to pass "
                        "(default 0.05)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--top", type=int, default=3,
                   help="Contexts to list per split in each run's detailed "
                        "listing (default 3)")
    p.add_argument("--composition", default=None,
                   help="split_composition.py --output, for each summarised "
                        "split's clades")
    p.add_argument("--summary_top", type=int, default=0,
                   help="Summarise only the N splits passing both runs with "
                        "the lowest all-branch ctrl_p (default 0: all)")
    p.add_argument("--summary_contexts", type=int, default=2,
                   help="Contexts per split in the summary (default 2)")
    p.add_argument("--one_site", type=float, default=0.2,
                   help="Top-site share at which a context is flagged as one "
                        "site (default 0.2)")
    p.add_argument("--min_class", type=int, default=80,
                   help="A class with fewer called mutations than this in "
                        "either row is marked 'few' in the report "
                        "(default 80, five per context)")
    p.add_argument("--prefix", default="split_contexts",
                   help="Prefix of the output TSVs (default split_contexts)")
    p.add_argument("--selftest", action="store_true",
                   help="Run the self-tests and exit")
    return p.parse_args(argv)


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------

def _xlogx(np, a):
    a = np.asarray(a, dtype=np.float64)
    out = np.zeros_like(a)
    pos = a > 0
    out[pos] = a[pos] * np.log(a[pos])
    return out


def g_parts(np, A, B):
    """G statistics of 2 x n tables, rows A and B (arrays of shape (..., n)),
    n = 16 x classes (96 or 192).

    Returns (between, within_by_class) with shapes (...) and (..., classes).
    G = 2 * [sum O ln O - sum rows ln rows - sum cols ln cols + N ln N],
    so between + within.sum(-1) equals the 2 x n G exactly.
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    shp = A.shape[:-1]
    ncls = A.shape[-1] // PER
    A6 = A.reshape(shp + (ncls, PER))
    B6 = B.reshape(shp + (ncls, PER))
    C6 = A6 + B6
    rA = A6.sum(-1)
    rB = B6.sum(-1)
    col = C6.sum(-1)
    within = 2 * (_xlogx(np, A6).sum(-1) + _xlogx(np, B6).sum(-1)
                  - _xlogx(np, rA) - _xlogx(np, rB)
                  - _xlogx(np, C6).sum(-1) + _xlogx(np, col))
    nA = rA.sum(-1)
    nB = rB.sum(-1)
    n = nA + nB
    between = 2 * (_xlogx(np, rA).sum(-1) + _xlogx(np, rB).sum(-1)
                   - _xlogx(np, nA) - _xlogx(np, nB)
                   - _xlogx(np, col).sum(-1) + _xlogx(np, n))
    ### Rounding can leave -1e-12 on an identical pair of rows.
    return np.maximum(between, 0.0), np.maximum(within, 0.0)


def shuffle_test(np, Xs, Xc, permutations, rng, chunk=256):
    """Branch-shuffle null for G(within).

    Xs, Xc: (branches, contexts) counts of the split's and the comparator's
    branches; branches with no called context are dropped by the caller.
    Each shuffle reassigns which branches belong to the split, keeping the
    number of split branches, and recomputes G(within) per class.

    Returns dict with obs (classes,), null (permutations, classes).
    """
    X = np.vstack([Xs, Xc]).astype(np.float32)
    ns, n = Xs.shape[0], X.shape[0]
    T = X.sum(0, dtype=np.float64)
    A = Xs.sum(0, dtype=np.float64)
    _, obs = g_parts(np, A, T - A)
    null = np.empty((permutations, X.shape[1] // PER))
    done = 0
    while done < permutations:
        c = min(chunk, permutations - done)
        if ns * 2 <= n:
            pick = np.argpartition(rng.random((c, n)), ns - 1, axis=1)[:, :ns]
            ind = np.zeros((c, n), dtype=np.float32)
            np.put_along_axis(ind, pick, 1.0, axis=1)
            Ap = (ind @ X).astype(np.float64)
        else:
            ### Pick the smaller side and subtract, for speed.
            nc = n - ns
            pick = np.argpartition(rng.random((c, n)), nc - 1, axis=1)[:, :nc]
            ind = np.zeros((c, n), dtype=np.float32)
            np.put_along_axis(ind, pick, 1.0, axis=1)
            Ap = T - (ind @ X).astype(np.float64)
        _, w = g_parts(np, Ap, T - Ap)
        null[done:done + c] = w
        done += c
    return {"obs": obs, "null": null}


def perm_p(obs, null):
    return (1 + int((null >= obs - 1e-9).sum())) / (1 + len(null))


def bh(ps):
    """Benjamini-Hochberg q values, in input order (None stays None)."""
    idx = [i for i, p in enumerate(ps) if p is not None]
    m = len(idx)
    q = [None] * len(ps)
    order = sorted(idx, key=lambda i: ps[i])
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * m / rank)
        q[i] = prev
    return q


def cell_contributions(np, A, B):
    """Each context's share of G(within): 2 * sum over the column's two cells
    of O ln(O/E), with E from that class's row and column totals.  Never
    negative (log-sum inequality), and sums to G(within)."""
    A = np.asarray(A, dtype=np.float64).reshape(-1, PER)
    B = np.asarray(B, dtype=np.float64).reshape(-1, PER)
    rA = A.sum(1, keepdims=True)
    rB = B.sum(1, keepdims=True)
    col = A + B
    N = rA + rB
    with np.errstate(divide="ignore", invalid="ignore"):
        EA = np.where(N > 0, rA * col / N, 0)
        EB = np.where(N > 0, rB * col / N, 0)
        tA = np.where(A > 0, A * np.log(A / EA), 0)
        tB = np.where(B > 0, B * np.log(B / EB), 0)
    return (2 * (tA + tB)).reshape(-1)


def within_shares(np, A):
    A = np.asarray(A, dtype=np.float64).reshape(-1, PER)
    r = A.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(r > 0, A / r, 0.0).reshape(-1)


# --------------------------------------------------------------------------
# tree
# --------------------------------------------------------------------------

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    parts = []
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def read_spectra(path):
    """Ordered split IDs, Total_Mutations per split, and whether the file has
    folded columns (a --fold search writes them, e.g. CT|GA)."""
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            ci = hdr.index("Node_ID")
        except ValueError:
            raise SystemExit("error: --spectra has no Node_ID column")
        ti = hdr.index("Total_Mutations") if "Total_Mutations" in hdr else None
        folded = any("|" in h for h in hdr)
        ids, totals = [], {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > ci and f[ci]:
                ids.append(f[ci])
                if ti is not None:
                    totals[f[ci]] = float(f[ti])
    return ids, totals, folded


def branch_matrix(np, tree, contexts, group, cap, internal_only, n_ctx):
    """(nodes, n_ctx) called-context counts per branch, plus per-group totals of
    parsed SNVs (called or not) for the check against spectra.tsv.  Branches
    over the cap, outside every group or excluded by --branches are zero."""
    n = len(tree.node_ids)
    X = np.zeros((n, n_ctx), dtype=np.int32)
    snvs = Counter()
    for i in range(n):
        g = group[i]
        if g < 0:
            continue
        if cap is not None and tree.raw_len[i] > cap:
            continue
        snvs[g] += len(tree.types[i])
        if internal_only and not tree.children[i]:
            continue
        ctx = contexts[i]
        if len(ctx):
            c = np.frombuffer(ctx, dtype=np.int16)
            c = c[c >= 0]
            if len(c):
                X[i] = np.bincount(c, minlength=n_ctx)
    return X, snvs


def top_sites(tree, contexts, members, cap, internal_only, n_ctx):
    """Per context, (share, position) of the most frequent site among a
    group's counted mutations in it, or (0.0, None) if it has none.  Counts
    the same branches as branch_matrix, so the shares are of the counts the
    tests use."""
    sites = [Counter() for _ in range(n_ctx)]
    for i in members:
        if cap is not None and tree.raw_len[i] > cap:
            continue
        if internal_only and not tree.children[i]:
            continue
        for p, c in zip(tree.positions[i], contexts[i]):
            if c >= 0:
                sites[c][p] += 1
    out = []
    for cnt in sites:
        tot = sum(cnt.values())
        if not tot:
            out.append((0.0, None))
            continue
        pos, k = cnt.most_common(1)[0]
        out.append((k / tot, pos))
    return out


def comparators(tree, group, roots, args):
    """split index -> comparator index (or None for the top split)."""
    index = {nid: i for i, nid in enumerate(tree.node_ids)}
    out = {}
    if args.against == "clade":
        for r in roots:
            p = tree.parent[r]
            out[r] = group[p] if p >= 0 and group[p] >= 0 else None
        return out
    from split_confidence import parse_splits_log
    with open(args.splits_log) as fh:
        _chi, parent, _order = parse_splits_log(fh)
    for r in roots:
        par = parent.get(tree.node_ids[r])
        out[r] = index.get(par) if par and index.get(par) in roots else None
    return out


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def fmt_p(p):
    if p is None:
        return "-"
    return f"{p:.2g}" if p >= 0.001 else f"{p:.1e}"


def run(np, tree, contexts, stats, group, roots, comp, totals, args, out,
        internal_only=False):
    """One run: every split on all branches, or on internal branches only.
    Starts its own random stream from --seed, so each run reproduces alone."""
    cap = args.max_branch_length
    sch = tn.scheme(args.fold)
    ncls = sch.ncls
    X, snvs = branch_matrix(np, tree, contexts, group, cap, internal_only,
                            sch.n)
    members = {}
    for i in range(len(tree.node_ids)):
        if group[i] >= 0:
            members.setdefault(group[i], []).append(i)
    mismatched = []
    for r in (roots if not internal_only else ()):
        nid = tree.node_ids[r]
        if nid in totals and abs(totals[nid] - snvs[r]) > 0.5:
            mismatched.append((nid, snvs[r], totals[nid]))
    if mismatched:
        print(f"  WARNING: {len(mismatched)} split(s) have a residual total "
              f"that differs from spectra.tsv, e.g. {mismatched[0][0]}: "
              f"{mismatched[0][1]} in the tree, {mismatched[0][2]:.0f} in "
              f"spectra.  Is this the tree and --max_branch_length the "
              f"search used?", file=sys.stderr)

    rng = np.random.default_rng(args.seed)
    rows = []
    ordered = sorted(roots, key=lambda r: tree.node_ids[r])
    for k, r in enumerate(ordered):
        c = comp.get(r)
        if c is None:
            continue
        ms = [i for i in members[r] if X[i].any()]
        mc = [i for i in members[c] if X[i].any()]
        Xs, Xc = X[ms], X[mc]
        A = Xs.sum(0, dtype=np.int64)
        B = Xc.sum(0, dtype=np.int64)
        row = {"node": tree.node_ids[r], "comparator": tree.node_ids[c],
               "split_mutations": int(A.sum()),
               "comparator_mutations": int(B.sum()),
               "split_branches": len(ms), "comparator_branches": len(mc),
               "A": A, "B": B}
        between, within = g_parts(np, A, B)
        row["g_between"] = float(between)
        row["g_within_cls"] = within
        row["g_within"] = float(within.sum())
        g_all = row["g_between"] + row["g_within"]
        row["within_share"] = row["g_within"] / g_all if g_all > 0 else 0.0
        if len(ms) < 2 or len(mc) < 2:
            row.update(p=None, null_mean=None, excess=None,
                       p_cls=[None] * ncls, excess_cls=[None] * ncls)
        else:
            res = shuffle_test(np, Xs, Xc, args.permutations, rng)
            tot_null = res["null"].sum(1)
            row["null_mean"] = float(tot_null.mean())
            row["excess"] = (row["g_within"] / row["null_mean"]
                             if row["null_mean"] > 0 else None)
            row["p"] = perm_p(row["g_within"], tot_null)
            row["p_cls"] = [perm_p(within[j], res["null"][:, j])
                            for j in range(ncls)]
            nm = res["null"].mean(0)
            row["excess_cls"] = [within[j] / nm[j] if nm[j] > 0 else None
                                 for j in range(ncls)]
        contrib = cell_contributions(np, A, B)
        row["contrib"] = contrib
        sa, sb = within_shares(np, A), within_shares(np, B)
        row["share_s"], row["share_c"] = sa, sb
        row["sites"] = top_sites(tree, contexts, members[r], cap,
                                 internal_only, sch.n)
        tops = []
        for ci in np.argsort(-contrib)[:args.top]:
            if contrib[ci] <= 0:
                break
            ts, pos = row["sites"][int(ci)]
            tops.append({"ctx": int(ci), "g": float(contrib[ci]),
                         "s": float(sa[ci]), "c": float(sb[ci]),
                         "site": ts, "pos": pos})
        row["tops"] = tops
        rA = A.reshape(ncls, PER).sum(1)
        rB = B.reshape(ncls, PER).sum(1)
        row["few"] = [bool(rA[j] < args.min_class or rB[j] < args.min_class)
                      for j in range(ncls)]
        rows.append(row)
        print(f"  {'internal' if internal_only else 'all'} "
              f"[{k + 1}/{len(ordered)}] {row['node']} vs "
              f"{row['comparator']}: G(within) {row['g_within']:.1f}, "
              f"p {fmt_p(row['p'])}", file=sys.stderr)

    qs = bh([r["p"] for r in rows])
    for r, q in zip(rows, qs):
        r["q"] = q
    ctrl = run_controls(np, tree, X, group, roots, members, rows, args, rng)
    if ctrl is None:
        for r in rows:
            r.setdefault("control_p", None)
            r.setdefault("control_n", None)
            r.setdefault("control_median", None)
    add_control_q(rows, args)
    report(rows, stats, args, out, ctrl, internal_only)
    return rows


def add_control_q(rows, args):
    """ctrl_q, the BH value of ctrl_p within the run, and for each split that
    passes the shuffle, whether its control count lets it reach
    --control_alpha at all: with ctrl_p at its floor, 1/(controls + 1),
    and every other split's as observed."""
    ps = [r.get("control_p") for r in rows]
    for r, q in zip(rows, bh(ps)):
        r["control_q"] = q
    for k, r in enumerate(rows):
        r["control_floor"] = None
        r["unreachable"] = False
        if r.get("control_n") is None:
            continue
        floor = 1.0 / (r["control_n"] + 1)
        r["control_floor"] = floor
        best = list(ps)
        best[k] = floor
        r["unreachable"] = bh(best)[k] > args.control_alpha


def clade_sizes(np, tree, X, group):
    """Called mutations in each node's clade, stopping at nested splits:
    the size the clade would have if it were carved out of its group."""
    sub = X.sum(1, dtype=np.int64)
    for i in range(len(tree.node_ids) - 1, -1, -1):
        p = tree.parent[i]
        if p >= 0 and group[i] >= 0 and group[p] == group[i]:
            sub[p] += sub[i]
    return sub


def clade_members(tree, group, v):
    out, stack = [], [v]
    while stack:
        i = stack.pop()
        out.append(i)
        stack.extend(c for c in tree.children[i] if group[c] == group[v])
    return out


def run_controls(np, tree, X, group, roots, members, rows, args, rng):
    """Size-matched random clades, tested exactly as the splits are.

    The branch shuffle asks whether a split differs from its comparator.  If
    context varies between sublineages everywhere on the tree, any clade
    would, so each split is also set against clades that are not splits:
    a node, its clade within the same residual group (nested splits carved
    out, as for a split) and the rest of that group as comparator.  For each
    split, --controls clades within --control_window fold of its size are
    drawn; control_p is the share of all drawn controls in that window whose
    excess reaches the split's.
    """
    if args.controls <= 0:
        return None
    tested = [r for r in rows if r["excess"] is not None]
    if not tested:
        return None
    sub = clade_sizes(np, tree, X, group)
    min_rest = min(r["split_mutations"] for r in tested)
    cand = [i for i in range(len(tree.node_ids))
            if group[i] >= 0 and i not in roots
            and sub[group[i]] - sub[i] >= min_rest]
    cand.sort(key=lambda i: sub[i])
    sizes = [int(sub[i]) for i in cand]
    import bisect
    drawn = {}
    tried = set()
    w = args.control_window
    for k, r in enumerate(tested):
        n = r["split_mutations"]
        lo = bisect.bisect_left(sizes, n / w)
        hi = bisect.bisect_right(sizes, n * w)
        pool = [cand[j] for j in range(lo, hi) if cand[j] not in tried]
        if not pool:
            continue
        pick = rng.choice(len(pool), size=min(args.controls, len(pool)),
                          replace=False)
        for j in pick:
            v = pool[int(j)]
            tried.add(v)
            clade = set(clade_members(tree, group, v))
            ms = [i for i in clade if X[i].any()]
            mc = [i for i in members[group[v]] if i not in clade
                  and X[i].any()]
            if len(ms) < 2 or len(mc) < 2:
                continue
            res = shuffle_test(np, X[ms], X[mc], args.control_permutations,
                               rng)
            tot = res["null"].sum(1)
            obs = float(res["obs"].sum())
            nm = float(tot.mean())
            drawn[v] = {"size": int(sub[v]), "excess": obs / nm if nm else 0,
                        "p": perm_p(obs, tot)}
        print(f"  controls [{k + 1}/{len(tested)}]: {len(drawn)} drawn",
              file=sys.stderr)
    ctl = list(drawn.values())
    for r in rows:
        r["control_n"] = r["control_median"] = r["control_p"] = None
        r["control_pool_left"] = None
        if r["excess"] is None:
            continue
        n = r["split_mutations"]
        lo = bisect.bisect_left(sizes, n / w)
        hi = bisect.bisect_right(sizes, n * w)
        ### Clades of this size on the tree not yet drawn: none left means
        ### more --controls cannot add any.
        r["control_pool_left"] = sum(1 for j in range(lo, hi)
                                     if cand[j] not in tried)
        win = [c["excess"] for c in ctl if n / w <= c["size"] <= n * w]
        if win:
            r["control_n"] = len(win)
            r["control_median"] = float(np.median(win))
            r["control_p"] = ((1 + sum(1 for e in win if e >= r["excess"]))
                              / (1 + len(win)))
    return ctl


def passes(r, args):
    """Both nulls: the branch shuffle (q) and the random clades (ctrl_q)."""
    return r.get("q") is not None and r["q"] < args.max_q \
        and r.get("control_q") is not None \
        and r["control_q"] <= args.control_alpha


def report(rows, stats, args, out, ctrl=None, internal_only=False):
    sch = tn.scheme(args.fold)
    L, ncls = sch.class_labels, sch.ncls
    print("=" * 72, file=out)
    print("TRINUCLEOTIDE CONTEXT WITHIN CLASSES, PER SPLIT: "
          + ("INTERNAL BRANCHES ONLY" if internal_only else "ALL BRANCHES"),
          file=out)
    print("=" * 72, file=out)
    print(f"  context classes: {sch.name}", file=out)
    tot = sum(stats.values()) or 1
    print(f"  contexts called: {stats['called']:,} of {tot:,} SNVs "
          f"({100 * stats['called'] / tot:.2f}%); unknown "
          f"{stats['unknown']:,}, doublet {stats['doublet']:,}", file=out)
    if stats["ref_mismatch"]:
        frac = stats["ref_mismatch"] / tot
        print(f"  own reference base disagrees with the parent genome: "
              f"{stats['ref_mismatch']:,} ({100 * frac:.2f}%)", file=out)
        if frac > 0.01:
            print("  WARNING: over 1%: --reference is probably not the "
                  "tree's root genome, and every context is suspect.",
                  file=out)
    print(f"  comparator: {'clade parent' if args.against == 'clade' else 'search subtree'}"
          f"; branches: {'internal only' if internal_only else 'all'}; "
          f"{args.permutations} branch shuffles "
          f"per split", file=out)
    print("  G(within) is the context shift inside classes; within% is its "
          "share of the", file=out)
    print(f"  2 x {sch.n} G; excess is G(within) over its shuffled mean (about 1 "
          "under the null).", file=out)
    print("  q is Benjamini-Hochberg over the splits.  Top class is the "
          "class with the", file=out)
    print("  largest G(within); '(few)' marks a class with fewer than "
          f"{args.min_class} mutations in a row.", file=out)
    tested = [r for r in rows if r["p"] is not None]
    hits = [r for r in tested if r["q"] is not None and r["q"] < args.max_q]
    beyond = [r for r in hits if passes(r, args)]
    print(f"\n  {len(tested)} splits tested; {len(hits)} with q < "
          f"{args.max_q:g}.", file=out)
    if ctrl:
        cp = sum(1 for c in ctrl if c["p"] <= 0.05)
        ex = sorted(c["excess"] for c in ctrl)
        sx = sorted(r["excess"] for r in tested if r["excess"] is not None)
        print(f"  {len(ctrl)} random non-split clades, size-matched to the "
              f"splits: {cp} ({100 * cp / len(ctrl):.0f}%) with p <= 0.05; "
              f"median excess {ex[len(ex) // 2]:.2f}, against "
              f"{sx[len(sx) // 2]:.2f} for the splits.", file=out)
        print("  ctrl_p is the share of controls within "
              f"{args.control_window:g}-fold of a split's size whose excess "
              "reaches the split's.", file=out)
        print("  If controls are often significant, context varies between "
              "clades everywhere,", file=out)
        print("  and only ctrl_p says a split differs more than any clade "
              "would.  ctrl_q is", file=out)
        print("  ctrl_p corrected (Benjamini-Hochberg) for the splits tested.",
              file=out)
        n_raw = sum(1 for r in hits if r.get("control_p") is not None
                    and r["control_p"] <= args.control_alpha)
        print(f"  {len(beyond)} splits have q < {args.max_q:g} and ctrl_q <= "
              f"{args.control_alpha:g} ({n_raw} with uncorrected ctrl_p <= "
              f"{args.control_alpha:g}).", file=out)
        stuck = [r for r in hits if r.get("unreachable")]
        if stuck:
            print(f"  {len(stuck)} splits with q < {args.max_q:g} have too few "
                  f"controls in their size window to reach ctrl_q <= "
                  f"{args.control_alpha:g}", file=out)
            print("  however strong their shift; ctrl_p can't fall below "
                  "1/(controls + 1).", file=out)
            more = [r for r in stuck if r.get("control_pool_left")]
            used = [r for r in stuck if not r.get("control_pool_left")]

            def lst(rs):
                return ", ".join(f"{r['node']} ({r['control_n']}, floor "
                                 f"{r['control_floor']:.2g})" for r in rs)
            if used:
                print(f"    every clade of their size on the tree is already "
                      f"drawn (only a wider --control_window can add more): "
                      f"{lst(used)}", file=out)
            if more:
                print(f"    more clades of their size exist (raise "
                      f"--controls): {lst(more)}", file=out)

    head = (f"  {'node':11s} {'vs':11s} {'split_n':>8s} {'comp_n':>8s} "
            f"{'G_betw':>8s} {'G_within':>8s} {'within%':>7s} {'excess':>6s} "
            f"{'p':>7s} {'q':>7s} {'ctrl_p':>6s} {'ctrl_q':>6s}  top class")
    print("\n" + head, file=out)
    for r in sorted(rows, key=lambda r: (r["p"] is None, r["p"] or 1,
                                         -r["g_within"])):
        j = int(max(range(ncls), key=lambda j: r["g_within_cls"][j]))
        exc = f"{r['excess']:.2f}" if r["excess"] is not None else "-"
        topc = (f"{L[j]} p {fmt_p(r['p_cls'][j])}"
                f"{' (few)' if r['few'][j] else ''}")
        print(f"  {r['node']:11s} {r['comparator']:11s} "
              f"{r['split_mutations']:8,d} {r['comparator_mutations']:8,d} "
              f"{r['g_between']:8.1f} {r['g_within']:8.1f} "
              f"{100 * r['within_share']:6.1f}% {exc:>6s} "
              f"{fmt_p(r['p']):>7s} {fmt_p(r['q']):>7s} "
              f"{fmt_p(r.get('control_p')):>6s} "
              f"{fmt_p(r.get('control_q')):>6s}  {topc}", file=out)

    if ctrl:
        show, what = beyond, (f"for splits with q < {args.max_q:g} and "
                              f"ctrl_q <= {args.control_alpha:g}")
    else:
        show, what = hits, f"for splits with q < {args.max_q:g}"
    if not show:
        show = sorted(tested, key=lambda r: (r["p"], -(r["excess"] or 0)))[:5]
        what = "for the five lowest p (none pass)"
    print(f"\n  Contexts carrying G(within), {what}:", file=out)
    print("  share = the context's share of its class, split vs comparator; "
          "site = share of", file=out)
    print("  the split's count in that context from its single most "
          "frequent site.", file=out)
    for r in show:
        cls = ", ".join(f"{L[j]} {fmt_p(r['p_cls'][j])}"
                        f"{'*' if r['few'][j] else ''}" for j in range(ncls))
        print(f"\n  {r['node']} (q {fmt_p(r['q'])}, ctrl_p "
              f"{fmt_p(r.get('control_p'))}, ctrl_q "
              f"{fmt_p(r.get('control_q'))}); per-class p: {cls}", file=out)
        for t in r["tops"]:
            arrow = "+" if t["s"] > t["c"] else "-"
            flag = "  <- one site" if t["site"] >= args.one_site else ""
            print(f"    {sch.labels[t['ctx']]:9s} {arrow} "
                  f"G {t['g']:7.1f}  share {t['s']:.3f} vs {t['c']:.3f}  "
                  f"site {t['site']:.2f} at {t['pos']}{flag}", file=out)
    if any(any(r["few"]) for r in show):
        print("\n  * fewer than --min_class mutations in that class in a row.",
              file=out)


def read_composition(path):
    with open(path) as fh:
        return {r["node"]: r for r in csv.DictReader(fh, delimiter="\t")}


def clade_label(r):
    """split_composition.py's top clades as 'name count, ...', plus flags."""
    parts = []
    for item in r.get("composition", "").split(","):
        if ":" in item:
            name, n = item.rsplit(":", 1)
            parts.append(f"{name} {n}")
    label = ", ".join(parts)
    flags = r.get("flags", "")
    return f"{label} ({flags})" if flags else label


def summarize(all_rows, int_rows, args, sch, comp=None):
    """The splits passing both nulls on both runs, as summary rows.

    Returns (rows, n_both, only_all, only_internal, one_site): rows sorted by
    all-branch excess, after --summary_top, if set, keeps the lowest
    all-branch ctrl_p;
    one_site lists (node, context, share, site) for summarised contexts whose
    top site carries --one_site or more.
    """
    by_int = {r["node"]: r for r in int_rows}
    ok_all = {r["node"] for r in all_rows if passes(r, args)}
    ok_int = {r["node"] for r in int_rows if passes(r, args)}
    both = ok_all & ok_int
    chosen = sorted((r for r in all_rows if r["node"] in both),
                    key=lambda r: r["control_p"])
    if args.summary_top:
        chosen = chosen[:args.summary_top]
    chosen.sort(key=lambda r: -r["excess"])
    out, one_site = [], []
    for a in chosen:
        i = by_int[a["node"]]
        order = sorted(range(sch.n), key=lambda c: -a["contrib"][c])
        listed = [c for c in order[:args.summary_contexts]
                  if a["contrib"][c] > 0]
        sites = [a["sites"][c] for c in listed]
        one_site += [(a["node"], sch.labels[c], sh, pos)
                     for c, (sh, pos) in zip(listed, sites)
                     if sh >= args.one_site]
        node = a["node"]
        out.append({
            "split": node,
            "comparator": a["comparator"],
            "clade": clade_label(comp[node]) if comp and node in comp else "",
            "excess_all": f"{a['excess']:.2f}",
            "excess_internal": f"{i['excess']:.2f}",
            "ctrl_p_all": f"{a['control_p']:.2g}",
            "ctrl_p_internal": f"{i['control_p']:.2g}",
            "ctrl_q_all": f"{a['control_q']:.2g}",
            "ctrl_q_internal": f"{i['control_q']:.2g}",
            "q_all": f"{a['q']:.2g}",
            "q_internal": f"{i['q']:.2g}",
            "contexts": "; ".join(
                f"{sch.labels[c]} {a['share_s'][c]:.3f} vs "
                f"{a['share_c'][c]:.3f}" for c in listed),
            "top_site_shares": "; ".join(f"{sh:.3f}" for sh, _p in sites),
            "max_site_share": f"{max(sh for sh, _p in sites):.3f}"
            if sites else "",
        })
    return (out, len(both), sorted(ok_all - ok_int), sorted(ok_int - ok_all),
            one_site)


SUMMARY_COLUMNS = ["split", "comparator", "clade", "excess_all",
                   "excess_internal", "ctrl_p_all", "ctrl_p_internal",
                   "ctrl_q_all", "ctrl_q_internal", "q_all", "q_internal", "contexts", "top_site_shares",
                   "max_site_share"]


def report_summary(summary, args, out):
    rows, n_both, only_all, only_int, one_site = summary
    print("=" * 72, file=out)
    print("SUMMARY: SPLITS THAT PASS BOTH NULLS ON BOTH RUNS", file=out)
    print("=" * 72, file=out)
    print(f"  {n_both} splits have q < {args.max_q:g} and ctrl_q <= "
          f"{args.control_alpha:g} on all branches AND on internal branches "
          f"only;", file=out)
    print(f"  {len(rows)} listed"
          f"{f' (the {args.summary_top} with the lowest all-branch ctrl_p)' if args.summary_top else ''}"
          f", by all-branch excess.", file=out)
    if only_all:
        print(f"  pass on all branches only: {', '.join(only_all)}", file=out)
    if only_int:
        print(f"  pass on internal branches only: {', '.join(only_int)}",
              file=out)
    print("  contexts: the class shares, split vs comparator, of the contexts "
          "carrying", file=out)
    print("  most of the all-branch G(within); sites: each one's top-site "
          "share.", file=out)
    print(file=out)
    if rows:
        cols = ["split", "excess_all", "excess_internal", "ctrl_p_all",
                "ctrl_p_internal", "clade", "contexts", "top_site_shares"]
        head = ["split", "excess all/int", "ctrl_q all/int", "clade"]
        print("  " + "\t".join(head), file=out)
        for r in rows:
            print(f"  {r['split']}\t{r['excess_all']} / "
                  f"{r['excess_internal']}\t{r['ctrl_q_all']} / "
                  f"{r['ctrl_q_internal']}\t{r['clade'] or '-'}", file=out)
            print(f"      {r['contexts']}  (sites {r['top_site_shares']})",
                  file=out)
        print(file=out)
    if one_site:
        print(f"  Contexts above where one site carries >= {args.one_site:g} "
              f"of the split's count:", file=out)
        for n, c, x, pos in one_site:
            print(f"    {n} {c}: {x:.2f} at {pos}", file=out)
    elif rows:
        big = max(float(r["max_site_share"]) for r in rows
                  if r["max_site_share"])
        print(f"  No context above has one site at >= {args.one_site:g}; the "
              f"largest top-site share is {big:.3f}.", file=out)
    print(file=out)


def write_summary(rows, path):
    with open(path, "w") as fh:
        fh.write("\t".join(SUMMARY_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in SUMMARY_COLUMNS) + "\n")


def write_tsv(rows, path, sch):
    L = sch.class_labels
    cols = ["node", "comparator", "split_mutations", "comparator_mutations",
            "split_branches", "comparator_branches", "g_between", "g_within",
            "within_share", "null_mean", "excess", "p", "q",
            "control_n", "control_median", "control_p", "control_q"]
    cols += [f"g_within_{c}" for c in L] + [f"p_{c}" for c in L]
    cols += ["top_contexts"]

    def f(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = [r.get(c) for c in cols[:17]]
            vals += [float(x) for x in r["g_within_cls"]] + list(r["p_cls"])
            vals.append(",".join(
                f"{sch.labels[t['ctx']]}{'+' if t['s'] > t['c'] else '-'}"
                f":{t['g']:.1f}" for t in r["tops"]))
            fh.write("\t".join(f(v) for v in vals) + "\n")


def write_contexts(rows, path, sch):
    with open(path, "w") as fh:
        fh.write("node\tcomparator\tclass\tcontext\tsplit_count\t"
                 "comparator_count\tsplit_class_share\t"
                 "comparator_class_share\tg_contribution\t"
                 "split_top_site\tsplit_top_site_share\n")
        for r in rows:
            for ci in range(sch.n):
                share, pos = r["sites"][ci]
                fh.write(f"{r['node']}\t{r['comparator']}\t"
                         f"{sch.class_labels[ci // PER]}\t"
                         f"{sch.labels[ci]}\t{int(r['A'][ci])}\t"
                         f"{int(r['B'][ci])}\t{r['share_s'][ci]:.6f}\t"
                         f"{r['share_c'][ci]:.6f}\t{r['contrib'][ci]:.4f}\t"
                         f"{'' if pos is None else pos}\t{share:.4f}\n")


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def selftest():
    import numpy as np
    rng = np.random.default_rng(7)

    # ---- the decomposition is exact --------------------------------
    for rep in range(40):
        n = 96 if rep % 2 else 192
        A = rng.integers(0, 50, n)
        B = rng.integers(0, 50, n)
        A[rng.integers(0, n, 10)] = 0
        b, w = g_parts(np, A, B)
        T = np.vstack([A, B]).astype(float)
        E = T.sum(1, keepdims=True) * T.sum(0, keepdims=True) / T.sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            g96 = 2 * np.where(T > 0, T * np.log(T / E), 0).sum()
        assert abs(b + w.sum() - g96) < 1e-6 * max(1, g96), (b, w.sum(), g96)
        assert abs(cell_contributions(np, A, B).sum() - w.sum()) < 1e-6
        assert (cell_contributions(np, A, B) >= -1e-9).all()
    print("ok: G(between) + G(within) equals the 2 x 96 and 2 x 192 G, and "
          "the per-context contributions sum to G(within)")

    ### A pure class shift has no within-class part, a pure context shift
    ### no between-class part.
    base = rng.integers(20, 60, 96).astype(float)
    cls_shift = base.copy()
    cls_shift[2 * PER:3 * PER] *= 3
    b, w = g_parts(np, base, cls_shift)
    assert w.sum() < 1e-9 and b > 10, (b, w)
    ctx_shift = base.copy()
    ctx_shift[2 * PER] += 15          # base is at least 20 everywhere
    ctx_shift[2 * PER + 1] -= 15
    b, w = g_parts(np, base * 10, ctx_shift * 10)
    assert b < 1e-6 and w[2] > 1 and w[[0, 1, 3, 4, 5]].sum() < 1e-9, (b, w)
    print("ok: a proportional class shift gives G(within) 0; a shift between "
          "two contexts of one class gives G(between) 0")

    # ---- calibration of the branch shuffle --------------------------
    probs = rng.dirichlet(np.ones(96) * 5)

    def branches(n, p, mean=15):
        sizes = rng.poisson(mean, n) + 1
        return np.array([rng.multinomial(s, p) for s in sizes])

    def classshift(p, k, f):
        q = p.copy()
        q[k * PER:(k + 1) * PER] *= f
        return q / q.sum()

    reps = 150
    null_p, prop_p, ctx_p = [], [], []
    shifted = probs.copy()
    ### The class's commonest context, so power does not hang on the draw.
    shifted[2 * PER + int(np.argmax(probs[2 * PER:3 * PER]))] *= 3.0
    shifted /= shifted.sum()
    for _ in range(reps):
        Xc = branches(300, probs)
        null_p.append(perm_p(*_obs_null(np, branches(120, probs), Xc, rng)))
        prop_p.append(perm_p(*_obs_null(
            np, branches(120, classshift(probs, 2, 3.0)), Xc, rng)))
        ctx_p.append(perm_p(*_obs_null(np, branches(120, shifted), Xc, rng)))
    fp = np.mean(np.array(null_p) <= 0.05)
    fpp = np.mean(np.array(prop_p) <= 0.05)
    pw = np.mean(np.array(ctx_p) <= 0.05)
    print(f"   p <= 0.05 in {reps} trials: null {fp:.3f}, proportional 3x "
          f"class shift {fpp:.3f}, one context 3x {pw:.3f}")
    assert fp < 0.10 and fpp < 0.10, (fp, fpp)
    assert pw > 0.5, pw
    print("ok: the shuffle is calibrated under the null and under a "
          "proportional class shift, and detects a one-context shift")

    # ---- the whole pipeline on a constructed tree -------------------
    class N(object):
        def __init__(self, i, muts, kids=()):
            self.id = i
            self.mutations = list(muts)
            self.children = list(kids)

    ref = "AACGTTACGTAACGTTAACC"
    ### a is a split whose leaves all have C>T at ACG (pos 3); root's other
    ### leaves have C>T at ACC (pos 19, A at 18... T at 20?).
    leaves_a = [N(f"a{i}", ["C3T"]) for i in range(6)]
    leaves_r = [N(f"r{i}", ["C19T"]) for i in range(6)]
    a = N("a", [], leaves_a)
    root = N("root", [], [a] + leaves_r)
    tree = extract(root)
    contexts, stats = tn.call_contexts(tree.parent, tree.children,
                                       tree.positions, tree.types, ref)
    group, roots, missing = residual_groups(tree, ["root", "a"])

    class Args(object):
        pass
    args = Args()
    args.against = "clade"
    args.max_branch_length = None
    args.permutations = 200
    args.seed = 1
    args.top = 2
    args.min_class = 1
    args.controls = 0
    args.control_alpha = 0.05
    args.max_q = 0.05
    args.one_site = 0.2
    args.summary_top = 0
    args.summary_contexts = 2
    args.fold = True
    comp = comparators(tree, group, roots, args)
    idx = {nid: i for i, nid in enumerate(tree.node_ids)}
    assert comp[idx["a"]] == idx["root"] and comp[idx["root"]] is None
    import io
    rows = run(np, tree, contexts, stats, group, roots, comp, {}, args,
               io.StringIO())
    assert len(rows) == 1
    r = rows[0]
    ### All C>T in both rows, so the class proportions are identical...
    assert r["g_between"] < 1e-9, r["g_between"]
    ### ...and the whole difference is context.
    assert r["g_within"] > 5 and r["within_share"] > 0.99
    want = {tn.context_index("A", "C", "G", "T"),
            tn.context_index("A", "C", "C", "T")}
    assert {t["ctx"] for t in r["tops"]} == want, r["tops"]
    assert r["p"] < 0.05, r["p"]
    ### Every split mutation is the same site, so its top site holds all.
    acg = tn.context_index("A", "C", "G", "T")
    assert r["sites"][acg] == (1.0, 3), r["sites"][acg]
    assert r["sites"][tn.context_index("A", "C", "C", "T")] == (0.0, None)
    print("ok: on a constructed tree the split's context shift is all "
          "within-class, names the two contexts and is significant, and "
          "the top site's share is counted per context")

    ### The same tree with 192 forward-strand classes finds the same two
    ### contexts, as forward-strand C>T.
    s192 = tn.scheme(False)
    c192, _st = s192.call(tree.parent, tree.children, tree.positions,
                          tree.types, ref)
    args.fold = False
    r = run(np, tree, c192, stats, group, roots, comp, {}, args,
            io.StringIO())[0]
    assert r["g_between"] < 1e-9 and r["within_share"] > 0.99
    assert {t["ctx"] for t in r["tops"]} == {
        tn.stranded_context_index("A", "C", "G", "T"),
        tn.stranded_context_index("A", "C", "C", "T")}, r["tops"]
    print("ok: without --fold the same tree gives the same result in the 192 "
          "forward-strand classes")

    ### A strand-only difference: the split's leaves carry C>T at ACG and
    ### the rest G>A at CGT, the same folded context on the other strand.
    ### Folded, nothing differs; unfolded, the class (type) differs.
    ref2 = "AACGTTAAAACGTTAA"
    la = [N(f"a{i}", ["C3T"]) for i in range(6)]
    lr = [N(f"r{i}", ["G12A"]) for i in range(6)]
    tree2 = extract(N("root", [], [N("a", [], la)] + lr))
    g2, roots2, _m = residual_groups(tree2, ["root", "a"])
    comp2 = comparators(tree2, g2, roots2, args)
    for fold in (True, False):
        args.fold = fold
        c2, st2 = tn.scheme(fold).call(tree2.parent, tree2.children,
                                        tree2.positions, tree2.types, ref2)
        r = run(np, tree2, c2, st2, g2, roots2, comp2, {}, args,
                io.StringIO())[0]
        if fold:
            assert r["g_between"] < 1e-9 and r["g_within"] < 1e-9, r
        else:
            assert r["g_between"] > 5 and r["g_within"] < 1e-9, r
    args.fold = True
    print("ok: a difference only in strand is invisible folded and a class "
          "difference in the 192 classes")

    ### Clade sizes stop at nested splits: root's clade is its own residual.
    X, _ = branch_matrix(np, tree, contexts, group, None, False, 96)
    sz = clade_sizes(np, tree, X, group)
    assert sz[idx["root"]] == 6 and sz[idx["a"]] == 6, sz
    assert sz[idx["a0"]] == 1
    assert sorted(clade_members(tree, group, idx["root"])) == sorted(
        [idx["root"]] + [idx[f"r{i}"] for i in range(6)])
    print("ok: control clades stop at nested splits, as a split's residual "
          "does")

    ### The internal-only run drops the leaves, leaving nothing to test.
    rows = run(np, tree, contexts, stats, group, roots, comp, {}, args,
               io.StringIO(), internal_only=True)
    assert rows[0]["split_mutations"] == 0
    print("ok: the internal-only run counts no terminal branch")

    _selftest_summary(np)
    print("all self-tests passed")
    return 0


def _selftest_summary(np):
    sch = tn.scheme(True)

    class A(object):
        pass
    args = A()
    args.max_q, args.control_alpha, args.one_site = 0.05, 0.05, 0.2
    args.summary_top, args.summary_contexts = 0, 2

    def row(node, excess, q, ctrl, contrib=None, sites=None):
        c = np.zeros(sch.n)
        for k, v in (contrib or {}).items():
            c[k] = v
        st = [(0.0, None)] * sch.n
        for k, v in (sites or {}).items():
            st[k] = v
        sh = np.full(sch.n, 0.1)
        return {"node": node, "comparator": "r", "excess": excess, "q": q,
                "control_p": ctrl, "control_q": ctrl, "contrib": c,
                "share_s": sh,
                "share_c": sh / 2, "sites": st}
    all_rows = [
        row("a", 5.0, 0.001, 0.003, {40: 9.0, 41: 3.0, 42: 0.0},
            {40: (0.01, 100), 41: (0.25, 17)}),
        row("b", 7.0, 0.001, 0.020, {1: 1.0}, {1: (0.0, 5)}),
        row("c", 3.0, 0.001, 0.040, {1: 1.0}, {1: (0.0, 5)}),
        row("d", 9.0, 0.001, 0.003),               # fails internal
        row("e", 2.0, 0.200, 0.003),               # fails q
        row("f", 2.0, 0.001, 0.300),               # fails ctrl_p
    ]
    int_rows = [row("a", 2.5, 0.001, 0.005), row("b", 3.0, 0.010, 0.040),
                row("c", 1.5, 0.020, 0.050), row("d", 1.0, 0.500, 0.600),
                row("e", 1.0, 0.001, 0.003), row("f", 1.0, 0.001, 0.003)]
    comp = {"a": {"composition": "L1:6,L2:4", "flags": "mixed"}}
    rows, nboth, only_all, only_int, one = summarize(all_rows, int_rows, args,
                                                     sch, comp)
    assert [r["split"] for r in rows] == ["b", "a", "c"], rows
    assert nboth == 3 and only_all == ["d"] and only_int == ["e", "f"]
    a = rows[1]
    assert a["contexts"] == (f"{sch.labels[40]} 0.100 vs 0.050; "
                             f"{sch.labels[41]} 0.100 vs 0.050"), a
    assert a["clade"] == "L1 6, L2 4 (mixed)" and rows[0]["clade"] == ""
    assert a["top_site_shares"] == "0.010; 0.250" and \
        a["max_site_share"] == "0.250"
    assert one == [("a", sch.labels[41], 0.25, 17)], one
    args.summary_top = 2
    rows = summarize(all_rows, int_rows, args, sch)[0]
    assert [r["split"] for r in rows] == ["b", "a"], rows
    ### ctrl_q and the reachability check.  Five splits; the first has
    ### only 9 controls, so its floor 0.1 can never pass at 0.05.
    rs = [{"control_p": p, "control_n": n} for p, n in
          ((0.1, 9), (0.002, 400), (0.004, 400), (0.5, 400), (None, None))]
    add_control_q(rs, args)
    want_q = [0.4 / 3, 0.008, 0.008, 0.5, None]   # 0.1 is 3rd of 4
    assert all((a["control_q"] is None and b is None)
               or abs(a["control_q"] - b) < 1e-12
               for a, b in zip(rs, want_q)), [r["control_q"] for r in rs]
    assert [r["unreachable"] for r in rs] == [True, False, False, False,
                                              False], rs
    print("ok: ctrl_q is BH over the run's ctrl_p, and a split whose control "
          "count can't reach --control_alpha is flagged")
    print("ok: the summary keeps splits passing both nulls on both runs, "
          "sorts by excess, lists contexts by G with their site shares, and "
          "flags one site; --summary_top keeps the lowest ctrl_p")


def _obs_null(np, Xs, Xc, rng):
    res = shuffle_test(np, Xs, Xc, 199, rng)
    return res["obs"].sum(), res["null"].sum(1)


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    for opt in ("input_tree", "spectra", "reference"):
        if not getattr(args, opt):
            print(f"error: --{opt} is required (or pass --selftest)",
                  file=sys.stderr)
            return 2
    if args.against == "search" and not args.splits_log:
        print("error: --against search needs --splits_log", file=sys.stderr)
        return 2
    try:
        import numpy as np
    except ImportError:
        print("error: numpy is required", file=sys.stderr)
        return 1
    try:
        import bte
    except ImportError:
        print("error: could not import bte; run this in the environment you "
              "use for spectrumSplits.py", file=sys.stderr)
        return 1

    split_ids, totals, folded = read_spectra(args.spectra)
    print(f"{len(split_ids)} splits in {args.spectra}", file=sys.stderr)
    if folded != args.fold:
        print(f"  WARNING: {args.spectra} looks like a search run "
              f"{'with' if folded else 'without'} --fold, but this run is "
              f"{'with' if args.fold else 'without'} it.  Match the search: "
              f"96 folded classes for a --fold search, 192 forward-strand "
              f"ones for an unfolded one.", file=sys.stderr)
    sch = tn.scheme(args.fold)
    ref = read_fasta(args.reference)
    masked = tn.read_mask(args.mask.split(",") if args.mask else [], len(ref))
    print(f"Reading {args.input_tree} ...", file=sys.stderr)
    mat = bte.MATree(args.input_tree)
    tree = extract(mat.root)
    del mat
    group, roots, missing = residual_groups(tree, split_ids)
    if missing:
        print(f"  WARNING: {len(missing)} split root(s) not in this tree, "
              f"e.g. {missing[0]}; node IDs are assigned when a tree is read, "
              f"so use the tree the search ran on.", file=sys.stderr)
    print(f"Calling contexts, {sch.name} ...", file=sys.stderr)
    contexts, stats = sch.call(tree.parent, tree.children, tree.positions,
                               tree.types, ref, masked)
    comp = comparators(tree, group, roots, args)
    comp_tab = read_composition(args.composition) if args.composition \
        else None
    out = sys.stdout
    runs = {}
    for internal_only in (False, True):
        runs[internal_only] = run(np, tree, contexts, stats, group, roots,
                                  comp, totals, args, out, internal_only)
        print(file=out)
    summary = summarize(runs[False], runs[True], args, sch, comp_tab)
    report_summary(summary, args, out)
    pre = args.prefix
    paths = []
    for internal_only, tag in ((False, ""), (True, "_internal")):
        paths.append(f"{pre}{tag}.tsv")
        write_tsv(runs[internal_only], paths[-1], sch)
        paths.append(f"{pre}{tag}_{sch.n}.tsv")
        write_contexts(runs[internal_only], paths[-1], sch)
    paths.append(f"{pre}_table.tsv")
    write_summary(summary[0], paths[-1])
    print(f"Wrote {', '.join(paths)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
