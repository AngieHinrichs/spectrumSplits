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

What it cannot tell you: whether a search scoring contexts would find splits
this one did not.  That needs the search to score contexts.

Usage:
    split_contexts.py --input_tree tree.pruned.masked.pb.gz \\
        --spectra spectra.tsv --fold --max_branch_length 250 \\
        --reference root_genome.fa --mask selection_loci.bed,alignment_loci.bed \\
        --output split_contexts.tsv --output_contexts split_contexts_96.tsv \\
        > split_contexts.txt

    split_contexts.py --selftest
"""

import argparse
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
    p.add_argument("--branches", choices=["all", "internal"], default="all",
                   help="internal: count only internal branches, to keep "
                        "per-genome error out (as in step 9)")
    p.add_argument("--permutations", type=int, default=1000,
                   help="Branch shuffles per split (default 1000)")
    p.add_argument("--controls", type=int, default=10,
                   help="Random non-split clades drawn per split, within "
                        "--control_window of its size, as a reference for "
                        "how much any clade differs (default 10; 0 skips)")
    p.add_argument("--control_window", type=float, default=2.0,
                   help="Fold range of size a control may differ from a "
                        "split by (default 2)")
    p.add_argument("--control_permutations", type=int, default=200,
                   help="Branch shuffles per control clade (default 200)")
    p.add_argument("--control_alpha", type=float, default=0.05,
                   help="ctrl_p at or below which a split's contexts are "
                        "listed in detail (default 0.05)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--top", type=int, default=3,
                   help="Contexts to list per split (default 3)")
    p.add_argument("--min_class", type=int, default=80,
                   help="A class with fewer called mutations than this in "
                        "either row is marked 'few' in the report "
                        "(default 80, five per context)")
    p.add_argument("--output", default=None, help="Per-split TSV")
    p.add_argument("--output_contexts", default=None,
                   help="Per-split, per-context TSV")
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


def run(np, tree, contexts, stats, group, roots, comp, totals, args, out):
    cap = args.max_branch_length
    sch = tn.scheme(args.fold)
    ncls = sch.ncls
    X, snvs = branch_matrix(np, tree, contexts, group, cap,
                            args.branches == "internal", sch.n)
    members = {}
    for i in range(len(tree.node_ids)):
        if group[i] >= 0:
            members.setdefault(group[i], []).append(i)
    mismatched = []
    for r in roots:
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
                                 args.branches == "internal", sch.n)
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
        print(f"  [{k + 1}/{len(ordered)}] {row['node']} vs "
              f"{row['comparator']}: G(within) {row['g_within']:.1f}, "
              f"p {fmt_p(row['p'])}", file=sys.stderr)

    qs = bh([r["p"] for r in rows])
    for r, q in zip(rows, qs):
        r["q"] = q
    ctrl = run_controls(np, tree, X, group, roots, members, rows, args, rng)
    report(rows, stats, args, out, ctrl)
    return rows


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
    w = args.control_window
    for k, r in enumerate(tested):
        n = r["split_mutations"]
        lo = bisect.bisect_left(sizes, n / w)
        hi = bisect.bisect_right(sizes, n * w)
        pool = [cand[j] for j in range(lo, hi) if cand[j] not in drawn]
        if not pool:
            continue
        pick = rng.choice(len(pool), size=min(args.controls, len(pool)),
                          replace=False)
        for j in pick:
            v = pool[int(j)]
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
        if r["excess"] is None:
            continue
        n = r["split_mutations"]
        win = [c["excess"] for c in ctl if n / w <= c["size"] <= n * w]
        if win:
            r["control_n"] = len(win)
            r["control_median"] = float(np.median(win))
            r["control_p"] = ((1 + sum(1 for e in win if e >= r["excess"]))
                              / (1 + len(win)))
    return ctl


def report(rows, stats, args, out, ctrl=None):
    sch = tn.scheme(args.fold)
    L, ncls = sch.class_labels, sch.ncls
    print("=" * 72, file=out)
    print("TRINUCLEOTIDE CONTEXT WITHIN CLASSES, PER SPLIT", file=out)
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
          f"; branches: {args.branches}; {args.permutations} branch shuffles "
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
    hits = [r for r in tested if r["q"] is not None and r["q"] < 0.05]
    beyond = [r for r in hits if r.get("control_p") is not None
              and r["control_p"] <= args.control_alpha]
    print(f"\n  {len(tested)} splits tested; {len(hits)} with q < 0.05.",
          file=out)
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
              "would.", file=out)
        print(f"  {len(beyond)} splits have q < 0.05 and ctrl_p <= "
              f"{args.control_alpha:g}.", file=out)

    head = (f"  {'node':11s} {'vs':11s} {'split_n':>8s} {'comp_n':>8s} "
            f"{'G_betw':>8s} {'G_within':>8s} {'within%':>7s} {'excess':>6s} "
            f"{'p':>7s} {'q':>7s} {'ctrl_p':>6s}  top class")
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
              f"{fmt_p(r.get('control_p')):>6s}  {topc}", file=out)

    if ctrl:
        show, what = beyond, (f"for splits with q < 0.05 and ctrl_p <= "
                              f"{args.control_alpha:g}")
    else:
        show, what = hits, "for splits with q < 0.05"
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
              f"{fmt_p(r.get('control_p'))}); per-class p: {cls}", file=out)
        for t in r["tops"]:
            arrow = "+" if t["s"] > t["c"] else "-"
            flag = "  <- one site" if t["site"] >= 0.2 else ""
            print(f"    {sch.labels[t['ctx']]:9s} {arrow} "
                  f"G {t['g']:7.1f}  share {t['s']:.3f} vs {t['c']:.3f}  "
                  f"site {t['site']:.2f} at {t['pos']}{flag}", file=out)
    if any(any(r["few"]) for r in show):
        print("\n  * fewer than --min_class mutations in that class in a row.",
              file=out)


def write_tsv(rows, path, sch):
    L = sch.class_labels
    cols = ["node", "comparator", "split_mutations", "comparator_mutations",
            "split_branches", "comparator_branches", "g_between", "g_within",
            "within_share", "null_mean", "excess", "p", "q",
            "control_n", "control_median", "control_p"]
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
            vals = [r.get(c) for c in cols[:16]]
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
    args.branches = "all"
    args.permutations = 200
    args.seed = 1
    args.top = 2
    args.min_class = 1
    args.controls = 0
    args.control_alpha = 0.05
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

    ### --branches internal drops the leaves, leaving nothing to test.
    args.branches = "internal"
    rows = run(np, tree, contexts, stats, group, roots, comp, {}, args,
               io.StringIO())
    assert rows[0]["split_mutations"] == 0
    print("ok: --branches internal counts no terminal branch")
    print("all self-tests passed")
    return 0


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
    rows = run(np, tree, contexts, stats, group, roots, comp, totals, args,
               sys.stdout)
    if args.output:
        write_tsv(rows, args.output, sch)
        print(f"Per-split table written to {args.output}", file=sys.stderr)
    if args.output_contexts:
        write_contexts(rows, args.output_contexts, sch)
        print(f"Per-context table written to {args.output_contexts}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
