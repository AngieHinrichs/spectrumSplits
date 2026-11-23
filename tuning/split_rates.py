#!/usr/bin/env python3
"""Is a split a mutator, a spectrum-shifter, or both?

spectrumSplits.py finds clades whose mutation SPECTRUM differs.  A spectrum
can change in two ways that point to different mechanisms: a clade can gain
mutations of some classes on top of everything else (a class-specific
mutator, such as a lost repair pathway), or it can trade some classes for
others at about the same overall rate.  The spectrum alone can't tell these
apart.  The RATE can, so this script measures each split's rate against its
relatives, and says which classes carry any difference.

There is no clock.  On a mutation-annotated tree a branch's length IS its
mutation count, so "mutations per unit time" has to come from a comparison in
which time is held fixed.  That is the relative-rate test: a split's clade and
its sister clades descend from one ancestor A, so if the genomes were sampled
at about the same time, every lineage from A to a tip has had the same time
to mutate.  A clade whose tips sit further from A has mutated faster.

For each split S (its whole clade, nested splits included):

  A          the nearest ancestor of S with at least --min_outgroup_tips tips
             outside S.  Usually S's parent; a step or two up when S's sisters
             are a few tips, so the outgroup isn't one genome.
  outgroup   the tips below A that are not in S.
  depth      a tip's distance from A, in mutations of each class.  mean_S and
             mean_out are the means over each side's tips.
  log2_rate  log2(mean_S / mean_out), on all classes.  0 is the same rate.

Tips share branches, so their depths are not independent and no analytic test
applies.  The null is every other clade on the tree, scored the same way:
each non-split node with A defined by the same rule, whose tip count is
within --control_window fold of the split's (or, when fewer than
--min_controls are, that many nearest in size).  ctrl_p is two-sided: the share
of those controls at least as far from the controls' median as the split.  q
is Benjamini-Hochberg over the splits.  A control is left out when one split
holds --control_overlap or more of its tips, since it would be scoring that
split again.  A control whose OUTGROUP is a fast split will look slow; the
two-sided test absorbs some of that, and the report gives the controls'
spread so it can be judged.

Classes.  The same is done for each mutation class (6 folded with --fold, 12
without, to match the search), with per-class controls.  "excess" is the
split's extra mutations per tip, mean_S - mean_out, and each class's share
of it says where any rate difference lies.  The summary's "excess" column
gives the three largest d_k = mean_S,k - mean_out,k as fractions of the
outgroup's mean depth, so they add up (over all classes) to the total rate
change.  "mix" summarises them:

  mix = 1 - |sum d_k| / sum |d_k|

0 when every class moves the same way (a change of rate, spectrum or not),
1 when gains and losses cancel (a spectrum shift at an unchanged rate).

The shape test asks whether a rate change is class-specific, as a repair
defect's would be, or spread over the classes as in random clades: each
class's log2 ratio minus the overall one, against the same in the
controls, two-sided, Bonferroni over the classes (shape_p).

Two runs, as in steps 9 and 11: all branches, then internal branches only.
Terminal branches carry per-genome sequencing error and quality differences
between studies, which inflate a clade's depth without any change in its
mutation rate.  A rate difference seen on all branches but not on internal
ones is more likely that than a mutator.

The reading.  The report opens with a READING section that applies the
decision rules and says what they found, in words: a conclusion, the
resolution random clades allow, what the dates did, and each split by
category with its reason:
  mutator, antimutator   q passes on both runs (and after dates, if any),
                         and the shape test passes: a class-specific change
  genome-wide faster/slower
                         the same, but spread over the classes alike
  sampling time          passes undated, not once dates are allowed for
  terminal only          passes on all branches, not on internal ones
  suggestive faster/slower
                         p but not q passes on both runs
  unchanged rate         anything else ("classes traded" when mix >= 0.5)
  untestable             too few controls for p to reach --max_q
A reason notes a comparison ancestor within --near_root nodes of the root.
The per-run tables and the summary table follow, under DETAILS.  In the
summary, "verdict" is the undated total-rate call on both runs, and
"classes" the classes whose own rate differs on both runs.

Dates (optional).  The test assumes S and its outgroup were sampled at
about the same time.  A clade sampled earlier (older isolates, ancient DNA)
has shorter tip depths and looks slower; one sampled later looks faster.
--dates reads collection dates from a TSV (--sample_column, --date_column)
and adds, on all branches only, using dated tips alone:

  clock          the tree's root-to-tip regression: mutations per year
                 against collection date over every dated tip, with r^2.
                 A slope at or below zero means the dates carry no clock
                 signal, and the report says so.
  adj_log2_rate  log2_rate with S's mean depth moved to the outgroup's mean
                 date along that clock:
                   mean_S - slope * (mean date_S - mean date_out)
                 Scored against controls adjusted the same way.  This is
                 the date-corrected answer to "faster or slower".
  slope_S, slope_out
                 root-to-tip slopes within S and within the outgroup:
                 mutations per year measured directly.  Their difference is
                 scored against the controls'.  It needs a spread of dates
                 on both sides, and is much noisier than adj_log2_rate.

Not on internal branches: a tip's date says how long its whole path ran,
and dropping the terminal branch drops the part of the path the sampling
date governs.  So read the date-adjusted verdict beside the undated
internal-branch run, which guards against terminal-branch error.  Dates may
be YYYY, YYYY-MM, YYYY-MM-DD (with or without a time after it), a decimal
year, or a range "A/B", which is used at its midpoint when its ends are at
most a year apart.  A partial date is placed at the middle of its interval;
--min_date_precision drops dates coarser than a year, month or day.  A
sample listed twice with different dates is dropped.

What it cannot tell you.  It measures the rate that survived selection:
purifying selection removes mutations from every clade alike, but a clade
under relaxed selection will look faster.

Inputs: the tree the search ran on, and its --output_spectrum TSV for the
split roots.  --max_branch_length should be the search's: branches over it
are left out, as the search left them out.  That shortens the depth of every
tip below such a branch, so the report says how many were dropped.

Outputs, named from --prefix (default split_rates):
  PREFIX.tsv           one row per split: both runs' total-rate numbers,
                       mix, verdict and classes, with --dates the
                       date-adjusted rate, slopes and date_verdict, and the
                       reading, shape_p, concentration and reason
  PREFIX_classes.tsv   one row per split, run and class
  PREFIX_controls.tsv  every control's tips and log2_rate, per run (and
                       adj_log2_rate and slope_diff with --dates)
The report goes to stdout.

Usage:
    split_rates.py --input_tree tree.pruned.masked.pb.gz \\
        --spectra spectra.tsv --fold --max_branch_length 250 \\
        --composition split_composition.tsv > split_rates.txt

    split_rates.py ... --dates metadata.tsv --sample_column biosample \\
        --date_column collection_date > split_rates.txt

    split_rates.py --selftest
"""

import argparse
import bisect
import csv
import datetime
import gzip
import io
import textwrap
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from permutation_null import MUT_TYPES, folded_index

PSEUDO = 0.5


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Test each split's mutation rate against its sister "
                    "clades, by class, against size-matched control clades.")
    p.add_argument("--input_tree", help="The tree the search ran on (protobuf)")
    p.add_argument("--spectra",
                   help="spectrumSplits.py --output_spectrum TSV; split "
                        "roots are read from it")
    p.add_argument("--fold", action="store_true",
                   help="Score the 6 strand-folded classes (match the search)")
    p.add_argument("--max_branch_length", type=int, default=None,
                   help="Leave out branches with more mutations than this, "
                        "as the search did (default none)")
    p.add_argument("--min_outgroup_tips", type=int, default=10,
                   help="Walk up from the split until the ancestor has at "
                        "least this many tips outside it (default 10)")
    p.add_argument("--control_window", type=float, default=2.0,
                   help="Controls are clades within this fold of the "
                        "split's tip count (default 2)")
    p.add_argument("--min_controls", type=int, default=40,
                   help="When fewer clades than this are within "
                        "--control_window of a split's size, use this many "
                        "nearest in size instead (default 40)")
    p.add_argument("--control_overlap", type=float, default=0.5,
                   help="Leave out a control when one split holds this "
                        "share of its tips or more (default 0.5)")
    p.add_argument("--max_q", type=float, default=0.05,
                   help="q at or below this passes (default 0.05)")
    p.add_argument("--dates", default=None,
                   help="Optional TSV of sample collection dates, with a "
                        "header line, gzipped if it ends in .gz; adds a date-adjusted rate and "
                        "root-to-tip clock slopes")
    p.add_argument("--sample_column", default=None,
                   help="--dates column holding the tree's sample names "
                        "(header name or 1-based number; default the first)")
    p.add_argument("--date_column", default=None,
                   help="--dates column holding the dates (header name or "
                        "1-based number): YYYY, YYYY-MM, YYYY-MM-DD or a "
                        "decimal year")
    p.add_argument("--min_date_precision", choices=["year", "month", "day"],
                   default="year",
                   help="Ignore dates less precise than this (default year)")
    p.add_argument("--min_dated_tips", type=int, default=10,
                   help="Score dates only where both sides have this many "
                        "dated tips (default 10)")
    p.add_argument("--near_root", type=int, default=5,
                   help="Note a split whose comparison ancestor is within "
                        "this many nodes of the root (default 5)")
    p.add_argument("--composition", default=None,
                   help="split_composition.py TSV, to label splits by clade")
    p.add_argument("--prefix", default="split_rates",
                   help="Prefix of the output TSVs (default split_rates)")
    p.add_argument("--selftest", action="store_true",
                   help="Run the built-in tests and exit")
    return p.parse_args(argv)


# --------------------------------------------------------------------------
# depths
# --------------------------------------------------------------------------

def branch_counts(np, tree, ncat, fold_map, cap, internal_only):
    """(nodes, ncat) class counts per branch.  Capped branches are zero, and
    with internal_only so are terminal ones."""
    n = len(tree.node_ids)
    L = np.zeros((n, ncat), dtype=np.float64)
    capped = 0
    for i in range(n):
        if internal_only and not tree.children[i]:
            continue
        if cap is not None and tree.raw_len[i] > cap:
            capped += 1
            continue
        for t in tree.types[i]:
            L[i, fold_map[t] if fold_map else t] += 1
    return L, capped


def tip_counts(np, tree):
    n = len(tree.node_ids)
    T = np.zeros(n, dtype=np.int64)
    for i in range(n - 1, -1, -1):
        if not tree.children[i]:
            T[i] += 1
        p = tree.parent[i]
        if p >= 0:
            T[p] += T[i]
    return T


def depth_sums(np, tree, L, T):
    """D[i]: summed distance from node i (its own branch excluded) to every
    tip below it, per class.  R[i]: distance from the root to node i, its own
    branch included.  Parents precede children, so one pass each way."""
    n = len(tree.node_ids)
    D = np.zeros_like(L)
    for i in range(n - 1, 0, -1):
        p = tree.parent[i]
        if p >= 0:
            D[p] += D[i] + T[i] * L[i]
    R = np.zeros_like(L)
    for i in range(1, n):
        p = tree.parent[i]
        if p >= 0:
            R[i] = R[p] + L[i]
    return D, R


def outgroup_ancestors(tree, T, min_out):
    """node -> its comparison ancestor A, or -1 when none has enough tips."""
    n = len(tree.node_ids)
    A = [-1] * n
    for v in range(n):
        a = tree.parent[v]
        while a >= 0 and T[a] - T[v] < min_out:
            a = tree.parent[a]
        A[v] = a
    return A


def rates(np, nodes, anc, T, D, R):
    """Per node: (mean_S, mean_out) per class, as two (len(nodes), ncat)
    arrays, and the outgroup's tip counts."""
    v = np.asarray(nodes, dtype=np.int64)
    a = np.asarray([anc[i] for i in nodes], dtype=np.int64)
    ts = T[v].astype(np.float64)
    to = (T[a] - T[v]).astype(np.float64)
    s_sum = D[v] + ts[:, None] * (R[v] - R[a])
    o_sum = D[a] - s_sum
    return s_sum / ts[:, None], o_sum / to[:, None], to.astype(np.int64)


def log2_rate(np, ms, mo, ts, to):
    """log2 of the mean depth ratio, with a half-mutation pseudocount on each
    side's sum so a class absent from one side stays finite."""
    return np.log2(((ms * ts[..., None] + PSEUDO) / ts[..., None])
                   / ((mo * to[..., None] + PSEUDO) / to[..., None]))


# --------------------------------------------------------------------------
# dates
# --------------------------------------------------------------------------

DATE_RE = re.compile(r"^(\d{4})(?:-(\d{1,2}|XX)(?:-(\d{1,2}|XX))?)?$",
                     re.IGNORECASE)
PRECISION = {"year": 0, "month": 1, "day": 2}


def parse_date(text):
    """(decimal year, precision) or None.  Precision 0 is a year, 1 a
    month, 2 a day; a partial date is placed at the middle of its interval.
    A decimal year (2010.25) counts as a day.  A time after the date
    (2010-03-01T12:00:00Z) is ignored.  A range (2001/2003, as NCBI
    BioSample writes them) is its midpoint, with precision -1, below a
    year, when its ends are more than a year apart."""
    t = text.strip()
    if "/" in t:
        ends = [parse_date(x) for x in t.split("/")]
        if len(ends) != 2 or None in ends:
            return None
        lo, hi = sorted(e[0] for e in ends)
        prec = min(e[1] for e in ends) if hi - lo <= 1 else -1
        return (lo + hi) / 2, prec
    t = re.sub(r"^(\d{4}-\d{2}-\d{2})[T ][\d:.]+Z?$", r"\1", t)
    m = DATE_RE.match(t)
    if m:
        y = int(m.group(1))
        mo, d = m.group(2), m.group(3)
        ylen = 366 if (y % 4 == 0 and y % 100 != 0) or y % 400 == 0 else 365
        if mo is None or mo.upper() == "XX":
            return y + 0.5, 0
        mo = int(mo)
        if not 1 <= mo <= 12:
            return None
        start = datetime.date(y, mo, 1).timetuple().tm_yday - 1
        nxt = datetime.date(y + (mo == 12), mo % 12 + 1, 1)
        days = (nxt - datetime.date(y, mo, 1)).days
        if d is None or d.upper() == "XX":
            return y + (start + days / 2) / ylen, 1
        d = int(d)
        if not 1 <= d <= days:
            return None
        return y + (start + d - 0.5) / ylen, 2
    try:
        f = float(t)
    except ValueError:
        return None
    if not math.isfinite(f):
        return None
    return f, 2


def column_index(header, spec, default=None):
    if spec is None:
        return default
    if spec in header:
        return header.index(spec)
    if spec.isdigit() and 1 <= int(spec) <= len(header):
        return int(spec) - 1
    raise SystemExit(f"error: no column {spec!r} in --dates; it has "
                     f"{', '.join(header)}")


def read_dates(lines, sample_col, date_col, min_precision):
    """sample -> decimal year, and counts of what was read and dropped."""
    lines = iter(lines)
    header = next(lines).rstrip("\n").split("\t")
    si = column_index(header, sample_col, 0)
    di = column_index(header, date_col)
    if di is None:
        raise SystemExit("error: --dates needs --date_column")
    need = PRECISION[min_precision]
    stats = {"rows": 0, "empty": 0, "unparsed": 0, "imprecise": 0,
             "conflicting": 0, "precision": [0, 0, 0], "examples": []}
    out, bad = {}, set()
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) <= max(si, di) or not f[si]:
            continue
        stats["rows"] += 1
        txt = f[di].strip()
        if not txt or txt.upper() in ("NA", "N/A", "NONE", "NULL", "MISSING",
                                      "UNKNOWN", "NOT COLLECTED",
                                      "NOT APPLICABLE", "?", "-"):
            stats["empty"] += 1
            continue
        got = parse_date(txt)
        if got is None:
            stats["unparsed"] += 1
            if len(stats["examples"]) < 5 and txt not in stats["examples"]:
                stats["examples"].append(txt)
            continue
        if got[1] < need:
            stats["imprecise"] += 1
            continue
        name = f[si].strip()
        if name in out and abs(out[name] - got[0]) > 1e-9:
            bad.add(name)
        out[name] = got[0]
        stats["precision"][got[1]] += 1
    for name in bad:
        del out[name]
    stats["conflicting"] = len(bad)
    return out, stats


def tip_dates(np, tree, dates):
    """Per node, the tip's decimal year, or nan."""
    a = np.full(len(tree.node_ids), np.nan)
    for i, nid in enumerate(tree.node_ids):
        if not tree.children[i] and nid in dates:
            a[i] = dates[nid]
    return a


def dated_sums(np, tree, R, td):
    """Per node, sums over its dated tips of 1, t, t^2, y, y*t and y^2,
    where t is the date and y the root-to-tip distance (all classes)."""
    n = len(tree.node_ids)
    S = np.zeros((n, 6))
    ok = ~np.isnan(td)
    y = R.sum(1)
    t = np.where(ok, td, 0.0)
    S[ok] = np.column_stack([np.ones(ok.sum()), t[ok], t[ok] ** 2, y[ok],
                             y[ok] * t[ok], y[ok] ** 2])
    for i in range(n - 1, 0, -1):
        p = tree.parent[i]
        if p >= 0 and S[i, 0]:
            S[p] += S[i]
    return S


def clock(np, row):
    """Root-to-tip regression from one row of sums: slope (mutations per
    year), r^2, tips and the date range's standard deviation."""
    n, st, stt, sy, syt, syy = row
    if n < 3:
        return None
    vt = stt / n - (st / n) ** 2
    vy = syy / n - (sy / n) ** 2
    cov = syt / n - st * sy / n ** 2
    if vt <= 1e-12:
        return None
    return {"slope": cov / vt, "r2": cov * cov / (vt * vy) if vy > 0 else 0,
            "n": int(n), "sd_years": math.sqrt(vt)}


def dated_stats(np, nodes, anc, S, R, b, min_n):
    """For each node and its ancestor A, on dated tips only: the mean
    depth from A on each side, their log2 ratio, the same with the split's
    side moved to the outgroup's mean date along the tree's clock b, and
    each side's root-to-tip slope.  nan where a side has fewer than min_n
    dated tips."""
    v = np.asarray(nodes, dtype=np.int64)
    a = np.asarray([anc[i] for i in nodes], dtype=np.int64)
    Ss = S[v]
    So = S[a] - Ss
    ra = R[a].sum(1)
    out = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        for side, X in (("S", Ss), ("out", So)):
            n = X[:, 0]
            tm = X[:, 1] / n
            ym = X[:, 3] / n
            vt = X[:, 2] / n - tm ** 2
            cov = X[:, 4] / n - tm * ym
            out["n_" + side] = n
            out["date_" + side] = tm
            out["depth_" + side] = ym - ra
            out["slope_" + side] = np.where(vt > 1e-12, cov / vt, np.nan)
        nS, nO = out["n_S"], out["n_out"]
        dO = (out["depth_out"] * nO + PSEUDO) / nO
        out["log2_rate"] = np.log2((out["depth_S"] * nS + PSEUDO) / nS / dO)
        adj = out["depth_S"] - b * (out["date_S"] - out["date_out"])
        num = adj * nS + PSEUDO
        out["adj_log2_rate"] = np.where(num > 0, np.log2(num / nS / dO),
                                        np.nan)
        out["slope_diff"] = out["slope_S"] - out["slope_out"]
    bad = (nS < min_n) | (nO < min_n)
    for k in ("log2_rate", "adj_log2_rate", "slope_diff"):
        out[k] = np.where(bad, np.nan, out[k])
    return out


# --------------------------------------------------------------------------
# tests
# --------------------------------------------------------------------------

def bh(ps):
    """Benjamini-Hochberg q values, in input order (None stays None)."""
    idx = [i for i, p in enumerate(ps) if p is not None]
    m = len(idx)
    q = [None] * len(ps)
    prev = 1.0
    for rank, i in reversed(list(enumerate(sorted(idx, key=lambda i: ps[i]),
                                           1))):
        prev = min(prev, ps[i] * m / rank)
        q[i] = prev
    return q


def control_pool(tree, T, anc, split_idx, overlap):
    """Nodes eligible as controls: not a split root, a comparison ancestor
    defined, at least 2 tips, and no split holding `overlap` of its tips."""
    bad = set(split_idx)
    for s in split_idx:
        a = tree.parent[s]
        while a >= 0 and T[s] >= overlap * T[a]:
            bad.add(a)
            a = tree.parent[a]
    return [v for v in range(len(tree.node_ids))
            if v not in bad and anc[v] >= 0 and T[v] >= 2]


def two_sided(np, x, null):
    """p that a draw from null is at least as far from its median as x."""
    if not len(null):
        return None, None
    med = float(np.median(null))
    d = abs(x - med)
    k = int(np.sum(np.abs(null - med) >= d - 1e-12))
    return (1 + k) / (1 + len(null)), med


def control_window(sizes, n, w, k):
    """Indices into the size-sorted controls for a split of n tips: those
    within w fold, or when fewer than k are, the k nearest in log size."""
    lo = bisect.bisect_left(sizes, n / w)
    hi = bisect.bisect_right(sizes, n * w)
    if hi - lo >= k or len(sizes) <= k:
        return (lo, hi, False) if hi - lo >= k else (0, len(sizes), True)
    ln = math.log(n)
    lo = hi = bisect.bisect_left(sizes, n)
    while hi - lo < k:
        if lo == 0:
            hi += 1
        elif hi == len(sizes):
            lo -= 1
        elif ln - math.log(sizes[lo - 1]) <= math.log(sizes[hi]) - ln:
            lo -= 1
        else:
            hi += 1
    return lo, hi, True


def run(np, tree, T, anc, L, split_idx, ctrl, args, cats, td=None):
    """Score splits and controls on one branch set.  Returns per-split rows,
    the controls' (tips, log2_rate, node, adj_log2_rate, slope_diff), and
    with tip dates td the tree's root-to-tip clock."""
    D, R = depth_sums(np, tree, L, T)
    scored = [s for s in split_idx if anc[s] >= 0]
    ms, mo, to = rates(np, scored, anc, T, D, R)
    ts = T[np.asarray(scored, dtype=np.int64)]
    cms, cmo, cto = rates(np, ctrl, anc, T, D, R)
    cts = T[np.asarray(ctrl, dtype=np.int64)]
    ### total over classes, then per class
    lr_tot = log2_rate(np, ms.sum(1, keepdims=True), mo.sum(1, keepdims=True),
                       ts, to)[:, 0]
    lr_cls = log2_rate(np, ms, mo, ts, to)
    c_tot = log2_rate(np, cms.sum(1, keepdims=True),
                      cmo.sum(1, keepdims=True), cts, cto)[:, 0]
    c_cls = log2_rate(np, cms, cmo, cts, cto)
    order = np.argsort(cts, kind="stable")
    sizes = cts[order]
    w = args.control_window
    rows = []
    for k, s in enumerate(scored):
        lo, hi, widened = control_window(sizes, ts[k], w, args.min_controls)
        win = order[lo:hi]
        p, med = two_sided(np, lr_tot[k], c_tot[win])
        ### Shape: does any class's rate change depart from the split's
        ### overall change by more than it does in random clades?
        dev = lr_cls[k] - lr_tot[k]
        ndev = c_cls[win] - c_tot[win][:, None]
        sp = [two_sided(np, dev[j], ndev[:, j])[0] for j in range(len(cats))]
        jm = min(range(len(cats)), key=lambda j: sp[j]) if win.size else 0
        shape_p = min(1.0, len(cats) * sp[jm]) if win.size else None
        d = ms[k] - mo[k]
        tot = float(d.sum())
        absum = float(np.abs(d).sum())
        cls = []
        for j, c in enumerate(cats):
            pj, mj = two_sided(np, lr_cls[k, j], c_cls[win, j])
            cls.append({"class": c, "mean_S": float(ms[k, j]),
                        "mean_out": float(mo[k, j]),
                        "excess": float(d[j]),
                        "share": float(d[j]) / tot if tot else None,
                        "log2_rate": float(lr_cls[k, j]),
                        "ctrl_median": mj, "ctrl_p": pj})
        rows.append({
            "node": tree.node_ids[s], "idx": s, "tips": int(ts[k]),
            "ancestor": tree.node_ids[anc[s]], "out_tips": int(to[k]),
            "mean_S": float(ms[k].sum()), "mean_out": float(mo[k].sum()),
            "excess": tot, "log2_rate": float(lr_tot[k]),
            "mix": 1 - abs(tot) / absum if absum else None,
            "ctrl_n": int(len(win)), "ctrl_median": med,
            "ctrl_widened": widened,
            "shape_p": shape_p, "shape_class": cats[jm],
            "ctrl_tips": (int(sizes[lo]), int(sizes[hi - 1])) if hi > lo
            else None,
            "ctrl_spread": (float(np.percentile(c_tot[win], 97.5)
                                  - np.percentile(c_tot[win], 2.5))
                            if len(win) else None),
            "ctrl_p": p, "classes": cls})
    for r, q in zip(rows, bh([r["ctrl_p"] for r in rows])):
        r["q"] = q
    for j in range(len(cats)):
        qs = bh([r["classes"][j]["ctrl_p"] for r in rows])
        for r, q in zip(rows, qs):
            r["classes"][j]["q"] = q
    ctl = [[t, lr, tree.node_ids[v], None, None]
           for t, lr, v in zip(cts.tolist(), c_tot.tolist(), ctrl)]
    if td is None:
        return rows, ctl, None
    S = dated_sums(np, tree, R, td)
    clk = clock(np, S[0])
    b = clk["slope"] if clk else 0.0
    ds = dated_stats(np, scored, anc, S, R, b, args.min_dated_tips)
    cd = dated_stats(np, ctrl, anc, S, R, b, args.min_dated_tips)
    for c, adj, sd in zip(ctl, cd["adj_log2_rate"].tolist(),
                          cd["slope_diff"].tolist()):
        c[3] = None if math.isnan(adj) else adj
        c[4] = None if math.isnan(sd) else sd

    def val(x):
        x = float(x)
        return None if math.isnan(x) else x

    for k, r in enumerate(rows):
        lo, hi, _ = control_window(sizes, ts[k], w, args.min_controls)
        win = order[lo:hi]
        d = {key: val(ds[key][k]) for key in ds}
        for key, null_key in (("adj", "adj_log2_rate"),
                              ("slope", "slope_diff")):
            null = cd[null_key][win]
            null = null[~np.isnan(null)]
            x = d[null_key]
            p, med = two_sided(np, x, null) if x is not None \
                else (None, None)
            d[key + "_ctrl_n"] = int(len(null))
            d[key + "_ctrl_median"] = med
            d[key + "_p"] = p
        r["dated"] = d
    for key in ("adj", "slope"):
        for r, q in zip(rows, bh([r["dated"][key + "_p"] for r in rows])):
            r["dated"][key + "_q"] = q
    return rows, ctl, clk


def passes(r, max_q):
    return r is not None and r["q"] is not None and r["q"] <= max_q


def verdict(a, i, max_q):
    if passes(a, max_q) and passes(i, max_q) and \
            (a["log2_rate"] > 0) == (i["log2_rate"] > 0):
        return "faster" if a["log2_rate"] > 0 else "slower"
    if passes(a, max_q) and not passes(i, max_q):
        return "terminal only"
    return "not different"


def class_calls(a, i, max_q):
    out = []
    for ca, ci in zip(a["classes"], i["classes"]):
        if ca["q"] is not None and ca["q"] <= max_q and ci["q"] is not None \
                and ci["q"] <= max_q and \
                (ca["log2_rate"] > 0) == (ci["log2_rate"] > 0):
            out.append(f"{ca['class']}{'+' if ca['log2_rate'] > 0 else '-'}")
    return ",".join(out)


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def fmt(v, spec=".3g"):
    if v is None:
        return "-"
    if isinstance(v, float) and spec == "p":
        return f"{v:.2g}"
    return format(v, spec)


def top_excess(r, n=3):
    """The classes whose depth differs most, as 'class+0.08': the class's
    difference mean_S - mean_out over the outgroup's total mean depth, so
    the values add up to the total rate change and compare between
    splits."""
    if not r["mean_out"]:
        return ""
    cls = sorted(r["classes"], key=lambda c: -abs(c["excess"]))[:n]
    return ",".join(f"{c['class']}{c['excess'] / r['mean_out']:+.2f}"
                    for c in cls)


def report_run(rows, label, capped, out):
    print(f"== {label} ==", file=out)
    if capped:
        print(f"  {capped} branches over --max_branch_length left out",
              file=out)
    hdr = ("node", "tips", "ancestor", "out_tips", "mean_S", "mean_out",
           "log2_rate", "ctrl_n", "ctrl_med", "ctrl_95%", "ctrl_p", "q",
           "mix")
    print("\t".join(hdr), file=out)
    for r in sorted(rows, key=lambda r: (r["ctrl_p"] is None, r["ctrl_p"])):
        print("\t".join([
            r["node"], str(r["tips"]), r["ancestor"], str(r["out_tips"]),
            fmt(r["mean_S"]), fmt(r["mean_out"]), f"{r['log2_rate']:+.3f}",
            str(r["ctrl_n"]), fmt(r["ctrl_median"], "+.3f"),
            fmt(r["ctrl_spread"], ".3f"), fmt(r["ctrl_p"], "p"),
            fmt(r["q"], "p"), fmt(r["mix"], ".2f")]), file=out)
    print(file=out)


def report_clock(clk, label, out):
    if clk is None:
        print(f"  root-to-tip clock, {label}: fewer than 3 dated tips or "
              f"one date; no date adjustment", file=out)
        return
    print(f"  root-to-tip clock, {label}: {clk['slope']:.4g} mutations per "
          f"year, r^2 {clk['r2']:.3f}, {clk['n']} dated tips, date SD "
          f"{clk['sd_years']:.1f} years", file=out)
    if clk["slope"] <= 0:
        print("    WARNING: the slope is not positive, so the dates carry no "
              "clock signal on this tree and the date adjustment is "
              "meaningless", file=out)


def report_dated(rows, label, clk, out):
    print(f"== {label}, dated tips only ==", file=out)
    report_clock(clk, label, out)
    hdr = ("node", "dated_S", "dated_out", "date_S", "date_out",
           "log2_rate", "adj_log2", "ctrl_n", "ctrl_med", "ctrl_p", "q",
           "slope_S", "slope_out", "ctrl_n", "ctrl_p", "q")
    print("\t".join(hdr), file=out)
    key = lambda r: (r["dated"]["adj_p"] is None, r["dated"]["adj_p"] or 0)
    for r in sorted(rows, key=key):
        d = r["dated"]
        print("\t".join([
            r["node"], str(int(d["n_S"] or 0)), str(int(d["n_out"] or 0)),
            fmt(d["date_S"], ".1f"), fmt(d["date_out"], ".1f"),
            fmt(d["log2_rate"], "+.3f"), fmt(d["adj_log2_rate"], "+.3f"),
            str(d["adj_ctrl_n"]), fmt(d["adj_ctrl_median"], "+.3f"),
            fmt(d["adj_p"], "p"), fmt(d["adj_q"], "p"),
            fmt(d["slope_S"], ".3g"), fmt(d["slope_out"], ".3g"),
            str(d["slope_ctrl_n"]), fmt(d["slope_p"], "p"),
            fmt(d["slope_q"], "p")]), file=out)
    print(file=out)


def dated_verdict(a, max_q):
    """faster, slower or not different on the all-branch date-adjusted
    rate, or '-' where it wasn't scored."""
    d = a.get("dated")
    if not d or d["adj_q"] is None:
        return "-"
    if d["adj_q"] <= max_q:
        return "faster" if d["adj_log2_rate"] > 0 else "slower"
    return "not different"


def summarize(all_rows, int_rows, args, comp, out):
    by_int = {r["node"]: r for r in int_rows}
    table = []
    for a in all_rows:
        i = by_int.get(a["node"])
        if i is None:
            continue
        table.append({
            "node": a["node"], "tips": a["tips"],
            "log2_rate": a["log2_rate"], "q": a["q"],
            "log2_rate_internal": i["log2_rate"], "q_internal": i["q"],
            "mix": a["mix"], "mix_internal": i["mix"],
            "verdict": verdict(a, i, args.max_q),
            "classes": class_calls(a, i, args.max_q),
            "excess_classes": top_excess(a),
            "clades": clade_label(comp[a["node"]])
            if comp and a["node"] in comp else ""})
        if "dated" in a:
            da = a["dated"]
            table[-1].update({
                "dated_tips": da["n_S"] and int(da["n_S"]),
                "dated_out": da["n_out"] and int(da["n_out"]),
                "adj_log2_rate": da["adj_log2_rate"], "q_adj": da["adj_q"],
                "slope_S": da["slope_S"], "slope_out": da["slope_out"],
                "q_slope": da["slope_q"],
                "date_verdict": dated_verdict(a, args.max_q)})
    order = {"faster": 0, "slower": 1, "terminal only": 2, "not different": 3}
    table.sort(key=lambda t: (order[t["verdict"]], -abs(t["log2_rate"])))
    print("== SUMMARY: rate against sister clades, both runs ==", file=out)
    print(f"verdict at q <= {args.max_q}; classes: own rate differs on both "
          f"runs; excess: each class's all-branch difference over the outgroup's mean depth", file=out)
    dated = bool(table) and "date_verdict" in table[0]
    if dated:
        print("dated, all branches only: dated tips on each side, the "
              "date-adjusted log2_rate and its q, root-to-tip slopes in "
              "mutations per year and the q of their difference; "
              "date_verdict from q_adj alone", file=out)
    hdr = ("node", "tips", "log2_rate", "q", "log2_int", "q_int", "mix",
           "verdict", "classes", "excess")
    if dated:
        hdr += ("dated_S", "dated_out", "adj_log2", "q_adj", "date_verdict",
                "slope_S", "slope_out", "q_slope")
    print("\t".join(hdr + ("clades",)), file=out)
    for t in table:
        f = [t["node"], str(t["tips"]), f"{t['log2_rate']:+.3f}",
             fmt(t["q"], "p"), f"{t['log2_rate_internal']:+.3f}",
             fmt(t["q_internal"], "p"), fmt(t["mix"], ".2f"), t["verdict"],
             t["classes"] or "-", t["excess_classes"]]
        if dated:
            f += [str(t["dated_tips"] or 0), str(t["dated_out"] or 0),
                  fmt(t["adj_log2_rate"], "+.3f"), fmt(t["q_adj"], "p"),
                  t["date_verdict"], fmt(t["slope_S"], ".3g"),
                  fmt(t["slope_out"], ".3g"), fmt(t["q_slope"], "p")]
        print("\t".join(f + [t["clades"]]), file=out)
    counts = {}
    for t in table:
        counts[t["verdict"]] = counts.get(t["verdict"], 0) + 1
    print("\n" + ", ".join(f"{v}: {counts.get(v, 0)}" for v in order),
          file=out)
    if dated:
        dc = {}
        for t in table:
            dc[t["date_verdict"]] = dc.get(t["date_verdict"], 0) + 1
        print("date-adjusted: " + ", ".join(
            f"{v}: {dc.get(v, 0)}" for v in ("faster", "slower",
                                              "not different", "-")),
            file=out)
    return table


def read_composition(path):
    with open(path) as fh:
        return {r["node"]: r for r in csv.DictReader(fh, delimiter="\t")}


def clade_label(r):
    parts = []
    for item in r.get("composition", "").split(","):
        if ":" in item:
            name, n = item.rsplit(":", 1)
            parts.append(f"{name} {n}")
    label = ", ".join(parts[:2])
    flags = r.get("flags", "")
    return f"{label} ({flags})" if flags else label


def write_outputs(table, runs, controls, cats, prefix):
    paths = [f"{prefix}.tsv", f"{prefix}_classes.tsv",
             f"{prefix}_controls.tsv"]
    cols = ("node", "tips", "log2_rate", "q", "log2_rate_internal",
            "q_internal", "mix", "mix_internal", "verdict", "classes",
            "excess_classes")
    if table and "date_verdict" in table[0]:
        cols += ("dated_tips", "dated_out", "adj_log2_rate", "q_adj",
                 "date_verdict", "slope_S", "slope_out", "q_slope")
    cols += ("reading", "shape_p", "concentration", "clades", "reason")
    with open(paths[0], "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for t in table:
            fh.write("\t".join("" if t.get(c) is None else
                               f"{t[c]:.6g}" if isinstance(t[c], float)
                               else str(t[c]) for c in cols) + "\n")
    with open(paths[1], "w") as fh:
        fh.write("node\tbranches\tclass\tmean_S\tmean_out\texcess\tshare\t"
                 "log2_rate\tctrl_median\tctrl_p\tq\n")
        for label, rows in runs:
            for r in rows:
                for c in r["classes"]:
                    fh.write("\t".join([r["node"], label, c["class"]] + [
                        "" if c[k] is None else f"{c[k]:.6g}"
                        for k in ("mean_S", "mean_out", "excess", "share",
                                  "log2_rate", "ctrl_median", "ctrl_p",
                                  "q")]) + "\n")
    with open(paths[2], "w") as fh:
        fh.write("branches\tnode\ttips\tlog2_rate\tadj_log2_rate\t"
                 "slope_diff\n")
        for label, ctl in controls:
            for tips, lr, node, adj, sd in ctl:
                fh.write(f"{label}\t{node}\t{tips}\t{lr:.6g}\t"
                         f"{'' if adj is None else f'{adj:.6g}'}\t"
                         f"{'' if sd is None else f'{sd:.6g}'}\n")
    return paths


# --------------------------------------------------------------------------
# reading: HOWTO's decision rules, applied
# --------------------------------------------------------------------------

def min_controls(max_q):
    """Fewest controls with which p = 1/(n+1) can reach max_q."""
    return max(1, math.ceil(1 / max_q) - 1)


def fold_text(lr):
    """A log2 ratio as a percentage change: +0.5 -> '41% faster'."""
    f = 2 ** lr
    if f >= 1:
        return f"{(f - 1) * 100:.0f}% faster"
    return f"{(1 - f) * 100:.0f}% slower"


def concentration(r):
    """How far a rate change's excess departs from the sisters' spectrum:
    half the summed |share of the excess - share of the spectrum| over the
    classes.  0 when every class gained in proportion (a genome-wide rate
    change), up to 1 when the excess lies in classes the sisters lack."""
    tot = r["excess"]
    out_tot = r["mean_out"]
    if not tot or not out_tot:
        return None
    return 0.5 * sum(abs(c["excess"] / tot - c["mean_out"] / out_tot)
                     for c in r["classes"])


def trade_text(r):
    """The classes a split gained and lost, for a redistribution."""
    gained = [c["class"] for c in sorted(r["classes"],
                                          key=lambda c: -c["excess"])
              if c["excess"] > 0][:2]
    lost = [c["class"] for c in sorted(r["classes"],
                                        key=lambda c: c["excess"])
            if c["excess"] < 0][:2]
    return (f"more {' and '.join(gained) or 'nothing'}, less "
            f"{' and '.join(lost) or 'nothing'}")


def specific(a, max_q):
    """Is a rate change class-specific: does one class's change depart from
    the overall change by more than in random clades of the same size?"""
    return a.get("shape_p") is not None and a["shape_p"] <= max_q


def shape_text(a):
    c = next(c for c in a["classes"] if c["class"] == a["shape_class"])
    return (f"{c['class']} changed {fold_text(c['log2_rate'])} against "
            f"{fold_text(a['log2_rate'])} overall, further apart than in "
            f"random clades (shape p {a['shape_p']:.2g}, over "
            f"{len(a['classes'])} classes)")


def classify(a, i, args, depth_a=None):
    """One split's reading: (category, reason).  a and i are its rows from
    the all-branch and internal-branch runs; a may carry 'dated'."""
    Q = args.max_q
    need = min_controls(Q)
    lr, p, q = a["log2_rate"], a["ctrl_p"], a["q"]
    near = (depth_a is not None and depth_a <= args.near_root)
    root_note = (f" Its comparison ancestor is {depth_a} node"
                 f"{'s' if depth_a != 1 else ''} below the root, so rooting "
                 f"and how well each lineage maps to the reference can move "
                 f"it." if near else "")
    widened = ""
    if a.get("ctrl_widened") and a.get("ctrl_tips"):
        lo_t, hi_t = a["ctrl_tips"]
        widened = (f" Few clades are within {args.control_window:g}x of its "
                   f"size, so it was compared with the {a['ctrl_n']} nearest "
                   f"({lo_t} to {hi_t} tips).")
    if a["ctrl_n"] < need:
        return ("untestable",
                f"only {a['ctrl_n']} clade{'s' if a['ctrl_n'] != 1 else ''} "
                f"on the tree {'is' if a['ctrl_n'] == 1 else 'are'} within "
                f"{args.control_window:g}x of its size, and it takes "
                f"{need} for p to reach {Q:g}. Its sister-clade ratio is "
                f"{fold_text(lr)} (all branches) and {fold_text(i['log2_rate'])} "
                f"(internal), with nothing to judge it against.")
    v = verdict(a, i, Q)
    d = a.get("dated")
    dated_ok = (d is not None and d.get("adj_q") is not None
                and d["adj_ctrl_n"] >= need)
    dv = dated_verdict(a, Q) if dated_ok else None
    half = (a["ctrl_spread"] or 0) / 2
    if v in ("faster", "slower"):
        if dv == "not different":
            return ("sampling time",
                    f"{v} than its sisters ({fold_text(lr)}, q {q:.2g}), but "
                    f"not once collection dates are allowed for "
                    f"(adjusted {fold_text(d['adj_log2_rate'])}, q "
                    f"{d['adj_q']:.2g}; mean dates {d['date_S']:.1f} against "
                    f"{d['date_out']:.1f}). It was sampled at a different "
                    f"time, not mutating at a different rate.")
        both = (f"{fold_text(lr)} on all branches (q {q:.2g}) and "
                f"{fold_text(i['log2_rate'])} on internal ones (q "
                f"{i['q']:.2g})")
        if dated_ok and dv == v:
            both += (f", and {fold_text(d['adj_log2_rate'])} after "
                     f"allowing for collection dates (q {d['adj_q']:.2g})")
        elif d is not None and not dated_ok:
            both += ", with too few dated tips to check against dates"
        if specific(a, Q):
            cat = "mutator" if v == "faster" else "antimutator"
            return (cat,
                    f"{v}: {both}. The change is class-specific: "
                    f"{shape_text(a)}. A candidate for step 12."
                    + root_note + widened)
        return (f"genome-wide {v}",
                f"{v}: {both}. No class departs from that overall change "
                f"more than in random clades (shape p "
                f"{fmt(a['shape_p'], 'p')}), which a repair defect would "
                f"not do. Suspect sampling dates"
                f"{'' if dated_ok else ' (not checked)'}, selection, or "
                f"tree structure." + root_note + widened)
    if v == "terminal only":
        return ("terminal only",
                f"{fold_text(lr)} on all branches (q {q:.2g}), but not on "
                f"internal branches ({fold_text(i['log2_rate'])}, q "
                f"{i['q']:.2g}). The difference is on terminal branches: "
                f"per-genome sequencing error or a noisier study, not "
                f"mutation.")
    ### not different: suggestive if p passes on both runs before correction
    if (p is not None and i["ctrl_p"] is not None and p <= Q
            and i["ctrl_p"] <= Q and (lr > 0) == (i["log2_rate"] > 0)):
        if specific(a, Q):
            kind = (f"class-specific: {shape_text(a)}; the pattern of "
                    f"{'a ' if lr > 0 else 'an anti'}mutator, short of "
                    f"significance")
        else:
            kind = (f"spread over the classes as in random clades (shape p "
                    f"{fmt(a['shape_p'], 'p')}), so not a repair defect")
        return (f"suggestive {'faster' if lr > 0 else 'slower'}",
                f"{fold_text(lr)} (p {p:.2g}, internal p "
                f"{i['ctrl_p']:.2g}) but q {q:.2g} after correcting for the "
                f"splits tested. The change is {kind}." + root_note
                + widened)
    res = f"resolves ±{fold_text(half).split()[0]}"
    if widened:
        res += f", against the {a['ctrl_n']} clades nearest its size"
    if a["mix"] is not None and a["mix"] >= 0.5:
        return ("unchanged rate",
                f"classes traded: {trade_text(a)} ({fold_text(lr)}, p "
                f"{p:.2g}, internal p {i['ctrl_p']:.2g}; {res}).")
    return ("unchanged rate", f"{fold_text(lr)}, p {p:.2g}, internal p "
                              f"{i['ctrl_p']:.2g}; {res}.")


HEADERS = {
    "unchanged rate": "Rate within the range of random clades of the same "
                      "size, or outside it on one run only: the change "
                      "against the sisters, its p on all and on internal "
                      "branches, and the "
                      "smallest change detectable at that size (half the "
                      "middle 95% of the controls). 'classes traded' marks "
                      "splits whose classes moved in opposite directions "
                      "(mix 0.5 or more): a spectrum shift at an unchanged "
                      "rate.",
    "suggestive faster": "p at or below {q:g} on both runs, but not q.",
    "suggestive slower": "p at or below {q:g} on both runs, but not q.",
}

ORDER = ["mutator", "antimutator", "genome-wide faster", "genome-wide slower",
         "sampling time", "terminal only", "suggestive faster",
         "suggestive slower", "unchanged rate", "untestable"]


def write_reading(table, runs, args, clocks, date_info, depth_a, out):
    """The report's first section: each split's category and why, then what
    the run as a whole says."""
    by_a = {r["node"]: r for r in runs[0][1]}
    by_i = {r["node"]: r for r in runs[1][1]}
    cats = {}
    for t in table:
        a, i = by_a[t["node"]], by_i[t["node"]]
        cat, why = classify(a, i, args, depth_a.get(t["node"]))
        t["reading"], t["reason"] = cat, why
        t["concentration"] = concentration(a)
        t["shape_p"] = a.get("shape_p")
        cats.setdefault(cat, []).append(t)
    wrap = textwrap.TextWrapper(width=78, subsequent_indent="    ")

    def para(text, indent=""):
        print(textwrap.fill(text, width=78, initial_indent=indent,
                            subsequent_indent=indent), file=out)

    print("=" * 78, file=out)
    print("READING: is each split a mutator, a spectrum-shifter, or both?",
          file=out)
    print("=" * 78, file=out)
    n = len(table)
    counts = ", ".join(f"{len(cats[c])} {c}" for c in ORDER if c in cats)
    para(f"{n} splits scored against their sister clades: {counts}.")
    print(file=out)

    ### the conclusion first
    movers = [t for c in ("mutator", "antimutator") for t in cats.get(c, [])]
    rate = [t for c in ("genome-wide faster", "genome-wide slower")
            for t in cats.get(c, [])]
    if movers:
        para("Conclusion: " + ", ".join(
            f"{t['node']} ({t['clades'] or 'unlabelled'}) is a class-specific "
            f"{t['reading']}" for t in movers) + ". Take "
            f"{'it' if len(movers) == 1 else 'them'} to step 12 to look for "
            f"the gene.")
    else:
        tested = n - len(cats.get("untestable", []))
        para(f"Conclusion: no split is a detectable class-specific mutator. "
             f"Of the {tested} that could be tested, "
             f"{len(cats.get('unchanged rate', []))} shifted their spectrum "
             f"at a rate within the range of random clades of their size"
             + (f", and {len(rate)} changed rate across every class alike, "
                f"which is not a repair defect" if rate else "")
             + ". Their spectra are redistributions, not added mutations, as "
             "far as this tree can resolve.")
        sugg = [t for c in ("suggestive faster", "suggestive slower")
                for t in cats.get(c, [])]
        near = [t for t in sugg if specific(by_a[t["node"]], args.max_q)]
        flat = [t for t in sugg if t not in near]
        if flat:
            rooty = [t for t in flat
                     if (depth_a.get(t["node"]) or 99) <= args.near_root]
            print(file=out)
            para(f"The other {len(flat)} differ in rate before correcting for the "
                 f"{n} splits tested ("
                 + ", ".join(f"{t['node']} {t['reading'].split()[1]}"
                             for t in flat)
                 + "), but none in a class-specific way: each spreads its "
                 "change over the classes as random clades do, which is not "
                 "what a repair defect does."
                 + (f" {len(rooty)} of them ("
                    + ", ".join(t["node"] for t in rooty)
                    + ") are compared next to the root, where rooting and "
                    "reference mapping can move the ratio."
                    if rooty else ""))
        if near:
            print(file=out)
            para(f"Short of significance, {len(near)} split"
                 f"{'s' if len(near) != 1 else ''} show"
                 f"{'' if len(near) != 1 else 's'} the class-specific "
                 f"pattern a mutator (or, slower, an antimutator) would "
                 f"make: "
                 + ", ".join(f"{t['node']} ({t['clades'] or 'unlabelled'})"
                             for t in near)
                 + ". They are the ones to recheck on a larger tree.")
    print(file=out)

    ### resolution
    widths = sorted(r["ctrl_spread"] / 2 for r in runs[0][1]
                    if r["ctrl_spread"] and r["ctrl_n"] >= min_controls(
                        args.max_q))
    if widths:
        mid = widths[len(widths) // 2]
        para(f"Resolution: random clades of a split's size vary by "
             f"±{widths[0]:.2f} to ±{widths[-1]:.2f} log2 (median "
             f"±{mid:.2f}, i.e. {fold_text(mid).split()[0]}) in their ratio "
             f"to their sisters. A rate change smaller than that can't be "
             f"told from chance here, so 'unchanged rate' means 'unchanged "
             f"within about {fold_text(mid).split()[0]}'.")
        print(file=out)

    ### dates
    if date_info is not None:
        clk = clocks[0]
        pr = date_info["precision"]
        kept = sum(pr)
        line = (f"Dates: {date_info['dated']} of {date_info['tips']} tips "
                f"dated ({100 * date_info['dated'] / date_info['tips']:.0f}%)"
                f"; {100 * pr[0] / kept:.0f}% of the dates kept are a year "
                f"only. " if kept else "Dates: none kept. ")
        if clk is None:
            line += "Too few dated tips for a clock; dates were not used."
        else:
            line += (f"The tree's root-to-tip clock is {clk['slope']:.3g} "
                     f"mutations per year with r² {clk['r2']:.3f}, over "
                     f"dates with an SD of {clk['sd_years']:.1f} years. ")
            if clk["slope"] <= 0:
                line += ("A slope at or below zero means the dates carry no "
                         "clock signal, and the adjustment did nothing. ")
            elif clk["r2"] < 0.05:
                line += ("Dates explain almost none of the variation in "
                         "root-to-tip distance, which lineage depth "
                         "dominates, so the slope is a rough clock; that is "
                         "enough to adjust for sampling-time differences of "
                         "a few years, not to measure rates. ")
            shifts = [(abs(r["dated"]["adj_log2_rate"]
                           - r["dated"]["log2_rate"]), r)
                      for r in runs[0][1]
                      if r.get("dated") and r["dated"]["adj_log2_rate"]
                      is not None]
            changed = [t for t in table if t["reading"] == "sampling time"]
            undated = [t["node"] for t in table
                       if t.get("date_verdict") == "-"]
            if shifts:
                big, br = max(shifts, key=lambda x: x[0])
                line += (f"Allowing for dates moves a split's ratio by at "
                         f"most {big:.3f} log2 ({br['node']}, mean dates "
                         f"{br['dated']['date_S']:.1f} against "
                         f"{br['dated']['date_out']:.1f}), and changes the "
                         f"reading of {len(changed)} split"
                         f"{'s' if len(changed) != 1 else ''}"
                         + (": " + ", ".join(t["node"] for t in changed)
                            if changed else "") + ". ")
            if undated:
                line += (f"{len(undated)} split"
                         f"{'s' if len(undated) != 1 else ''} had too few "
                         f"dated tips on one side to adjust ("
                         f"{', '.join(undated)}). ")
            passing = [r["node"] for r in runs[0][1] if r.get("dated")
                       and r["dated"].get("slope_q") is not None
                       and r["dated"]["slope_q"] <= args.max_q]
            line += ("No split's own root-to-tip slope differs from its "
                     "sisters' beyond the controls (these slopes rest on "
                     "few years of dates and are a check only)."
                     if not passing else
                     f"Own root-to-tip slopes differ from the sisters' for "
                     f"{', '.join(passing)}; read them as a check on the "
                     f"adjusted ratio, as they rest on few years of dates.")
        para(line)
        print(file=out)

    ### per split, by category
    for c in ORDER:
        if c not in cats:
            continue
        print(f"-- {c} ({len(cats[c])}) --", file=out)
        if c in HEADERS:
            para(HEADERS[c].format(q=args.max_q), "  ")
        for t in cats[c]:
            head = (f"{t['node']} ({t['clades'] or 'unlabelled'}; "
                    f"{t['tips']} tips): ")
            print(wrap.fill(head + t["reason"]), file=out)
        print(file=out)
    if "untestable" in cats:
        para("Untestable splits had fewer controls than it takes for p to "
             "reach --max_q even after widening to the --min_controls "
             "nearest in size: the tree has too few clades.")
        print(file=out)
    print("=" * 78, file=out)
    print("DETAILS", file=out)
    print("=" * 78, file=out)
    print(file=out)


def node_depths(tree):
    depth = [0] * len(tree.node_ids)
    for i in range(1, len(tree.node_ids)):
        p = tree.parent[i]
        depth[i] = depth[p] + 1 if p >= 0 else 0
    return depth


def analyse(np, tree, split_ids, args, out, comp=None, dates=None,
            dates_note=None, date_stats=None):
    final_out = out
    out = io.StringIO()
    if args.fold:
        fold_map, cats = folded_index()
    else:
        fold_map, cats = None, list(MUT_TYPES)
    ncat = len(cats)
    index = {nid: i for i, nid in enumerate(tree.node_ids)}
    missing = [s for s in split_ids if s not in index]
    if missing:
        print(f"  WARNING: {len(missing)} split root(s) not in this tree, "
              f"e.g. {missing[0]}; node IDs are assigned when a tree is "
              f"read, so use the tree the search ran on.", file=sys.stderr)
    split_idx = [index[s] for s in split_ids if s in index]
    T = tip_counts(np, tree)
    anc = outgroup_ancestors(tree, T, args.min_outgroup_tips)
    no_anc = [tree.node_ids[s] for s in split_idx if anc[s] < 0]
    ctrl = control_pool(tree, T, anc, split_idx, args.control_overlap)
    print(f"{len(split_idx)} splits, {int(T[0])} tips, "
          f"{len(ctrl)} control clades; outgroups of at least "
          f"{args.min_outgroup_tips} tips, controls within "
          f"{args.control_window}x of a split's tips, {len(cats)} classes",
          file=out)
    if no_anc:
        print(f"  not scored (no ancestor with enough outgroup tips): "
              f"{', '.join(no_anc)}", file=out)
    td = None
    if dates is not None:
        td = tip_dates(np, tree, dates)
        for line in dates_note or ():
            print(line, file=out)
        nd = int((~np.isnan(td)).sum())
        print(f"dates: {nd} of {int(T[0])} tips dated; a side is scored "
              f"with at least {args.min_dated_tips} dated tips", file=out)
        if nd:
            q = np.percentile(td[~np.isnan(td)], [0, 25, 50, 75, 100])
            print("  dates min, quartiles, max: "
                  + ", ".join(f"{x:.1f}" for x in q), file=out)
    print(file=out)
    runs, controls = [], []
    clocks = []
    for internal_only, label in ((False, "all"), (True, "internal")):
        L, capped = branch_counts(np, tree, ncat, fold_map,
                                  args.max_branch_length, internal_only)
        ### Dates only on all branches: a tip's date says how long its whole
        ### path ran, and dropping the terminal branch drops the part of the
        ### path that the sampling date governs.
        rows, ctl, clk = run(np, tree, T, anc, L, split_idx, ctrl, args,
                             cats, None if internal_only else td)
        report_run(rows, f"{label} branches", capped, out)
        if td is not None and not internal_only:
            report_dated(rows, f"{label} branches", clk, out)
        runs.append((label, rows))
        controls.append((label, ctl))
        clocks.append(clk)
    table = summarize(runs[0][1], runs[1][1], args, comp, out)
    depth = node_depths(tree)
    depth_a = {tree.node_ids[s]: depth[anc[s]] for s in split_idx
               if anc[s] >= 0}
    date_info = None
    if td is not None:
        date_info = {"dated": int((~np.isnan(td)).sum()), "tips": int(T[0]),
                     "precision": (date_stats or {}).get("precision",
                                                         [0, 0, 0])}
    write_reading(table, runs, args, clocks, date_info, depth_a, final_out)
    final_out.write(out.getvalue())
    return table, runs, controls, cats


# --------------------------------------------------------------------------
# selftest
# --------------------------------------------------------------------------

class _Tree(object):
    def __init__(self):
        self.node_ids, self.parent, self.children = [], [], []
        self.raw_len, self.types = [], []

    def add(self, name, parent, types):
        i = len(self.node_ids)
        self.node_ids.append(name)
        self.parent.append(parent)
        self.children.append([])
        if parent >= 0:
            self.children[parent].append(i)
        self.types.append(list(types))
        self.raw_len.append(len(types))
        return i


def selftest():
    import random
    try:
        import numpy as np
    except ImportError:
        print("error: numpy is required", file=sys.stderr)
        return 1
    fails = []

    def check(ok, label):
        print(f"  {'ok' if ok else 'FAIL'}: {label}")
        if not ok:
            fails.append(label)

    ### Depth sums against brute force on a small tree.
    t = _Tree()
    r = t.add("r", -1, [])
    a = t.add("a", r, [0, 0])
    b = t.add("b", r, [1])
    t.add("a1", a, [0])
    t.add("a2", a, [1, 1, 1])
    t.add("b1", b, [])
    t.add("b2", b, [0])
    t.add("b3", b, [0, 1])
    L, _ = branch_counts(np, t, 12, None, None, False)
    T = tip_counts(np, t)
    D, R = depth_sums(np, t, L, T)
    check(list(T[:3]) == [5, 2, 3], "tip counts")
    ### from r: a1 = 2+1 = 3, a2 = 2+3 = 5, b1 = 1, b2 = 2, b3 = 3
    check(D[0].sum() == 14 and D[1].sum() == 4 and D[2].sum() == 3,
          "depth sums")
    anc = outgroup_ancestors(t, T, 1)
    ms, mo, to = rates(np, [a], anc, T, D, R)
    check(anc[a] == r and abs(ms[0].sum() - 4) < 1e-9
          and abs(mo[0].sum() - 2) < 1e-9 and to[0] == 3,
          "clade a: mean depth 4 against 2 from the root")
    anc3 = outgroup_ancestors(t, T, 3)
    check(anc3[3] == r, "a1 walks up to the root for 3 outgroup tips")
    Li, _ = branch_counts(np, t, 12, None, None, True)
    Di, Ri = depth_sums(np, t, Li, T)
    msi, moi, _ = rates(np, [a], anc, T, Di, Ri)
    check(abs(msi[0].sum() - 2) < 1e-9 and abs(moi[0].sum() - 1) < 1e-9,
          "internal only drops terminal branches")
    Lc, capped = branch_counts(np, t, 12, None, 2, False)
    check(capped == 1 and Lc[4].sum() == 0, "cap drops the 3-mutation branch")

    ### A random tree with one mutator clade (class 0 tripled) and one
    ### spectrum-shifter (class 0 for class 1, same rate).
    rng = random.Random(3)
    nrng = np.random.default_rng(3)
    t = _Tree()
    root = t.add("node_1", -1, [])
    base = [0.2, 0.2] + [0.6 / 10] * 10

    def draw(n, w):
        return rng.choices(range(12), weights=w, k=n)

    def grow(parent, depth, w, name=None):
        k = 0
        frontier = [parent]
        nodes = []
        for d in range(depth):
            nxt = []
            for p in frontier:
                for c in range(2):
                    k += 1
                    m = sum(w) / sum(base)
                    i = t.add(name if (name and d == 0 and c == 0)
                              else f"n{len(t.node_ids)}", p,
                              draw(int(nrng.poisson(6 * m)), w))
                    nxt.append(i)
            frontier = nxt
            nodes += nxt
        return frontier

    mut = [1.0, 0.2] + [0.6 / 10] * 10
    shift = [0.35, 0.05] + [0.6 / 10] * 10
    hubs = []
    for h in range(40):
        hubs.append(t.add(f"hub{h}", root, draw(3, base)))
    ### Each split hangs from a hub, beside a normal sister of the same size.
    ms_ = t.add("node_mut", hubs[0], draw(4, mut))
    grow(ms_, 6, mut)
    grow(t.add("sis0", hubs[0], draw(4, base)), 6, base)
    ss_ = t.add("node_shift", hubs[1], draw(4, shift))
    grow(ss_, 6, shift)
    grow(t.add("sis1", hubs[1], draw(4, base)), 6, base)
    for h in hubs[2:]:
        grow(h, 7, base)

    class A(object):
        fold = False
        max_branch_length = None
        min_outgroup_tips = 10
        control_window = 2.0
        control_overlap = 0.5
        min_controls = 40
        max_q = 0.05
        near_root = 5

    import io
    buf = io.StringIO()
    table, runs, _, _ = analyse(np, t, ["node_mut", "node_shift"], A(), buf)
    v = {x["node"]: x for x in table}
    check(v["node_mut"]["verdict"] == "faster",
          f"mutator called faster ({v['node_mut']['verdict']}, "
          f"log2 {v['node_mut']['log2_rate']:+.2f})")
    check(v["node_mut"]["classes"].startswith("AC+"),
          f"mutator's class is AC ({v['node_mut']['classes']})")
    check(v["node_shift"]["verdict"] == "not different",
          f"shifter called not different ({v['node_shift']['verdict']}, "
          f"log2 {v['node_shift']['log2_rate']:+.2f})")
    check(v["node_shift"]["mix"] > 0.5 > v["node_mut"]["mix"],
          f"mix: shifter {v['node_shift']['mix']:.2f} high, mutator "
          f"{v['node_mut']['mix']:.2f} low")
    check(all(r["ctrl_n"] > 10 for _, rows in runs for r in rows),
          "both splits have controls")
    check(v["node_mut"]["reading"] == "mutator"
          and "AC" in v["node_mut"]["reason"],
          f"reading: the mutator is a class-specific mutator in AC "
          f"({v['node_mut']['reading']})")
    check(v["node_shift"]["reading"] == "unchanged rate"
          and "classes traded" in v["node_shift"]["reason"],
          f"reading: the shifter traded classes at an unchanged rate "
          f"({v['node_shift']['reading']}: {v['node_shift']['reason']})")
    check(buf.getvalue().startswith("=" * 78 + "\nREADING")
          and "node_mut" in buf.getvalue().split("Conclusion:")[1]
          .split("\n\n")[0],
          "the report opens with the reading, and its conclusion names "
          "the mutator")
    ### classify() on hand-made rows: too few controls, terminal only
    cls = [{"class": "x", "excess": 1.0, "mean_out": 5.0}]
    base = {"log2_rate": 0.5, "ctrl_p": 0.5, "q": 0.5, "ctrl_n": 3,
            "ctrl_spread": 0.5, "mix": 0.0, "excess": 1.0, "mean_out": 5.0,
            "classes": cls}
    check(classify(base, dict(base), A())[0] == "untestable",
          "reading: 3 controls is untestable at q 0.05")
    a_ = dict(base, ctrl_n=100, ctrl_p=0.001, q=0.01)
    check(classify(a_, dict(base, ctrl_n=100), A())[0] == "terminal only",
          "reading: passes on all branches only is terminal only")
    check(control_window([1, 2, 3, 50, 60, 70, 80], 4, 2.0, 3)
          == (0, 3, True)
          and control_window([1, 2, 3, 50, 60, 70, 80], 60, 2.0, 3)
          == (3, 7, False),
          "control window widens to the nearest only when short")
    check(abs(bh([0.01, 0.04, 0.03])[1] - 0.04) < 1e-12, "bh")
    _selftest_dates(np, check, A)
    print("selftest:", "FAILED" if fails else "passed")
    return 1 if fails else 0


def _selftest_dates(np, check, Args):
    pd = parse_date
    check(pd("2010") == (2010.5, 0) and pd("2010-XX-XX") == (2010.5, 0),
          "a bare year is its middle")
    m = pd("2010-07")
    check(m[1] == 1 and abs(m[0] - (2010 + 196.5 / 365)) < 1e-9,
          "a month is its middle")
    check(pd("2012-03-01")[1] == 2
          and abs(pd("2012-03-01")[0] - (2012 + 60.5 / 366)) < 1e-9,
          "a day, in a leap year")
    check(pd("2010-03-01T00:00:00Z") == pd("2010-03-01")
          and pd("2010/2011") == (2011.0, 0) and pd("2001/2003")[1] == -1,
          "timestamps, and ranges at their midpoint")
    check(pd("2010.25") == (2010.25, 2) and pd("2010-13") is None
          and pd("2010-02-30") is None and pd("soon") is None,
          "decimal years, and bad dates refused")
    lines = ["run\tsample\twhen\n", "r1\tA\t2001-05-02\n",
             "r2\tB\t2003\n", "r3\tC\tNA\n", "r4\tD\t2004-06\n",
             "r5\tD\t2005-06\n", "r6\tE\tlast spring\n",
             "r7\tF\t2002-02-02\n", "r8\tF\t2002-02-02\n"]
    got, st = read_dates(lines, "sample", "3", "month")
    check(sorted(got) == ["A", "F"] and st["imprecise"] == 1
          and st["empty"] == 1 and st["unparsed"] == 1
          and st["conflicting"] == 1,
          f"read_dates: columns by name and number, precision, conflicts "
          f"({sorted(got)}, {st})")

    ### A timed tree: lineages mutate at `rate` per year along branches of
    ### known duration.  hub0 holds a clade sampled 30 years early beside a
    ### normal sister; hub1 a mutator (1.8x) beside a normal sister.
    rng = np.random.default_rng(5)
    rate = 0.5
    t = _Tree()
    dates = {}
    root = t.add("node_1", -1, [])

    def branch(years, r):
        return list(rng.integers(0, 12, size=rng.poisson(r * years)))

    def clade(parent, start, r, lo, hi, name=None, levels=6, step=10):
        top = t.add(name or f"n{len(t.node_ids)}", parent,
                    branch(step, r))
        frontier = [top]
        for lev in range(1, levels):
            nxt = []
            for p_ in frontier:
                for _ in range(2):
                    nxt.append(t.add(f"n{len(t.node_ids)}", p_,
                                     branch(step, r)))
            frontier = nxt
        when = start + levels * step
        for p_ in frontier:
            for _ in range(2):
                d = rng.uniform(lo, hi)
                tip = t.add(f"s{len(t.node_ids)}", p_,
                            branch(d - when, r))
                dates[t.node_ids[tip]] = d
        return top

    hubs = [t.add(f"hub{h}", root, branch(100, rate)) for h in range(40)]
    clade(hubs[0], 100, rate, 165, 205, name="node_old")
    clade(hubs[0], 100, rate, 195, 235)
    clade(hubs[1], 100, 1.8 * rate, 195, 235, name="node_fast")
    clade(hubs[1], 100, rate, 195, 235)
    for h in hubs[2:]:
        clade(h, 100, rate, 175, 235, levels=7)

    args = Args()
    args.min_dated_tips = 10
    args.fold = False
    import io
    buf = io.StringIO()
    table, runs, ctl, _ = analyse(np, t, ["node_old", "node_fast"], args,
                                  buf, dates=dates)
    v = {x["node"]: x for x in table}
    o, f = v["node_old"], v["node_fast"]
    check(o["log2_rate"] < -0.15 and abs(o["adj_log2_rate"]) < 0.25,
          f"an early-sampled clade looks slow ({o['log2_rate']:+.2f}) "
          f"until adjusted ({o['adj_log2_rate']:+.2f})")
    check(o["date_verdict"] == "not different"
          and o["reading"] not in ("antimutator", "genome-wide slower"),
          f"... and is not called slower with dates ({o['date_verdict']}, "
          f"{o['reading']})")
    check(f["reading"] in ("genome-wide faster", "suggestive faster",
                           "terminal only"),
          f"reading: a mutator of every class alike is not class-specific "
          f"({f['reading']})")
    check(f["date_verdict"] == "faster" and f["adj_log2_rate"] > 0.5,
          f"the mutator stays faster with dates ({f['date_verdict']}, "
          f"{f['adj_log2_rate']:+.2f})")
    check(f["slope_S"] > f["slope_out"],
          f"the mutator's clock is steeper ({f['slope_S']:.2f} against "
          f"{f['slope_out']:.2f})")
    m = re.search(r"all branches: ([\d.]+) mutations per year", buf.getvalue())
    check(m is not None and abs(float(m.group(1)) - rate) < 0.15,
          f"the tree's clock is about {rate} per year "
          f"({m.group(1) if m else None})")


# --------------------------------------------------------------------------

def read_spectra_ids(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if "Node_ID" not in hdr:
            raise SystemExit("error: --spectra has no Node_ID column")
        ci = hdr.index("Node_ID")
        folded = any("|" in h for h in hdr)
        ids = [f[ci] for f in (line.rstrip("\n").split("\t") for line in fh)
               if len(f) > ci and f[ci]]
    return ids, folded


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    for opt in ("input_tree", "spectra"):
        if not getattr(args, opt):
            print(f"error: --{opt} is required (or pass --selftest)",
                  file=sys.stderr)
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
    from permutation_null import extract
    split_ids, folded = read_spectra_ids(args.spectra)
    if folded != args.fold:
        print(f"  WARNING: {args.spectra} looks like a search run "
              f"{'with' if folded else 'without'} --fold, but this run is "
              f"{'with' if args.fold else 'without'} it.", file=sys.stderr)
    print(f"Reading {args.input_tree} ...", file=sys.stderr)
    mat = bte.MATree(args.input_tree)
    tree = extract(mat.root)
    del mat
    comp = read_composition(args.composition) if args.composition else None
    dates = dates_note = st = None
    if args.dates:
        opener = gzip.open if args.dates.endswith(".gz") else open
        with opener(args.dates, "rt") as fh:
            dates, st = read_dates(fh, args.sample_column, args.date_column,
                                   args.min_date_precision)
        pr = st["precision"]
        dates_note = [
            f"{args.dates}: {st['rows']} rows, {len(dates)} dates kept "
            f"(day {pr[2]}, month {pr[1]}, year {pr[0]}); {st['empty']} "
            f"empty, {st['unparsed']} unparsed, {st['imprecise']} less "
            f"precise than {args.min_date_precision}, {st['conflicting']} "
            f"samples with conflicting dates"]
        if st["examples"]:
            dates_note.append(f"  unparsed, e.g. "
                              f"{', '.join(sorted(set(st['examples'])))}")
        for line in dates_note:
            print(line, file=sys.stderr)
    elif args.date_column:
        print("error: --date_column needs --dates", file=sys.stderr)
        return 2
    table, runs, controls, cats = analyse(np, tree, split_ids, args,
                                          sys.stdout, comp, dates,
                                          dates_note, st)
    paths = write_outputs(table, runs, controls, cats, args.prefix)
    print(f"Wrote {', '.join(paths)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
