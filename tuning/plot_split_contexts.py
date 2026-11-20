#!/usr/bin/env python3
"""Figures for split_contexts.py's results.

Reads the TSVs one split_contexts.py run writes (--prefix) and draws three
figures:

  PREFIX_overview   which splits stand out.  Excess -- G(within) over its
                    shuffled mean -- against the split's size, one panel per
                    run (all branches, internal branches only), over the
                    cloud of random non-split clades the ctrl_p null draws
                    from.  Splits in the summary (both nulls, both runs) are
                    filled and labelled; splits that pass the shuffle but
                    have too few clades of their size to be judged are
                    crosses; the rest are hollow.  The dashed line is
                    excess 1, the shuffle's expectation.

  PREFIX_profiles   what the summary splits are.  One row per split, the
                    contexts in signature order and grouped by class: each
                    bar is the change in the context's share of its class,
                    split minus comparator, on all branches -- the quantity
                    G(within) tests, so a class shifted in proportion draws
                    flat.  Dots are the same on internal branches only.  The
                    two contexts carrying most of G(within) are named, and a
                    triangle marks a context where one site carries
                    --one_site or more of the split's count (with at least
                    --site_min_count mutations): a hotspot, not a context.

  PREFIX_heatmap    patterns shared between splits.  Every split tested,
                    by all-branch excess, against every context.  The colour
                    is sign(share change) x sqrt(the context's G
                    contribution): a signed deviance, which weighs a shift by
                    the counts behind it, so small noisy splits do not
                    dominate.  Cells contributing less than --min_cell_g
                    (default 3.84, nominal p 0.05 for one degree of
                    freedom) are drawn neutral, so the colour is the cells
                    that clear noise.  Red: more of the class in the split than in
                    its comparator; blue: less.  Summary splits are in bold.

Without --fold the contexts are the 192 forward-strand classes; each folded
pair of types shares a colour, and the purine-reference type (G>A, A>G, ...)
is hatched.

Needs only matplotlib, so it runs outside the bte environment.

Usage:
    plot_split_contexts.py --prefix split_contexts
    plot_split_contexts.py --prefix split_contexts --format png --output figs/sc
    plot_split_contexts.py --selftest
"""

import argparse
import csv
import math
import os
import sys

### Reference palette (light mode): chart chrome and ink, six categorical
### slots validated as an adjacent set (one per folded class, in the
### classes' order), and the blue/red diverging pair with a gray midpoint.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CLASS_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                "#008300"]
BLUE, MID, RED = "#2a78d6", "#f0efec", "#e34948"
HIGHLIGHT = "#2a78d6"

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Plot split_contexts.py's results: an overview, context "
                    "profiles of the summary splits, and a heatmap.")
    p.add_argument("--prefix", default="split_contexts",
                   help="The --prefix split_contexts.py wrote with "
                        "(default split_contexts)")
    p.add_argument("--output", default=None,
                   help="Prefix for the figures (default: --prefix)")
    p.add_argument("--format", default="pdf",
                   help="pdf (default), png or svg")
    p.add_argument("--max_q", type=float, default=0.05,
                   help="Branch-shuffle q a split must be below to count as "
                        "passing the shuffle in the overview (default 0.05, "
                        "as in split_contexts.py)")
    p.add_argument("--one_site", type=float, default=0.2,
                   help="Top-site share marked as a hotspot in the profiles "
                        "(default 0.2)")
    p.add_argument("--site_min_count", type=int, default=20,
                   help="Smallest count in a context for its top-site share "
                        "to be marked (default 20)")
    p.add_argument("--min_cell_g", type=float, default=3.84,
                   help="Heatmap cells whose G contribution is below this "
                        "are drawn neutral (default 3.84, chi-squared 1 df "
                        "at p 0.05; nominal, as the mutations are not "
                        "independent)")
    p.add_argument("--selftest", action="store_true",
                   help="Draw the figures from constructed files and exit")
    return p.parse_args(argv)


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def context_file(prefix, tag=""):
    for n in (96, 192):
        path = f"{prefix}{tag}_{n}.tsv"
        if os.path.exists(path):
            return path, n
    raise SystemExit(f"error: no {prefix}{tag}_96.tsv or _192.tsv; is "
                     f"--prefix the one split_contexts.py wrote with?")


def read_contexts(path):
    """node -> list of context rows, in file (signature) order."""
    out = {}
    for r in read_tsv(path):
        out.setdefault(r["node"], []).append(r)
    return out


def load(prefix):
    data = {}
    data["all"] = read_tsv(f"{prefix}.tsv")
    data["internal"] = read_tsv(f"{prefix}_internal.tsv")
    path, n = context_file(prefix)
    data["n"] = n
    data["ctx_all"] = read_contexts(path)
    data["ctx_int"] = read_contexts(context_file(prefix, "_internal")[0])
    data["table"] = read_tsv(f"{prefix}_table.tsv")
    cpath = f"{prefix}_controls.tsv"
    data["controls"] = read_tsv(cpath) if os.path.exists(cpath) else None
    return data


# --------------------------------------------------------------------------
# context layout
# --------------------------------------------------------------------------

def class_style(ctx_rows):
    """Per context: (colour, hatched, class label) and the block boundaries.

    Folded files have 6 classes like CT|GA; forward-strand files 12 types
    like CT and GA, which share their folded pair's colour, the
    purine-reference one hatched.
    """
    labels = []
    for r in ctx_rows:
        if r["class"] not in labels:
            labels.append(r["class"])
    folded = all("|" in c for c in labels)
    fold_order = ["CA|GT", "CG|GC", "CT|GA", "TA|AT", "TC|AG", "TG|AC"]
    style = []
    for r in ctx_rows:
        c = r["class"]
        if folded:
            k, hatch = fold_order.index(c), False
        else:
            ref, alt = c[0], c[1]
            if ref in "GA":
                ref, alt = COMPLEMENT[ref], COMPLEMENT[alt]
                hatch = True
            else:
                hatch = False
            k = fold_order.index(f"{ref}{alt}|{COMPLEMENT[ref]}"
                                 f"{COMPLEMENT[alt]}")
        style.append((CLASS_COLORS[k], hatch, c))
    bounds = [i for i in range(1, len(ctx_rows))
              if ctx_rows[i]["class"] != ctx_rows[i - 1]["class"]]
    return style, labels, bounds


def trinuc(label):
    """A[C>T]G -> ACG, the trinucleotide the context sits in."""
    return label[0] + label[2] + label[6]


def share_diff(r):
    return float(r["split_class_share"]) - float(r["comparator_class_share"])


def signed_deviance(r, min_g=0.0):
    g = max(float(r["g_contribution"]), 0.0)
    if g < min_g:
        return 0.0
    d = share_diff(r)
    return math.copysign(math.sqrt(g), d) if d else 0.0


def quantile(vals, q):
    v = sorted(vals)
    if not v:
        return 0.0
    i = min(len(v) - 1, max(0, int(round(q * (len(v) - 1)))))
    return v[i]


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

def style_axes(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=7)


def label_points(fig, ax, points, fontsize=6.5):
    """Label points without overlap: for each, highest first, the first of
    a few offsets whose box clears the labels already placed, with a
    hairline leader when the label sits away from its point."""
    fig.canvas.draw()
    to_px = ax.transData.transform
    pt = fontsize * fig.dpi / 72.0
    placed = []
    offsets = [(5, 2), (5, -10), (-5, 2), (-5, -10), (8, 12), (8, -20),
               (-8, 12), (-8, -20), (12, 22), (12, -30), (-12, 22),
               (-12, -30)]
    for x, y, text in sorted(points, key=lambda p: -p[1]):
        px, py = to_px((x, y))
        w, h = 0.66 * pt * len(text), 1.7 * pt
        for dx, dy in offsets:
            sx, sy = dx * fig.dpi / 72.0, dy * fig.dpi / 72.0
            x0 = px + sx if dx > 0 else px + sx - w
            box = (x0, py + sy, x0 + w, py + sy + h)
            if all(box[2] < b[0] or box[0] > b[2] or box[3] < b[1]
                   or box[1] > b[3] for b in placed):
                break
        placed.append(box)
        far = abs(dx) > 5 or dy not in (2, -10)
        ax.annotate(text, (x, y), xytext=(dx, dy),
                    textcoords="offset points", fontsize=fontsize,
                    color=INK, ha="left" if dx > 0 else "right",
                    va="bottom",
                    arrowprops=dict(arrowstyle="-", color=AXIS,
                                    linewidth=0.5, shrinkA=0, shrinkB=3)
                    if far else None)


def overview(plt, data, args, path):
    summary = {r["split"] for r in data["table"]}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, run, title in ((axes[0], "all", "All branches"),
                           (axes[1], "internal", "Internal branches only")):
        style_axes(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        if data["controls"]:
            cs = [c for c in data["controls"] if c["branches"] == run]
            ax.scatter([float(c["mutations"]) for c in cs],
                       [float(c["excess"]) for c in cs], s=7, color=MUTED,
                       alpha=0.35, linewidths=0,
                       label="random non-split clade (control)")
        ax.axhline(1.0, color=INK2, linewidth=0.8, linestyle=(0, (4, 3)),
                   label="excess 1: the shuffle's expectation")
        rows = [r for r in data[run] if num(r["excess"]) is not None
                and num(r["split_mutations"])]
        other, stuck, top = [], [], []
        for r in rows:
            if r["node"] in summary:
                top.append(r)
            elif (num(r["q"]) is not None and num(r["q"]) < args.max_q
                  and r.get("control_reachable") == "no"):
                stuck.append(r)
            else:
                other.append(r)

        def xy(rs):
            return ([float(r["split_mutations"]) for r in rs],
                    [float(r["excess"]) for r in rs])
        ax.scatter(*xy(other), s=26, facecolors="none", edgecolors=INK2,
                   linewidths=0.8, label="split")
        if stuck:
            ax.scatter(*xy(stuck), s=30, marker="x", color=INK2,
                       linewidths=1.0,
                       label="split passing the shuffle, too few clades "
                             "of its size to judge")
        ax.scatter(*xy(top), s=46, color=HIGHLIGHT, edgecolors=SURFACE,
                   linewidths=1.0, zorder=3,
                   label="passes both nulls on both runs")
        ax.set_xscale("log")
        label_points(fig, ax, [(float(r["split_mutations"]),
                                float(r["excess"]), r["node"]) for r in top])
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("called mutations in the split's residual",
                      fontsize=8, color=INK2)
    axes[0].set_ylabel("excess: G(within) / shuffled mean", fontsize=8,
                       color=INK2)
    if not data["controls"]:
        axes[0].text(0.02, 0.98, f"no {args.prefix}_controls.tsv: rerun "
                     "split_contexts.py to draw the controls",
                     transform=axes[0].transAxes, fontsize=7, color=INK2,
                     va="top")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7,
               frameon=False, labelcolor=INK2)
    fig.suptitle("Trinucleotide context within classes: which splits "
                 "differ more than any clade of their size", fontsize=11,
                 color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def profiles(plt, data, args, path):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    rows = data["table"]
    if not rows:
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.axis("off")
        ax.text(0.5, 0.5, "No split passes both nulls on both runs.",
                ha="center", va="center", color=INK2)
        fig.savefig(path, facecolor=SURFACE)
        plt.close(fig)
        return
    n = data["n"]
    fig, axes = plt.subplots(len(rows), 1,
                             figsize=(12, 1.55 * len(rows) + 1.3),
                             sharex=True, squeeze=False)
    axes = [a[0] for a in axes]
    fig.patch.set_facecolor(SURFACE)
    first = data["ctx_all"][rows[0]["split"]]
    style, labels, bounds = class_style(first)
    xs = list(range(n))
    for ax, t in zip(axes, rows):
        node = t["split"]
        ca = data["ctx_all"][node]
        ci = {r["context"]: r for r in data["ctx_int"].get(node, [])}
        style_axes(ax)
        d_all = [share_diff(r) for r in ca]
        d_int = [share_diff(ci[r["context"]]) if r["context"] in ci else None
                 for r in ca]
        for x, d, (col, hatch, _c) in zip(xs, d_all, style):
            ax.bar(x, d, width=0.78, color=col, linewidth=0,
                   hatch="////" if hatch else None,
                   edgecolor=SURFACE if hatch else col)
        ax.scatter([x for x, d in zip(xs, d_int) if d is not None],
                   [d for d in d_int if d is not None], s=5, color=INK,
                   zorder=3, linewidths=0)
        ax.axhline(0, color=AXIS, linewidth=0.6)
        for b in bounds:
            ax.axvline(b - 0.5, color=GRID, linewidth=0.6)
        lim = max([abs(d) for d in d_all]
                  + [abs(d) for d in d_int if d is not None] + [1e-3]) * 1.3
        ax.set_ylim(-lim, lim)
        ax.set_xlim(-0.8, n - 0.2)
        ### The contexts the summary names, and one-site hotspots.
        named = [c.split(" ")[0] for c in t["contexts"].split("; ") if c]
        last = {}
        for x, r, d in zip(xs, ca, d_all):
            if r["context"] in named:
                up = d >= 0
                ### A neighbour's label on the same side: lift this one.
                lift = 9 if up in last and x - last[up] <= 6 else 0
                last[up] = x
                ax.annotate(r["context"], (x, d),
                            xytext=(0, (3 + lift) if up else -(3 + lift)),
                            textcoords="offset points", ha="center",
                            va="bottom" if up else "top", fontsize=6.5,
                            color=INK,
                            arrowprops=dict(arrowstyle="-", color=AXIS,
                                            linewidth=0.5, shrinkA=0,
                                            shrinkB=0) if lift else None)
            share = num(r.get("split_top_site_share"))
            if share is not None and share >= args.one_site \
                    and int(r["split_count"]) >= args.site_min_count:
                ax.scatter([x], [lim * 0.88], marker="v", s=18, color=INK,
                           zorder=4)
        ax.set_title(f"{node}   {t['clade'] or ''}   excess "
                     f"{t['excess_all']} / {t['excess_internal']}   ctrl_q "
                     f"{t['ctrl_q_all']} / {t['ctrl_q_internal']}",
                     fontsize=7.5, color=INK, loc="left", pad=3)
        ax.tick_params(axis="y", labelsize=6)
    ### Class names over the first row's blocks: identity is never colour
    ### alone.
    starts = [0] + bounds
    ends = bounds + [n]
    for lab, s0, e0 in zip(labels, starts, ends):
        col, hatch, _c = style[s0]
        axes[0].text((s0 + e0 - 1) / 2, 1.32, lab.replace("|", " | "),
                     transform=axes[0].get_xaxis_transform(), ha="center",
                     fontsize=7.5, color=INK)
        axes[0].add_patch(plt.Rectangle(
            (s0 - 0.5, 1.2), e0 - s0, 0.07,
            transform=axes[0].get_xaxis_transform(), clip_on=False,
            facecolor=col, hatch="////" if hatch else None,
            edgecolor=SURFACE if hatch else col, linewidth=0))
    axes[-1].set_xticks(xs)
    axes[-1].set_xticklabels([trinuc(r["context"]) for r in first],
                             rotation=90, fontsize=4.5, family="monospace",
                             color=INK2)
    fig.text(0.005, 0.5, "change in the context's share of its class, "
             "split minus comparator", rotation=90, va="center",
             fontsize=7.5, color=INK2)
    legend = [Patch(color=AXIS, label="bars, coloured by class: all branches"),
              Line2D([], [], marker="o", linestyle="", color=INK,
                     markersize=3, label="dots: internal branches only"),
              Line2D([], [], marker="v", linestyle="", color=INK,
                     markersize=4,
                     label=f"one site carries >= {args.one_site:g} of the "
                           f"split's count")]
    fig.legend(handles=legend, loc="lower center", ncol=3, fontsize=7,
               frameon=False, labelcolor=INK2)
    fig.suptitle("Where in the class: context profiles of the splits that "
                 "pass both nulls on both runs", fontsize=11, color=INK,
                 x=0.01, ha="left")
    fig.tight_layout(rect=(0.015, 0.04, 1, 0.94), h_pad=1.4)
    height = fig.get_size_inches()[1]
    fig.subplots_adjust(top=1 - 1.2 / height)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def heatmap(plt, data, args, path):
    from matplotlib.colors import LinearSegmentedColormap
    summary = {r["split"] for r in data["table"]}
    rows = sorted((r for r in data["all"] if num(r["excess"]) is not None
                   and r["node"] in data["ctx_all"]),
                  key=lambda r: -float(r["excess"]))
    if not rows:
        return
    n = data["n"]
    first = data["ctx_all"][rows[0]["node"]]
    style, labels, bounds = class_style(first)
    mat = [[signed_deviance(c, args.min_cell_g)
            for c in data["ctx_all"][r["node"]]] for r in rows]
    vmax = max(quantile([abs(v) for row in mat for v in row if v], 0.99),
               math.sqrt(args.min_cell_g) * 2)
    cmap = LinearSegmentedColormap.from_list("div", [BLUE, MID, RED])
    h = 0.15 * len(rows) + 2.2
    fig = plt.figure(figsize=(12, h))
    fig.patch.set_facecolor(SURFACE)
    top = 1 - 0.9 / h
    ax = fig.add_axes([0.09, 0.9 / h, 0.8, top - 0.9 / h - 0.25 / h])
    strip = fig.add_axes([0.09, top - 0.18 / h, 0.8, 0.12 / h])
    cax = fig.add_axes([0.905, 0.9 / h, 0.012, 0.3])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax,
                   interpolation="nearest")
    for b in bounds:
        ax.axvline(b - 0.5, color=SURFACE, linewidth=1.2)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["node"] for r in rows], fontsize=5.5)
    for lab, r in zip(ax.get_yticklabels(), rows):
        if r["node"] in summary:
            lab.set_fontweight("bold")
            lab.set_color(INK)
        else:
            lab.set_color(INK2)
    ax.set_xticks(range(n))
    ax.set_xticklabels([trinuc(c["context"]) for c in first], rotation=90,
                       fontsize=4.5, family="monospace", color=INK2)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ### Class strip with names: identity is never colour alone.
    starts = [0] + bounds
    ends = bounds + [n]
    strip.set_xlim(-0.5, n - 0.5)
    strip.set_ylim(0, 1)
    strip.axis("off")
    for lab, s0, e0 in zip(labels, starts, ends):
        col, hatch, _c = style[s0]
        strip.add_patch(plt.Rectangle((s0 - 0.5, 0), e0 - s0, 1, facecolor=col,
                                      hatch="////" if hatch else None,
                                      edgecolor=SURFACE if hatch else col,
                                      linewidth=0))
        strip.text((s0 + e0 - 1) / 2, 1.4, lab.replace("|", " | "),
                   ha="center", va="bottom", fontsize=7.5, color=INK)
    cb = fig.colorbar(im, cax=cax)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=6, colors=INK2)
    cb.set_label("sign of share change x sqrt(G contribution)\n"
                 "red: more of its class in the split; neutral: G "
                 f"contribution < {args.min_cell_g:g}", fontsize=6.5,
                 color=INK2)
    fig.text(0.01, 1 - 0.3 / h, "Context shifts across all splits, by "
             "all-branch excess (bold: passes both nulls on both runs)",
             fontsize=11, color=INK, va="top")
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def _write(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _fixture(prefix, n=96):
    import random
    rng = random.Random(3)
    bases = "ACGT"
    if n == 96:
        classes = [("CA|GT", "C", "A"), ("CG|GC", "C", "G"),
                   ("CT|GA", "C", "T"), ("TA|AT", "T", "A"),
                   ("TC|AG", "T", "C"), ("TG|AC", "T", "G")]
    else:
        classes = [(r + a, r, a) for r, a in
                   ("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT",
                    "TA", "TC", "TG")]
    nodes = [f"node_{i}" for i in range(1, 7)]
    split_cols = ["node", "comparator", "split_mutations", "excess", "q",
                  "control_reachable"]
    for tag in ("", "_internal"):
        _write(f"{prefix}{tag}.tsv", split_cols, [
            {"node": nd, "comparator": "node_0",
             "split_mutations": 10000 * (i + 1), "excess": 5 - 0.6 * i,
             "q": 0.002 if i < 4 else 0.3,
             "control_reachable": "no" if i == 3 else "yes"}
            for i, nd in enumerate(nodes)])
        rows = []
        for nd in nodes:
            for cl, r, a in classes:
                for left in bases:
                    for right in bases:
                        s = rng.random() / 8
                        rows.append({
                            "node": nd, "comparator": "node_0", "class": cl,
                            "context": f"{left}[{r}>{a}]{right}",
                            "split_count": rng.randint(0, 200),
                            "comparator_count": 100,
                            "split_class_share": s,
                            "comparator_class_share": 1 / 16,
                            "g_contribution": rng.random() * 9,
                            "split_top_site_share": 0.5 if s > 0.124
                            else 0.01})
        _write(f"{prefix}{tag}_{n}.tsv",
               ["node", "comparator", "class", "context", "split_count",
                "comparator_count", "split_class_share",
                "comparator_class_share", "g_contribution",
                "split_top_site_share"], rows)
    lab = f"A[{classes[0][1]}>{classes[0][2]}]A"
    _write(f"{prefix}_table.tsv",
           ["split", "clade", "excess_all", "excess_internal", "ctrl_q_all",
            "ctrl_q_internal", "contexts"],
           [{"split": nd, "clade": "L1 10", "excess_all": "5.0",
             "excess_internal": "3.0", "ctrl_q_all": "0.02",
             "ctrl_q_internal": "0.03",
             "contexts": f"{lab} 0.100 vs 0.050"} for nd in nodes[:2]])
    _write(f"{prefix}_controls.tsv",
           ["branches", "node", "group", "mutations", "excess", "p"],
           [{"branches": b, "node": f"c{i}", "group": "node_0",
             "mutations": rng.randint(5000, 90000),
             "excess": rng.random() * 2, "p": rng.random()}
            for b in ("all", "internal") for i in range(50)])


def selftest():
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    assert trinuc("A[C>T]G") == "ACG"
    r = {"split_class_share": "0.1", "comparator_class_share": "0.3",
         "g_contribution": "4"}
    assert abs(share_diff(r) + 0.2) < 1e-12 and signed_deviance(r) == -2.0
    assert signed_deviance(r, 4.5) == 0.0 and signed_deviance(r, 3.84) == -2.0
    rows = [{"class": "CT"}, {"class": "GA"}, {"class": "GA"}]
    style, labels, bounds = class_style(rows)
    assert labels == ["CT", "GA"] and bounds == [1]
    assert style[0][0] == style[1][0] == CLASS_COLORS[2]
    assert style[0][1] is False and style[1][1] is True
    print("ok: share change, signed deviance, and colours: a folded pair "
          "shares one, the purine-reference type hatched")

    with tempfile.TemporaryDirectory() as d:
        for n in (96, 192):
            prefix = os.path.join(d, f"sc{n}")
            _fixture(prefix, n)
            args = parse_args(["--prefix", prefix])
            data = load(prefix)
            assert data["n"] == n
            for name, fn in (("overview", overview), ("profiles", profiles),
                             ("heatmap", heatmap)):
                path = f"{prefix}_{name}.pdf"
                fn(plt, data, args, path)
                assert os.path.getsize(path) > 1000, (n, name)
            ### Without controls, the overview still draws and says so.
            os.remove(f"{prefix}_controls.tsv")
            overview(plt, load(prefix), args, f"{prefix}_nc.pdf")
    print("ok: all three figures draw from constructed 96- and 192-class "
          "outputs, and the overview draws without controls")
    print("all self-tests passed")
    return 0


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("error: matplotlib is required", file=sys.stderr)
        return 1
    data = load(args.prefix)
    out = args.output or args.prefix
    written = []
    for name, fn in (("overview", overview), ("profiles", profiles),
                     ("heatmap", heatmap)):
        path = f"{out}_{name}.{args.format}"
        fn(plt, data, args, path)
        written.append(path)
    if not data["controls"]:
        print(f"note: no {args.prefix}_controls.tsv, so the overview has no "
              f"control clades; rerun split_contexts.py to get it",
              file=sys.stderr)
    print(f"Wrote {', '.join(written)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
