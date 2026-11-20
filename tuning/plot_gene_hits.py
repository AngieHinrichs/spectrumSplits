#!/usr/bin/env python3
"""Figure for step 12: one gene's spectrum shift, against every gene and its null.

Two panels:

  A  Every gene in the list, both effect sets: the shift in one class
     (x, percentage points, hit child minus siblings) against the omnibus
     -log10 p (y).  Nonsynonymous hits in colour, the synonymous control in
     grey.  A dashed line marks 0.05 over the number of genes, a dotted one
     the smallest p the permutations can give.  The control cloud should sit
     below the line; a gene beyond it is the result.

  B  The chosen gene against its own permutation null, one row per test:
     nonsynonymous hits, nonsynonymous internal hits only (if given), and the
     synonymous control.  Histogram of the permuted shifts, with the observed
     shift marked.

Needs only matplotlib.  The null draws come from gene_hit_contrasts.py
--null_draws, run with the same settings and seed as the summary it goes
with, so the histograms are the nulls the p-values were computed from:

  gene_hit_contrasts.py ... --permutations 10000 \\
      --null_draws null_draws.tsv --null_draws_genes dnaE2 \\
      --output_events hit_events.tsv --output_summary hit_summary.tsv
  gene_hit_contrasts.py ... --permutations 10000 --top_branch exclude \\
      --null_draws null_draws_internal.tsv --null_draws_genes dnaE2 \\
      --output_events hit_events_internal.tsv \\
      --output_summary hit_summary_internal.tsv

  plot_gene_hits.py --summary hit_summary.tsv --null null_draws.tsv \\
      --internal_summary hit_summary_internal.tsv \\
      --internal_null null_draws_internal.tsv \\
      --gene dnaE2 --class 'CT|GA' --label udgB,ogt,lexA \\
      --output dnaE2_figure.png

The output format follows the extension (.png, .pdf, .svg).
"""

import argparse
import csv
import math
import sys

### Reference palette (light mode): one hue for the result, neutral grey for
### the control, text in ink rather than series colour.
BLUE = "#2a78d6"
GREY = "#9b9a95"
NULL_FILL = "#d9d8d3"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a8984"
SURFACE = "#ffffff"
EFFECTS = (("nonsynonymous", "nonsynonymous hits", BLUE),
           ("synonymous", "synonymous hits (control)", GREY))


def read_summary(path):
    with open(path) as fh:
        return {(r["effect_set"], r["gene"]): r
                for r in csv.DictReader(fh, delimiter="\t")}


def read_null(path, gene, cls):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] == gene:
                out.setdefault(r["effect_set"], []).append(
                    float(r[f"diff_{cls}"]))
    return out


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def class_label(cls):
    """CT|GA -> C>T / G>A."""
    return " / ".join(f"{p[0]}>{p[1]}" for p in cls.split("|"))


def fmt_p(p):
    return f"{p:.1g}" if p < 0.001 else f"{p:.2g}"


def label_points(ax, points, gene):
    """Name points to their right, nudging a label down past any it would
    overlap (a greedy pass in screen coordinates, top to bottom)."""
    fig = ax.figure
    fig.canvas.draw()
    to_px = ax.transData.transform
    boxes = []
    for x, y, g in sorted(points, key=lambda p: -p[1]):
        px, py = to_px((x, y))
        dy = 0
        while any(abs(px - bx) < 60 and abs(py + dy - by) < 11
                  for bx, by in boxes):
            dy -= 11
        boxes.append((px, py + dy))
        ax.annotate(g, (x, y), xytext=(6, dy * 72 / fig.dpi),
                    textcoords="offset points", va="center", fontsize=9,
                    color=INK if g == gene else INK_2,
                    fontweight="bold" if g == gene else "normal",
                    fontstyle="italic",
                    bbox=dict(boxstyle="square,pad=0.1", fc=SURFACE,
                              ec="none"),
                    arrowprops=None if dy == 0 else dict(
                        arrowstyle="-", color=INK_3, lw=0.6,
                        shrinkA=0, shrinkB=3))


def panel_a(ax, summary, cls, gene, labels, n_perm):
    genes = sorted({g for (_e, g) in summary if g != "ALL"})
    n_genes = len(genes)
    handles = {}
    placed = []
    for es, name, colour in reversed(EFFECTS):
        xs, ys, names = [], [], []
        for g in genes:
            r = summary.get((es, g))
            if not r:
                continue
            d, p = num(r.get(f"diff_{cls}")), num(r.get("omnibus_p"))
            if d is None or p is None:
                continue
            xs.append(100 * d)
            ys.append(-math.log10(p))
            names.append(g)
        handles[es] = ax.scatter(xs, ys, s=38, color=colour,
                                 edgecolors=SURFACE, linewidths=0.8,
                                 zorder=3, label=name)
        if es == "nonsynonymous":
            placed = [(x, y, g) for x, y, g in zip(xs, ys, names)
                      if g == gene or g in labels]
    bar = -math.log10(0.05 / n_genes)
    lines = [(ax.axhline(bar, color=INK_2, lw=1, ls=(0, (4, 3)), zorder=1),
              f"p = 0.05 / {n_genes} genes")]
    top = bar
    if n_perm:
        floor = math.log10(n_perm + 1)
        lines.append((ax.axhline(floor, color=INK_3, lw=0.8,
                                 ls=(0, (1, 2)), zorder=1),
                      f"smallest p with {n_perm:,} permutations"))
        top = max(top, floor)
    ax.set_ylim(0, max(ax.get_ylim()[1], top) + 0.35)
    label_points(ax, placed, gene)
    ax.axvline(0, color=INK_3, lw=0.8, zorder=1)
    ax.set_xlabel(f"{class_label(cls)} shift below the hit "
                  f"(percentage points)")
    ax.set_ylabel("omnibus  −log$_{10}$ p")
    series = [(handles[e], n) for e, n, _c in EFFECTS if e in handles]
    ### Outside the plot, so nothing in it competes with the points.
    ax.legend([h for h, _n in series + lines], [n for _h, n in series + lines],
              loc="upper left", bbox_to_anchor=(0, -0.13), ncol=2,
              frameon=False, fontsize=8.5, columnspacing=1.5,
              handlelength=2.2)
    ax.set_title(f"A   All {n_genes} genes", loc="left", fontsize=10.5,
                 color=INK, fontweight="bold", pad=16)


def panel_b(axes, rows, cls):
    lo = min(min(min(r[1]), r[2]) for r in rows)
    hi = max(max(max(r[1]), r[2]) for r in rows)
    pad = 0.08 * (hi - lo)
    ### The observed value and its p sit to the right of the line, clear of
    ### the histogram, so the axis runs on far enough to hold them.
    right = max(hi + pad, max(r[2] for r in rows) + 0.45 * (hi - lo))
    edges_n = 60
    width = (hi - lo + 2 * pad) / edges_n
    bins = [lo - pad + i * width for i in range(edges_n + 1)]
    for ax, (label, null, obs, p, colour, k, eff) in zip(axes, rows):
        ax.hist(null, bins=bins, color=NULL_FILL, edgecolor=SURFACE,
                linewidth=0.5)
        ax.axvline(obs, color=colour, lw=2)
        ax.annotate(f"observed {obs:+.2f}\np = {fmt_p(p)}", (obs, 1.0),
                    xycoords=("data", "axes fraction"), xytext=(6, -2),
                    textcoords="offset points", color=INK, fontsize=8.5,
                    va="top", ha="left")
        ax.set_xlim(lo - pad, right)
        events = f"{k:,} events" + (f", {eff:,.0f} effective" if eff else "")
        ax.text(0.0, 1.02, f"{label}  ({events})", transform=ax.transAxes,
                fontsize=8.5, color=INK_2, va="bottom")
        ax.set_yticks([])
        for side in ("left", "right", "top"):
            ax.spines[side].set_visible(False)
    axes[-1].set_xlabel(f"{class_label(cls)} shift (percentage points)")


def main():
    p = argparse.ArgumentParser(
        description="Plot one gene's spectrum shift against all genes and "
                    "its own permutation null.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", required=True,
                   help="hit_summary.tsv from gene_hit_contrasts.py")
    p.add_argument("--null", required=True,
                   help="--null_draws TSV from the same run")
    p.add_argument("--internal_summary", default=None,
                   help="hit_summary.tsv from a --top_branch exclude run")
    p.add_argument("--internal_null", default=None,
                   help="--null_draws TSV from that run")
    p.add_argument("--gene", required=True)
    p.add_argument("--class", dest="cls", required=True,
                   help="Class to plot, as in the summary columns, e.g. "
                        "'CT|GA'")
    p.add_argument("--label", default="",
                   help="Comma-separated other genes to name in panel A")
    p.add_argument("--title", default=None,
                   help="Figure title (default names the gene and class)")
    p.add_argument("--output", default="gene_hits.png")
    args = p.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("error: matplotlib is needed; run this in an environment "
              "that has it", file=sys.stderr)
        return 1
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.edgecolor": INK_3, "axes.labelcolor": INK_2,
        "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
    })

    summary = read_summary(args.summary)
    nulls = read_null(args.null, args.gene, args.cls)
    if not nulls:
        print(f"error: no draws for {args.gene} in {args.null}; was it in "
              f"--null_draws_genes?", file=sys.stderr)
        return 1
    n_perm = len(nulls.get("nonsynonymous", []))

    def row(summ, nl, es, label, colour):
        r = summ.get((es, args.gene))
        if not r or es not in nl:
            return None
        return (label, [100 * x for x in nl[es]],
                100 * float(r[f"diff_{args.cls}"]),
                float(r[f"p_{args.cls}"]), colour, int(r["informative"]),
                num(r.get("effective_events")))

    rows = [row(summary, nulls, "nonsynonymous", "nonsynonymous, all hits",
                BLUE)]
    if args.internal_summary and args.internal_null:
        rows.append(row(read_summary(args.internal_summary),
                        read_null(args.internal_null, args.gene, args.cls),
                        "nonsynonymous", "nonsynonymous, internal hits only",
                        BLUE))
    rows.append(row(summary, nulls, "synonymous",
                    "synonymous, all hits (control)", GREY))
    rows = [r for r in rows if r]

    fig = plt.figure(figsize=(10, 6.4))
    grid = fig.add_gridspec(len(rows), 2, width_ratios=[1.15, 1],
                            wspace=0.22, hspace=0.95, top=0.85, bottom=0.2)
    ax_a = fig.add_subplot(grid[:, 0])
    labels = {g for g in args.label.split(",") if g}
    panel_a(ax_a, summary, args.cls, args.gene, labels, n_perm)
    axes_b = [fig.add_subplot(grid[i, 1]) for i in range(len(rows))]
    for ax in axes_b[1:]:
        ax.sharex(axes_b[0])
    panel_b(axes_b, rows, args.cls)
    axes_b[0].set_title(f"B   {args.gene} against its permutation null "
                        f"(grey)",
                        loc="left", fontsize=10.5, color=INK,
                        fontweight="bold", pad=16)
    title = args.title or (f"Coding changes in {args.gene} are followed by "
                           f"more {class_label(args.cls)} in the subtree "
                           f"below, relative to sibling subtrees")
    fig.suptitle(title, x=0.06, ha="left", fontsize=11.5, color=INK, y=0.97)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
