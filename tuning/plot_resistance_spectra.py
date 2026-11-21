#!/usr/bin/env python3
"""Figures for resistance_spectra.py's results.

Reads the TSVs one resistance_spectra.py run writes (--prefix) and draws
three figures:

  PREFIX_rate         does a split's spectrum change how much resistance it
                      acquires?  One row per split, on one log axis of
                      events relative to the tree-wide expectation
                      (E_global): a hollow circle for what the split's
                      spectrum alone predicts (spectrum_factor), a filled
                      one for what it got (events / E_global), joined.
                      Coloured by the split's main lineage (--composition).
                      Rows are sorted by the all-branch observed rate, in the
                      same order on the internal-branch panel.  A split with
                      no events is a left-pointing triangle at the axis's
                      left edge.

  PREFIX_composition  does the spectrum tilt which variants arise?  Each
                      tested split's llr against llr_expected, what its
                      spectrum alone would give, one panel per reading.  On
                      the dashed diagonal a split leans its spectrum's way
                      by as much as the spectrum predicts; above it, further
                      than the spectrum explains; below zero, against it.
                      Splits at llr q below --max_q are filled and
                      labelled.  Both axes are symmetric-log: linear within
                      1 of zero, logarithmic beyond.

  PREFIX_example      what an llr is made of, for one split and drug
                      (--example_drug, --example_split, --example_reading):
                      for its --top_variants commonest variants, the rest of
                      the tree's share (q), the share the split's spectrum
                      predicts (q'), and the split's own share.  llr sums
                      log(q'/q) over the split's events, so it is positive
                      where the split's own bars follow q' away from q.

--drug picks the drug of the first two figures (default 'any drug', every
drug pooled).

Needs only matplotlib, so it runs outside the bte environment.

Usage:
    plot_resistance_spectra.py --prefix resistance_spectra \\
        --composition split_composition.tsv
    plot_resistance_spectra.py --selftest
"""

import argparse
import csv
import math
import os
import sys
from collections import Counter

### Reference palette (light mode), as in plot_split_contexts.py.  Five
### categorical slots for the lineages, validated as an adjacent set, gray
### for 'other'.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
HIGHLIGHT = "#2a78d6"
LINEAGE_GROUPS = ["L1", "L2", "L3", "L4", "M. africanum (L5, L6, La)"]
### q, q' and observed in the example.
REST, PREDICTED, OBSERVED = MUTED, "#eb6834", "#2a78d6"
ANY = "any drug"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Plot resistance_spectra.py's results: rate, "
                    "composition, and one worked example.")
    p.add_argument("--prefix", default="resistance_spectra",
                   help="The --prefix resistance_spectra.py wrote with "
                        "(default resistance_spectra)")
    p.add_argument("--composition", default=None,
                   help="split_composition.py --output, to colour splits by "
                        "lineage and name it")
    p.add_argument("--drug", default=ANY,
                   help="Drug of the rate and composition figures "
                        "(default 'any drug': all pooled)")
    p.add_argument("--max_q", type=float, default=0.05,
                   help="llr q below which a split is filled and labelled "
                        "(default 0.05)")
    p.add_argument("--example_drug", default="Streptomycin")
    p.add_argument("--example_split", default="node_14336")
    p.add_argument("--example_reading", choices=["all", "internal"],
                   default="internal")
    p.add_argument("--top_variants", type=int, default=8,
                   help="Variants in the example (default 8)")
    p.add_argument("--format", default="pdf",
                   help="Output format: pdf, png or svg (default pdf)")
    p.add_argument("--output", default=None,
                   help="Prefix of the figure files (default --prefix)")
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


def lineage_of(composition):
    """The main entry of split_composition.py's composition column
    ('lineage2.2.1:10684,lineage2.2.1.1:2938') shortened ('L2.2.1'), and
    its colour group."""
    best, top = None, -1
    for part in (composition or "").split(","):
        name, _, k = part.rpartition(":")
        try:
            k = int(k)
        except ValueError:
            continue
        if k > top and name:
            best, top = name, k
    if not best or best == "(none)":
        return "", "other"
    short = best.replace("lineage", "L")
    head = short.split(".")[0]
    if head in ("L1", "L2", "L3", "L4"):
        return short, head
    if head in ("L5", "L6") or head.startswith("La"):
        return short, LINEAGE_GROUPS[4]
    return short, "other"


def colour(group):
    return SLOTS[LINEAGE_GROUPS.index(group)] if group in LINEAGE_GROUPS \
        else MUTED


def load(prefix, composition=None):
    data = {"rows": read_tsv(f"{prefix}.tsv"), "variants": [],
            "lineage": {}}
    vpath = f"{prefix}_variants.tsv"
    if os.path.exists(vpath):
        data["variants"] = read_tsv(vpath)
    if composition:
        for r in read_tsv(composition):
            data["lineage"][r["node"]] = lineage_of(r.get("composition"))
    return data


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


def label_points(fig, ax, points, obstacles=(), fontsize=6.5):
    """Label points without overlap, as in plot_split_contexts.py, and also
    clear of every point in `obstacles` and inside the axes."""
    fig.canvas.draw()
    to_px = ax.transData.transform
    pt = fontsize * fig.dpi / 72.0
    r = 4 * fig.dpi / 72.0
    placed = [(px - r, py - r, px + r, py + r)
              for px, py in (to_px(o) for o in obstacles)]
    frame = ax.get_window_extent()
    offsets = [(5, 2), (5, -10), (-5, 2), (-5, -10), (8, 12), (8, -20),
               (-8, 12), (-8, -20), (12, 22), (12, -30), (-12, 22),
               (-12, -30), (20, 34), (-20, 34), (20, -42), (-20, -42)]
    for x, y, text in sorted(points, key=lambda p: -p[1]):
        px, py = to_px((x, y))
        w, h = 0.66 * pt * len(text), 1.6 * pt
        for dx, dy in offsets:
            sx, sy = dx * fig.dpi / 72.0, dy * fig.dpi / 72.0
            x0 = px + sx if dx > 0 else px + sx - w
            box = (x0, py + sy, x0 + w, py + sy + h)
            inside = (box[0] >= frame.x0 and box[2] <= frame.x1
                      and box[1] >= frame.y0 and box[3] <= frame.y1)
            if inside and all(box[2] < b[0] or box[0] > b[2]
                              or box[3] < b[1] or box[1] > b[3]
                              for b in placed):
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


def drug_rows(data, drug, reading):
    return {r["split"]: r for r in data["rows"]
            if r["drug"] == drug and r["reading"] == reading}


def rate(plt, data, args, path):
    from matplotlib.lines import Line2D
    by = {rd: drug_rows(data, args.drug, rd) for rd in ("all", "internal")}
    splits = list(by["all"]) or list(by["internal"])

    def observed(r):
        e = num(r.get("e_global")) if r else None
        return num(r["events"]) / e if e else None

    order = sorted(splits, key=lambda s: (observed(by["all"].get(s)) or 0.0))
    vals = [v for rd in by for s in splits for v in (
        observed(by[rd].get(s)), num((by[rd].get(s) or {}).get(
            "spectrum_factor"))) if v]
    lo = min(vals) / 2 if vals else 0.1
    hi = max(vals) * 1.5 if vals else 10
    n = len(order)
    fig, axes = plt.subplots(1, 2, figsize=(10, 1.6 + 0.13 * n),
                             sharey=True)
    fig.patch.set_facecolor(SURFACE)
    groups_seen = set()
    for ax, rd, title in ((axes[0], "all", "All branches"),
                          (axes[1], "internal", "Internal branches only")):
        style_axes(ax)
        ax.grid(axis="x", color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.axvline(1.0, color=INK2, linewidth=0.8, linestyle=(0, (4, 3)))
        for y, s in enumerate(order):
            r = by[rd].get(s)
            if not r:
                continue
            group = data["lineage"].get(s, ("", "other"))[1]
            groups_seen.add(group)
            c = colour(group)
            f = num(r["spectrum_factor"])
            o = observed(r)
            if f and o:
                ax.plot([f, o], [y, y], color=AXIS, linewidth=1.0, zorder=1)
            if f:
                ax.scatter([f], [y], s=18, facecolors=SURFACE,
                           edgecolors=INK2, linewidths=0.8, zorder=2)
            if o:
                ax.scatter([o], [y], s=22, color=c, edgecolors=SURFACE,
                           linewidths=0.6, zorder=3)
            elif o is not None or num(r["events"]) == 0:
                ax.scatter([lo * 1.05], [y], s=22, marker="<", color=c,
                           zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(lo, hi)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("resistance events / tree-wide expectation (E_global)",
                      fontsize=8, color=INK2)
    axes[0].set_yticks(range(n))
    axes[0].set_yticklabels(
        [f"{s}  {data['lineage'].get(s, ('', ''))[0]}".rstrip()
         for s in order], fontsize=5.5, color=INK2)
    axes[0].set_ylim(-0.8, n - 0.2)
    handles = [
        Line2D([], [], marker="o", linestyle="", markerfacecolor=SURFACE,
               markeredgecolor=INK2, markersize=5,
               label="predicted by the split's spectrum alone"),
        Line2D([], [], marker="o", linestyle="", color=INK2, markersize=5,
               label="observed"),
        Line2D([], [], color=INK2, linewidth=0.8, linestyle=(0, (4, 3)),
               label="1: the tree-wide expectation")]
    for g in LINEAGE_GROUPS + ["other"]:
        if g in groups_seen and data["lineage"]:
            handles.append(Line2D([], [], marker="o", linestyle="",
                                  color=colour(g), markersize=5, label=g))
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7,
               frameon=False, labelcolor=INK2)
    drug = "all drugs pooled" if args.drug == ANY else args.drug
    fig.suptitle(f"Resistance supply: what each split's spectrum predicts, "
                 f"and what it got ({drug})", fontsize=11, color=INK,
                 x=0.01, ha="left",
                 y=1 - 0.08 / fig.get_size_inches()[1], va="top")
    legend_h = 0.3 + 0.16 * math.ceil(len(handles) / 4)
    height = fig.get_size_inches()[1]
    fig.tight_layout(rect=(0, legend_h / height, 1, 1 - 0.3 / height))
    fig.subplots_adjust(top=1 - 0.55 / height)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def composition(plt, data, args, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    allv = [num(r[k]) for r in data["rows"] if r["drug"] == args.drug
            for k in ("llr", "llr_expected") if num(r[k]) is not None]
    top = max([abs(v) for v in allv] + [1.0]) * 1.4
    for ax, rd, title in ((axes[0], "all", "All branches"),
                          (axes[1], "internal", "Internal branches only")):
        style_axes(ax)
        ax.grid(color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.set_xscale("symlog", linthresh=1.0)
        ax.set_yscale("symlog", linthresh=1.0)
        ax.plot([0, top], [0, top], color=INK2, linewidth=0.8,
                linestyle=(0, (4, 3)), label="llr = llr_expected")
        ax.axhline(0, color=AXIS, linewidth=0.8)
        rows = [r for r in drug_rows(data, args.drug, rd).values()
                if num(r["llr"]) is not None]
        hit = [r for r in rows if num(r["llr_q"]) is not None
               and num(r["llr_q"]) < args.max_q]
        rest = [r for r in rows if r not in hit]

        def xy(rs):
            return ([num(r["llr_expected"]) for r in rs],
                    [num(r["llr"]) for r in rs])
        ax.scatter(*xy(rest), s=24, facecolors="none", edgecolors=INK2,
                   linewidths=0.8, label="tested split")
        ax.scatter(*xy(hit), s=40, color=HIGHLIGHT, edgecolors=SURFACE,
                   linewidths=1.0, zorder=3,
                   label=f"llr q < {args.max_q}")
        ax.set_xlim(0, top)
        ax.set_ylim(-top, top)
        label_points(fig, ax, [
            (num(r["llr_expected"]), num(r["llr"]),
             f"{r['split']} {data['lineage'].get(r['split'], ('', ''))[0]}"
             .rstrip()) for r in hit],
            obstacles=list(zip(*xy(rows))))
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("llr_expected: what the split's spectrum alone "
                      "would give", fontsize=8, color=INK2)
    axes[0].set_ylabel("llr: lean of the split's variants towards its "
                       "spectrum", fontsize=8, color=INK2)
    for ax, text, y in ((axes[0], "leans further than its spectrum "
                         "explains", 0.97),
                        (axes[0], "leans against its spectrum", 0.03)):
        ax.text(0.02, y, text, transform=ax.transAxes, fontsize=6.5,
                color=MUTED, va="top" if y > 0.5 else "bottom")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7,
               frameon=False, labelcolor=INK2)
    drug = "all drugs pooled" if args.drug == ANY else args.drug
    fig.suptitle(f"Resistance variants: does a split's spectrum predict "
                 f"which it acquires? ({drug})", fontsize=11, color=INK,
                 x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def example(plt, data, args, path):
    rd, drug, split = args.example_reading, args.example_drug, \
        args.example_split
    vs = [r for r in data["variants"] if r["reading"] == rd
          and r["drug"] == drug and r["split"] == split]
    fig, ax = plt.subplots(figsize=(7.5, 1.4 + 0.42 * max(
        1, min(len(vs), args.top_variants))))
    fig.patch.set_facecolor(SURFACE)
    style_axes(ax)
    if not vs:
        why = ("no variants table: rerun resistance_spectra.py"
               if not data["variants"] else
               f"{split} has no tested {drug} events on {rd} branches")
        ax.text(0.5, 0.5, why, transform=ax.transAxes, ha="center",
                fontsize=8, color=INK2)
        ax.set_axis_off()
        fig.savefig(path, facecolor=SURFACE)
        plt.close(fig)
        return
    keys = ("rest_share", "spectrum_share", "observed_share")
    vs.sort(key=lambda r: -max(num(r[k]) for k in keys))
    vs = vs[:args.top_variants]
    n = len(vs)
    h = 0.26
    for k, (key, c, lab) in enumerate((
            ("rest_share", REST, "rest of the tree (q)"),
            ("spectrum_share", PREDICTED,
             "rest, reweighted by the split's spectrum (q')"),
            ("observed_share", OBSERVED, "the split's own events"))):
        ys = [i + (k - 1) * h for i in range(n)]
        ax.barh(ys, [num(r[key]) for r in vs], height=h * 0.92, color=c,
                label=lab, edgecolor=SURFACE, linewidth=0.5)
    for i, r in enumerate(vs):
        ax.text(num(r["observed_share"]), i + h, f"  {r['observed']}",
                va="center", fontsize=6, color=INK2)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["variant"] for r in vs], fontsize=7, color=INK)
    ax.invert_yaxis()
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"share of the {drug} events", fontsize=8, color=INK2)
    row = drug_rows(data, drug, rd).get(split, {})
    lin = data["lineage"].get(split, ("", ""))[0]
    stats = ""
    if num(row.get("llr")) is not None:
        stats = (f": llr {num(row['llr']):.1f}, expected "
                 f"{num(row['llr_expected']):.1f}, q {num(row['llr_q']):.2g}")
    ax.set_title(f"{drug}, {split}{' ' + lin if lin else ''}, {rd} branches, "
                 f"{vs[0]['split_events']} events{stats}", fontsize=9,
                 color=INK, loc="left")
    ax.legend(loc="lower right", fontsize=7, frameon=False,
              labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


FIGURES = (("rate", rate), ("composition", composition),
           ("example", example))


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def _write(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _fixture(prefix, comp_path):
    import random
    rng = random.Random(1)
    cols = ["reading", "drug", "split", "events", "mutations", "e_global",
            "e_spectrum", "spectrum_factor", "obs_over_spectrum", "llr",
            "llr_expected", "llr_p", "llr_q", "g", "g_p", "g_q",
            "top_variant"]
    rows, vrows, comp = [], [], []
    lineages = ["lineage1.1", "lineage2.2.1", "lineage3", "lineage4.5",
                "La1.8", "lineage7", "(none)"]
    for i in range(20):
        s = f"node_{i}"
        comp.append({"node": s, "composition":
                     f"{lineages[i % len(lineages)]}:{100 + i},x:3"})
        for rd in ("all", "internal"):
            for drug in ("any drug", "Streptomycin"):
                ev = 0 if i == 0 else rng.randint(3, 900)
                eg = rng.uniform(50, 500)
                f = rng.uniform(0.93, 1.1)
                tested = ev >= 5
                rows.append({
                    "reading": rd, "drug": drug, "split": s, "events": ev,
                    "mutations": 1000, "e_global": eg, "e_spectrum": eg * f,
                    "spectrum_factor": f, "obs_over_spectrum": ev / (eg * f),
                    "llr": rng.uniform(-20, 60) if tested else "",
                    "llr_expected": rng.uniform(0, 50) if tested else "",
                    "llr_q": rng.uniform(0, 0.2) if tested else ""})
    for k in range(10):
        vrows.append({"reading": "internal", "drug": "Streptomycin",
                      "split": "node_3", "variant": f"rpsL_p.Lys{k}Arg",
                      "observed": 10 - k, "split_events": 55,
                      "observed_share": (10 - k) / 55,
                      "rest_share": 0.1, "spectrum_share": 0.08 + k / 200})
    _write(f"{prefix}.tsv", cols, rows)
    _write(f"{prefix}_variants.tsv", list(vrows[0]), vrows)
    _write(comp_path, ["node", "composition"], comp)


def selftest():
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    assert lineage_of("lineage2.2.1:10684,lineage2.2.1.1:2938") == \
        ("L2.2.1", "L2")
    assert lineage_of("La1.8.1:1826,La1.7.1:674")[1] == LINEAGE_GROUPS[4]
    assert lineage_of("lineage4.5:1444,(none):31") == ("L4.5", "L4")
    assert lineage_of("(none):5") == ("", "other")
    assert lineage_of("lineage9:3")[1] == "other"
    print("ok: lineage names and colour groups")
    with tempfile.TemporaryDirectory() as d:
        prefix = os.path.join(d, "rs")
        comp = os.path.join(d, "comp.tsv")
        _fixture(prefix, comp)
        for extra in ([], ["--composition", comp]):
            args = parse_args(["--prefix", prefix, "--example_split",
                               "node_3"] + extra)
            data = load(prefix, args.composition)
            for name, fn in FIGURES:
                out = f"{prefix}_{name}{len(extra)}.pdf"
                fn(plt, data, args, out)
                assert os.path.getsize(out) > 1000, name
        ### Without the variants table, the example says so.
        os.remove(f"{prefix}_variants.tsv")
        args = parse_args(["--prefix", prefix])
        example(plt, load(prefix), args, f"{prefix}_nv.pdf")
        assert os.path.getsize(f"{prefix}_nv.pdf") > 500
    print("ok: all three figures draw, with and without lineages, and the "
          "example draws without a variants table")
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
    data = load(args.prefix, args.composition)
    out = args.output or args.prefix
    written = []
    for name, fn in FIGURES:
        path = f"{out}_{name}.{args.format}"
        fn(plt, data, args, path)
        written.append(path)
    if not data["variants"]:
        print(f"note: no {args.prefix}_variants.tsv, so the example is "
              f"empty; rerun resistance_spectra.py to get it",
              file=sys.stderr)
    print(f"Wrote {', '.join(written)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
