#!/usr/bin/env python3
"""The splits whose context shift survives both of split_contexts.py's runs.

Step 11 runs split_contexts.py twice, on all branches and on internal
branches only, and reads the two together.  This derives the table that
reading ends in: every split that passes both nulls -- q at or below --max_q
from the branch shuffle, and ctrl_p at or below --max_ctrl_p from the random
clades -- in BOTH runs, sorted by its all-branch excess.  --top N keeps the N
with the lowest all-branch ctrl_p, as the README's table does.

Per split it writes:

  clade              with --composition, split_composition.py's top clades
                     of the split's residual samples, and its flags (mixed)
  excess_all, excess_internal, ctrl_p_all, ctrl_p_internal, q_all,
  q_internal         from the two split_contexts.py --output tables
  contexts           with --contexts, the --per_split contexts carrying most
                     of G(within) on all branches, each as its share of its
                     class, split vs comparator, from split_contexts.py
                     --output_contexts; without it, split_contexts.py's own
                     top_contexts column (labels and G only)
  top_site_shares    with --contexts, for each of those contexts the share
                     of the split's count from its single most frequent
                     site, and max_site_share the largest of them.  A
                     context at --one_site or above is named on stderr: one
                     site standing in for a context

Usage:
    summarize_split_contexts.py --all split_contexts.tsv \\
        --internal split_contexts_internal.tsv \\
        --contexts split_contexts_96.tsv --composition split_composition.tsv \\
        --top 12 > split_contexts_table.tsv

    summarize_split_contexts.py --selftest
"""

import argparse
import csv
import io
import sys


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Tabulate the splits whose context shift passes both "
                    "nulls on all branches and on internal branches.")
    p.add_argument("--all", help="split_contexts.py --output, all branches")
    p.add_argument("--internal",
                   help="split_contexts.py --output, --branches internal")
    p.add_argument("--contexts", default=None,
                   help="split_contexts.py --output_contexts of the "
                        "all-branch run, for the class shares")
    p.add_argument("--composition", default=None,
                   help="split_composition.py --output, for each split's "
                        "clades")
    p.add_argument("--max_q", type=float, default=0.05,
                   help="Branch-shuffle q a split must be below, in both "
                        "runs (default 0.05)")
    p.add_argument("--max_ctrl_p", type=float, default=0.05,
                   help="Random-clade ctrl_p a split must be at or below, in "
                        "both runs (default 0.05)")
    p.add_argument("--top", type=int, default=0,
                   help="Keep only the N passing splits with the lowest "
                        "all-branch ctrl_p (default 0: all of them)")
    p.add_argument("--per_split", type=int, default=2,
                   help="Contexts listed per split (default 2)")
    p.add_argument("--one_site", type=float, default=0.2,
                   help="Top-site share at which a listed context is "
                        "flagged as one site (default 0.2, as in "
                        "split_contexts.py's report)")
    p.add_argument("--output", default=None,
                   help="TSV to write (default stdout)")
    p.add_argument("--selftest", action="store_true",
                   help="Run the self-tests and exit")
    return p.parse_args(argv)


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def p2(v):
    x = num(v)
    return "" if x is None else f"{x:.2g}"


def read_table(fh, key="node"):
    return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def passes(r, args):
    q, c = num(r.get("q")), num(r.get("control_p"))
    return q is not None and q < args.max_q \
        and c is not None and c <= args.max_ctrl_p


def read_contexts(fh):
    """node -> list of context rows, highest G contribution first."""
    out = {}
    for r in csv.DictReader(fh, delimiter="\t"):
        out.setdefault(r["node"], []).append(r)
    for rows in out.values():
        rows.sort(key=lambda r: -float(r["g_contribution"]))
    return out


def clade_label(r):
    comp = r.get("composition", "")
    parts = []
    for item in comp.split(","):
        if ":" in item:
            name, n = item.rsplit(":", 1)
            parts.append(f"{name} {n}")
    label = ", ".join(parts)
    flags = r.get("flags", "")
    return f"{label} ({flags})" if flags else label


def build(args, all_rows, int_rows, ctx=None, comp=None):
    both = [n for n, r in all_rows.items()
            if passes(r, args) and n in int_rows and passes(int_rows[n], args)]
    only_all = sorted(n for n, r in all_rows.items() if passes(r, args)
                      and n not in both)
    only_int = sorted(n for n, r in int_rows.items() if passes(r, args)
                      and n not in both)
    chosen = sorted(both, key=lambda n: num(all_rows[n]["control_p"]))
    if args.top:
        chosen = chosen[:args.top]
    chosen.sort(key=lambda n: -num(all_rows[n]["excess"]))
    out = []
    one_site = []
    for n in chosen:
        a, i = all_rows[n], int_rows[n]
        sites, max_site = "", ""
        if ctx is not None and n in ctx:
            listed = [c for c in ctx[n][:args.per_split]
                      if float(c["g_contribution"]) > 0]
            contexts = "; ".join(
                f"{c['context']} {float(c['split_class_share']):.3f} vs "
                f"{float(c['comparator_class_share']):.3f}" for c in listed)
            shares = [num(c.get("split_top_site_share")) for c in listed]
            if listed and all(x is not None for x in shares):
                sites = "; ".join(f"{x:.3f}" for x in shares)
                max_site = f"{max(shares):.3f}"
                one_site += [(n, c["context"], x, c.get("split_top_site"))
                             for c, x in zip(listed, shares)
                             if x >= args.one_site]
        else:
            contexts = ",".join(a.get("top_contexts", "").split(",")
                                [:args.per_split])
        out.append({
            "split": n,
            "comparator": a["comparator"],
            "clade": clade_label(comp[n]) if comp and n in comp else "",
            "excess_all": f"{num(a['excess']):.2f}",
            "excess_internal": f"{num(i['excess']):.2f}",
            "ctrl_p_all": p2(a["control_p"]),
            "ctrl_p_internal": p2(i["control_p"]),
            "q_all": p2(a["q"]),
            "q_internal": p2(i["q"]),
            "contexts": contexts,
            "top_site_shares": sites,
            "max_site_share": max_site,
        })
    return out, len(both), only_all, only_int, one_site


COLUMNS = ["split", "comparator", "clade", "excess_all", "excess_internal",
           "ctrl_p_all", "ctrl_p_internal", "q_all", "q_internal", "contexts",
           "top_site_shares", "max_site_share"]


def write(rows, fh):
    fh.write("\t".join(COLUMNS) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")


def selftest():
    hdr = "node\tcomparator\texcess\tq\tcontrol_p\ttop_contexts\n"
    all_tsv = hdr + "".join([
        "a\tr\t5.0\t0.001\t0.003\tX+:10,Y-:5,Z+:1\n",   # passes both
        "b\tr\t7.0\t0.001\t0.020\tX+:9\n",              # passes both
        "c\tr\t3.0\t0.001\t0.040\tX+:3\n",              # passes both, worst ctrl
        "d\tr\t9.0\t0.001\t0.003\tX+:3\n",              # fails internal
        "e\tr\t2.0\t0.200\t0.003\tX+:3\n",              # fails q
        "f\tr\t2.0\t0.001\t0.300\tX+:3\n",              # fails ctrl_p
    ])
    int_tsv = hdr + "".join([
        "a\tr\t2.5\t0.001\t0.005\t\n",
        "b\tr\t3.0\t0.010\t0.040\t\n",
        "c\tr\t1.5\t0.020\t0.050\t\n",
        "d\tr\t1.0\t0.500\t0.600\t\n",
        "e\tr\t1.0\t0.001\t0.003\t\n",
        "f\tr\t1.0\t0.001\t0.003\t\n",
    ])
    ctx_tsv = ("node\tcomparator\tclass\tcontext\tsplit_count\t"
               "comparator_count\tsplit_class_share\tcomparator_class_share\t"
               "g_contribution\tsplit_top_site\tsplit_top_site_share\n"
               "a\tr\tCT|GA\tA[C>T]G\t1\t1\t0.2\t0.1\t3.0\t17\t0.25\n"
               "a\tr\tCT|GA\tC[C>T]G\t1\t1\t0.3\t0.1\t9.0\t40\t0.01\n"
               "a\tr\tCT|GA\tG[C>T]G\t1\t1\t0.1\t0.1\t0.0\t5\t0.9\n")
    comp_tsv = ("node\ttier\tsamples\tresidual\tclades\tother_lineage\tflags\t"
                "composition\n"
                "a\tA\t10\t10\t2\t5\tmixed\tL1:6,L2:4\n")
    args = parse_args([])
    rows, nboth, only_all, only_int, one = build(
        args, read_table(io.StringIO(all_tsv)),
        read_table(io.StringIO(int_tsv)),
        read_contexts(io.StringIO(ctx_tsv)),
        read_table(io.StringIO(comp_tsv)))
    assert [r["split"] for r in rows] == ["b", "a", "c"], rows
    assert nboth == 3 and only_all == ["d"] and only_int == ["e", "f"]
    a = [r for r in rows if r["split"] == "a"][0]
    assert a["contexts"] == "C[C>T]G 0.300 vs 0.100; A[C>T]G 0.200 vs 0.100", a
    assert a["clade"] == "L1 6, L2 4 (mixed)", a["clade"]
    b = [r for r in rows if r["split"] == "b"][0]
    assert b["contexts"] == "X+:9" and b["clade"] == ""
    ### Site shares follow the listed contexts; the unlisted 0.9 is ignored.
    assert a["top_site_shares"] == "0.010; 0.250", a
    assert a["max_site_share"] == "0.250" and b["max_site_share"] == ""
    assert one == [("a", "A[C>T]G", 0.25, "17")], one
    print("ok: a split must pass both nulls in both runs; rows sort by "
          "all-branch excess; contexts come from the 96 table, highest G "
          "first, with their top-site shares, and one site is flagged")
    args.top = 2
    rows, _n, _a, _i, _o = build(args, read_table(io.StringIO(all_tsv)),
                             read_table(io.StringIO(int_tsv)))
    assert [r["split"] for r in rows] == ["b", "a"], rows
    print("ok: --top keeps the lowest all-branch ctrl_p, then sorts by excess")
    print("all self-tests passed")
    return 0


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    if not args.all or not args.internal:
        print("error: --all and --internal are required (or pass "
              "--selftest)", file=sys.stderr)
        return 2
    with open(args.all) as fh:
        all_rows = read_table(fh)
    with open(args.internal) as fh:
        int_rows = read_table(fh)
    ctx = comp = None
    if args.contexts:
        with open(args.contexts) as fh:
            ctx = read_contexts(fh)
    if args.composition:
        with open(args.composition) as fh:
            comp = read_table(fh)
    rows, nboth, only_all, only_int, one_site = build(
        args, all_rows, int_rows, ctx, comp)
    print(f"{nboth} splits have q < {args.max_q:g} and ctrl_p <= "
          f"{args.max_ctrl_p:g} on both runs; {len(rows)} listed"
          f"{f' (--top {args.top})' if args.top else ''}.", file=sys.stderr)
    if only_all:
        print(f"  pass on all branches only: {', '.join(only_all)}",
              file=sys.stderr)
    if only_int:
        print(f"  pass on internal branches only: {', '.join(only_int)}",
              file=sys.stderr)
    if ctx is not None:
        has_sites = any(r["max_site_share"] for r in rows)
        if rows and not has_sites:
            print("  no top-site shares in --contexts: rerun "
                  "split_contexts.py to get them", file=sys.stderr)
        elif one_site:
            print(f"  listed contexts where one site carries >= "
                  f"{args.one_site:g} of the split's count:", file=sys.stderr)
            for n, c, x, pos in one_site:
                print(f"    {n} {c}: {x:.2f} at {pos}", file=sys.stderr)
        elif has_sites:
            big = max(float(r["max_site_share"]) for r in rows
                      if r["max_site_share"])
            print(f"  no listed context has one site at >= "
                  f"{args.one_site:g}; the largest top-site share is "
                  f"{big:.3f}.", file=sys.stderr)
    if args.output:
        with open(args.output, "w") as fh:
            write(rows, fh)
    else:
        write(rows, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
