#!/usr/bin/env python3
"""The punchline of gene_hit_contrasts.py's summary, in a page of text.

hit_summary.tsv has a row per gene and effect set and twenty-odd columns.
This reads it, and optionally the summaries of the internal-hits-only run
and of a random-gene control (step 12), and answers in order:

  1. Is the comparison sound?  The synonymous control: how many genes reach
     omnibus p < 0.05, and how many class tests reach p < 0.05, against what
     chance predicts, with a binomial tail p for the excess.
  2. Is there signal?  The same counts for the nonsynonymous hits.  With
     --control, the list's rate against the control genes' (one-sided
     Fisher), and the control's pooled shift per class: the baseline a hit
     of any kind shows, which a single gene's shift should be read against.
  3. Which genes clear the bar, 0.05 over the number of genes tested?  For
     each, every class at p < 0.05 (or the best one): its shift and p with
     all hits, with internal hits only (--internal), and in the gene's own
     synonymous control, with events and effective events for each.
  4. Which reach p < 0.05 without clearing it.  These are leads, and chance
     alone gives about one in twenty.

It also says when the smallest p is the permutation floor (so more
permutations could only lower it), when the bar is below that floor (so no
gene could clear it), and when a test rests on few effective events.

The number of permutations is read from the p-values themselves, which are
(k + 1) / (permutations + 1).

Usage:
  summarize_gene_hits.py --summary hit_summary.tsv \\
      --internal hit_summary_internal.tsv \\
      --control hit_summary_control.tsv > hit_summary.txt

  summarize_gene_hits.py --selftest
"""

import argparse
import csv
import math
import sys

NONSYN = "nonsynonymous"
SYN = "synonymous"
FEW_EFFECTIVE = 10


def read_summary(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def classes_of(rows):
    if not rows:
        return []
    return [c[len("diff_"):] for c in rows[0] if c.startswith("diff_")]


def genes(rows, es):
    """Per-gene rows (not ALL) of one effect set that have a test result."""
    return [r for r in rows if r["effect_set"] == es and r["gene"] != "ALL"
            and num(r.get("omnibus_p")) is not None]


def pooled(rows, es):
    for r in rows:
        if r["effect_set"] == es and r["gene"] == "ALL":
            return r
    return None


def permutations_from(rows):
    """Smallest N for which every p is (integer)/(N + 1), or None."""
    ps = [num(r[c]) for r in rows for c in r
          if (c.startswith("p_") or c == "omnibus_p") and num(r[c])]
    if not ps:
        return None
    for n in (100, 1000, 10000, 100000, 1000000):
        if all(abs(p * (n + 1) - round(p * (n + 1))) < 1e-3 * max(1, p * n)
               for p in ps):
            return n
    return None


def binom_tail(k, n, p):
    """P(X >= k), X ~ Binomial(n, p)."""
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i)
               for i in range(k, n + 1))


def fisher_greater(a, b, c, d):
    """One-sided Fisher p that row 1's rate a/(a+b) exceeds c/(c+d)."""
    n1, n2, k = a + b, c + d, a + c

    def h(x):
        return math.comb(n1, x) * math.comb(n2, k - x) / math.comb(n1 + n2, k)
    return sum(h(x) for x in range(a, min(n1, k) + 1))


def counts(rows, es, cls):
    g = genes(rows, es)
    n = len(g)
    o = sum(1 for r in g if num(r["omnibus_p"]) < 0.05)
    tests = [(r, c) for r in g for c in cls if num(r.get(f"p_{c}")) is not None]
    k = sum(1 for r, c in tests if num(r[f"p_{c}"]) < 0.05)
    return n, o, len(tests), k


def best_class(r, cls):
    """(class, shift in points, p) for the class with the smallest p."""
    best = None
    for c in cls:
        p, d = num(r.get(f"p_{c}")), num(r.get(f"diff_{c}"))
        if p is None or d is None:
            continue
        key = (p, -abs(d))
        if best is None or key < best[0]:
            best = (key, c, 100 * d, p)
    return best[1:] if best else (None, None, None)


def fmt_p(p):
    if p is None:
        return "NA"
    return f"{p:.1e}" if p < 0.001 else f"{p:.3g}"


def eff_str(r):
    e = num(r.get("effective_events")) if r else None
    return f"{e:,.0f}" if e is not None else "NA"


def by_gene(rows, es):
    return {r["gene"]: r for r in rows if r["effect_set"] == es}


def count_line(label, n, o, t, k, out):
    eo, ek = 0.05 * n, 0.05 * t
    print(f"  {label:<14} genes at omnibus p < 0.05: {o} of {n} "
          f"(chance {eo:.1f}, P(>= {o}) = {fmt_p(binom_tail(o, n, 0.05))});  "
          f"class tests at p < 0.05: {k} of {t} (chance {ek:.1f}, "
          f"P = {fmt_p(binom_tail(k, t, 0.05))})", file=out)


def report(main, args, out, internal=None, control=None):
    cls = classes_of(main)
    n_perm = permutations_from(main)
    floor = 1 / (n_perm + 1) if n_perm else None

    print(f"GENE HIT CONTRASTS: {args.summary}", file=out)
    for es in (NONSYN, SYN):
        a = pooled(main, es)
        if a:
            print(f"  {es}: {len(genes(main, es))} genes tested, "
                  f"{int(num(a['informative']) or 0):,} informative events"
                  + (f" ({eff_str(a)} effective)"
                     if num(a.get("effective_events")) else ""), file=out)
    if n_perm:
        print(f"  {n_perm:,} permutations (read from the p-values): the "
              f"smallest possible p is {fmt_p(floor)}", file=out)

    # 1. control ------------------------------------------------------------
    print(file=out)
    print("1. IS THE COMPARISON SOUND?  Synonymous hits should look like "
          "chance.", file=out)
    n, o, t, k = counts(main, SYN, cls)
    if n:
        count_line("synonymous", n, o, t, k, out)
        bad = binom_tail(o, n, 0.05) < 0.01 or binom_tail(k, t, 0.05) < 0.01
        print("  -> " + ("FAILS: the synonymous hits show more than chance, "
                         "so the comparison is biased and the nonsynonymous "
                         "rows cannot be read either."
                         if bad else "passes: no more than chance."), file=out)
    else:
        print("  no synonymous rows; run with --effects including "
              "synonymous", file=out)

    # 2. signal -------------------------------------------------------------
    print(file=out)
    print("2. IS THERE SIGNAL?  Nonsynonymous hits against chance"
          + (" and against the control genes." if control else "."),
          file=out)
    n, o, t, k = counts(main, NONSYN, cls)
    count_line("nonsynonymous", n, o, t, k, out)
    if control:
        cn, co, ct, ck = counts(control, NONSYN, cls)
        count_line("control genes", cn, co, ct, ck, out)
        print(f"  list vs control, one-sided Fisher: genes p = "
              f"{fmt_p(fisher_greater(o, n - o, co, cn - co))}, class tests "
              f"p = {fmt_p(fisher_greater(k, t - k, ck, ct - ck))}", file=out)
        for es in (NONSYN, SYN):
            a = pooled(control, es)
            if not a:
                continue
            shifts = ", ".join(f"{c} {100 * num(a[f'diff_{c}']):+.2f}"
                               for c in cls if num(a.get(f"diff_{c}")) is not None)
            print(f"  control baseline, any {es} hit (pooled, omnibus p "
                  f"{fmt_p(num(a['omnibus_p']))}): {shifts} points", file=out)

    # 3. genes that clear the bar ------------------------------------------
    ns = genes(main, NONSYN)
    bar = 0.05 / len(ns) if ns else None
    print(file=out)
    print(f"3. GENES THAT CLEAR THE BAR: omnibus p < 0.05 / {len(ns)} = "
          f"{fmt_p(bar)}", file=out)
    if floor and bar is not None and bar < floor:
        print(f"  NOTE: the bar is below the smallest p {n_perm:,} "
              f"permutations can give ({fmt_p(floor)}), so no gene can clear "
              f"it.  Rerun with more --permutations.", file=out)
    syn = by_gene(main, SYN)
    inn = by_gene(internal, NONSYN) if internal else {}
    clear = sorted((r for r in ns if num(r["omnibus_p"]) < bar),
                   key=lambda r: num(r["omnibus_p"]))
    if not clear:
        print("  none", file=out)
    for r in clear:
        at_floor = floor and num(r["omnibus_p"]) <= floor * 1.0001
        print(f"  {r['gene']}: omnibus p {fmt_p(num(r['omnibus_p']))}"
              + (" (the permutation floor; the true p may be smaller)"
                 if at_floor else ""), file=out)
        shown = sorted((c for c in cls if num(r.get(f"p_{c}")) is not None
                        and num(r[f"p_{c}"]) < 0.05),
                       key=lambda c: num(r[f"p_{c}"]))
        if not shown:
            shown = [best_class(r, cls)[0]]
        ri, rs = inn.get(r["gene"]), syn.get(r["gene"])
        cols = [("all hits", r)] + ([("internal hits", ri)] if ri else []) \
            + ([("synonymous", rs)] if rs else [])
        table = [["  class"] + [name for name, _x in cols]]
        for c in shown:
            line = [f"  {c}"]
            for _name, x in cols:
                d, p = num(x.get(f"diff_{c}")), num(x.get(f"p_{c}"))
                line.append(f"{100 * d:+.2f} pts, p {fmt_p(p)}"
                            if d is not None else "NA")
            table.append(line)
        table.append(["  events"] + [
            f"{int(num(x['informative']) or 0):,} ({eff_str(x)} effective)"
            for _name, x in cols])
        widths = [max(len(t[i]) for t in table) for i in range(len(table[0]))]
        for t in table:
            print("  " + "   ".join(v.ljust(w) for v, w in zip(t, widths))
                  .rstrip(), file=out)
        e = num(rs.get("effective_events")) if rs else None
        if e is not None and e < FEW_EFFECTIVE:
            print(f"    its synonymous control has only {e:.0f} effective "
                  f"events: a weak control; the one across genes (section 1) "
                  f"carries the weight", file=out)
        e = num(r.get("effective_events"))
        if e is not None and e < FEW_EFFECTIVE:
            print(f"    WARNING: only {e:.0f} effective events; the result "
                  f"is about those few branches", file=out)

    # 4. leads --------------------------------------------------------------
    leads = sorted((r for r in ns if bar <= num(r["omnibus_p"]) < 0.05),
                   key=lambda r: num(r["omnibus_p"]))
    print(file=out)
    print(f"4. LEADS: omnibus p < 0.05 but above the bar ({len(leads)}; "
          f"chance alone gives about {0.05 * len(ns):.1f})", file=out)
    for r in leads:
        c, d, p = best_class(r, cls)
        extra = ""
        ri = inn.get(r["gene"])
        if ri and num(ri.get(f"diff_{c}")) is not None:
            extra = (f"; internal {100 * num(ri[f'diff_{c}']):+.2f}, p "
                     f"{fmt_p(num(ri[f'p_{c}']))}")
        print(f"  {r['gene']:<10} omnibus p {fmt_p(num(r['omnibus_p'])):<8} "
              f"{c} {d:+.2f} points (p {fmt_p(p)}), {eff_str(r)} effective"
              f"{extra}", file=out)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Summarize gene_hit_contrasts.py output.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", help="hit_summary.tsv to report on")
    p.add_argument("--internal", default=None,
                   help="the same from a --top_branch exclude run")
    p.add_argument("--control", default=None,
                   help="the same for a control gene list (control_genes.py)")
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    if not args.selftest and not args.summary:
        p.error("--summary is required")
    return args


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    main_rows = read_summary(args.summary)
    internal = read_summary(args.internal) if args.internal else None
    control = read_summary(args.control) if args.control else None
    report(main_rows, args, sys.stdout, internal, control)
    return 0


def selftest():
    import io
    import random
    ok = True

    def check(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    rng = random.Random(3)
    cls = ["A", "B"]
    n_perm = 1000

    def p_of(k):
        return f"{(k + 1) / (n_perm + 1):.6g}"

    def row(es, gene, omni_k, a_k, a_diff, eff=50.0, inf=100):
        return {"effect_set": es, "gene": gene, "events": str(inf),
                "informative": str(inf), "informative_terminal": "0",
                "effective_events": str(eff),
                "diff_A": str(a_diff), "diff_B": "0.0",
                "p_A": p_of(a_k), "p_B": p_of(500),
                "omnibus": "1", "omnibus_p": p_of(omni_k)}

    rows = []
    for i in range(40):
        g = f"g{i}"
        if g == "g0":
            rows.append(row(NONSYN, g, 0, 0, 0.01))       # clears the bar
        elif g == "g1":
            rows.append(row(NONSYN, g, 20, 20, 0.02))     # a lead
        else:
            k = rng.randrange(100, 1000)
            rows.append(row(NONSYN, g, k, k, 0.0))
        rows.append(row(SYN, g, rng.randrange(100, 1000), 600, 0.0,
                        eff=2.0 if g == "g0" else 40.0))
    rows.append(row(NONSYN, "ALL", 500, 500, 0.001, eff=900, inf=4000))
    rows.append(row(SYN, "ALL", 500, 500, 0.0, eff=800, inf=3000))
    check(permutations_from(rows) == 1000, "permutations read from p-values")
    args = parse_args(["--selftest"])
    args.summary = "test.tsv"
    buf = io.StringIO()
    report(rows, args, buf)
    text = buf.getvalue()
    check("passes: no more than chance" in text, "synonymous control passes")
    check("g0: omnibus p" in text and "permutation floor" in text,
          "g0 clears the bar, at the floor")
    check("a weak control" in text, "g0's synonymous control flagged weak")
    check("\n  g1 " in text and "LEADS" in text, "g1 listed as a lead")
    check("bar is below the smallest p" not in text,
          "0.05 / 40 is within reach of 1,000 permutations")
    more = rows + [row(NONSYN, f"h{i}", 600, 600, 0.0) for i in range(20)]
    buf = io.StringIO()
    report(more, args, buf)
    text = buf.getvalue()
    check("bar is below the smallest p" in text and "g0: omnibus" not in text,
          "0.05 / 60 is not, and nothing clears it")
    check(abs(binom_tail(0, 10, 0.05) - 1.0) < 1e-12
          and abs(fisher_greater(5, 5, 0, 10) - 0.01625) < 1e-4,
          "binomial and Fisher tails")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
