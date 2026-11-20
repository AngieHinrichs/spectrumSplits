#!/usr/bin/env python3
"""Match the splits of two searches by sample membership.

For the "Testing the method" section of tuning/METHODS.md.  Node IDs do not
survive a mask, a prune or any other change to the tree, so the splits of a
second search are matched to the first by the samples they contain.  Each
input is a node<TAB>sample table, one per tree, made by running matUtils
extract -I on every node in that run's spectra.tsv (step 10's loop).

Both sides are first reduced to a shared universe: the samples present in
both trees, read from their node_1 rows (the root holds every sample).  A
sample pruned from the second tree would otherwise count against every old
clade that held it.  Then each old split is matched to the new split with the
highest Jaccard index (shared samples over samples in either), and the
reverse match says which new splits have no old counterpart.

Reported per old split: its best new split, the Jaccard index, the samples
only in the old and only in the new, and whether the match is mutual.  The
summary counts exact matches (Jaccard 1), near matches (>= --near) and
splits with no counterpart (best Jaccard < --lost).

Usage:
  match_splits.py --old split_clade_samples.tsv \\
      --new split_clade_samples_nobcg.tsv --output split_matches.tsv

  match_splits.py --selftest
"""

import argparse
import collections
import sys

ROOT = "node_1"


def read_sets(lines):
    sets = collections.defaultdict(set)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] and not f[0].startswith("#"):
            sets[f[0]].add(f[1])
    return dict(sets)


def universe(old, new):
    """Samples in both trees, from their roots, or the union of all splits."""
    def every(sets):
        if ROOT in sets:
            return sets[ROOT]
        out = set()
        for s in sets.values():
            out |= s
        return out
    return every(old) & every(new)


def best_matches(a, b):
    """For each node in a: (best node in b, intersection, Jaccard)."""
    member = collections.defaultdict(list)
    for node, samples in b.items():
        for s in samples:
            member[s].append(node)
    out = {}
    for node, samples in a.items():
        inter = collections.Counter()
        for s in samples:
            for m in member.get(s, ()):
                inter[m] += 1
        best = (None, 0, 0.0)
        for m, k in inter.items():
            j = k / (len(samples) + len(b[m]) - k)
            if j > best[2] or (j == best[2] and best[0] and m < best[0]):
                best = (m, k, j)
        out[node] = best
    return out


def match(old, new):
    u = universe(old, new)
    o = {n: s & u for n, s in old.items() if n != ROOT}
    w = {n: s & u for n, s in new.items() if n != ROOT}
    o = {n: s for n, s in o.items() if s}
    w = {n: s for n, s in w.items() if s}
    fwd, back = best_matches(o, w), best_matches(w, o)
    rows = []
    for n in sorted(o, key=lambda n: -len(o[n])):
        m, k, j = fwd[n]
        rows.append({"old": n, "old_samples": len(o[n]), "new": m or "",
                     "new_samples": len(w[m]) if m else 0, "jaccard": j,
                     "only_old": len(o[n]) - k,
                     "only_new": (len(w[m]) - k) if m else 0,
                     "mutual": bool(m) and back[m][0] == n})
    unmatched_new = sorted(m for m in w if back[m][2] < 1.0
                           and not any(r["new"] == m for r in rows))
    return rows, unmatched_new, len(u), back, w


def report(rows, unmatched_new, n_universe, back, new_sets, args, out):
    exact = sum(1 for r in rows if r["jaccard"] == 1.0)
    near = sum(1 for r in rows if args.near <= r["jaccard"] < 1.0)
    lost = [r for r in rows if r["jaccard"] < args.lost]
    print(f"SPLIT MATCHING by sample membership, on the {n_universe:,} "
          f"samples present in both trees", file=out)
    print(f"  {len(rows)} old splits, {len(new_sets)} new", file=out)
    print(f"  Jaccard 1.000: {exact}", file=out)
    print(f"  Jaccard >= {args.near} but < 1: {near}", file=out)
    print(f"  in between: {len(rows) - exact - near - len(lost)}", file=out)
    print(f"  no counterpart (best Jaccard < {args.lost}): {len(lost)}"
          + (": " + " ".join(r["old"] for r in lost) if lost else ""),
          file=out)
    print(file=out)
    table = [["old", "samples", "best new", "samples", "Jaccard",
              "only old", "only new", "mutual"]]
    for r in rows:
        table.append([r["old"], f"{r['old_samples']:,}", r["new"],
                      f"{r['new_samples']:,}", f"{r['jaccard']:.3f}",
                      str(r["only_old"]), str(r["only_new"]),
                      "yes" if r["mutual"] else "no"])
    widths = [max(len(t[i]) for t in table) for i in range(len(table[0]))]
    for t in table:
        print("  " + "  ".join(v.ljust(w) for v, w in zip(t, widths)).rstrip(),
              file=out)
    new_only = [m for m in new_sets
                if not any(r["new"] == m and r["mutual"] for r in rows)]
    print(file=out)
    print(f"  new splits that are no old split's mutual match: {len(new_only)}"
          + (": " + " ".join(f"{m} (best {back[m][0] or '-'}, "
                             f"{back[m][2]:.3f})" for m in sorted(new_only))
             if new_only else ""), file=out)


def write_tsv(rows, path):
    cols = ["old", "old_samples", "new", "new_samples", "jaccard", "only_old",
            "only_new", "mutual"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = [f"{r[c]:.4f}" if c == "jaccard" else str(r[c])
                    for c in cols]
            fh.write("\t".join(vals) + "\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Match two searches' splits by sample membership.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--old", help="node<TAB>sample for the first search")
    p.add_argument("--new", help="node<TAB>sample for the second search")
    p.add_argument("--near", type=float, default=0.95)
    p.add_argument("--lost", type=float, default=0.5)
    p.add_argument("--output", default=None)
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    if not args.selftest and not (args.old and args.new):
        p.error("--old and --new are required")
    return args


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    with open(args.old) as fh:
        old = read_sets(fh)
    with open(args.new) as fh:
        new = read_sets(fh)
    rows, unmatched, n_u, back, w = match(old, new)
    if args.output:
        write_tsv(rows, args.output)
    report(rows, unmatched, n_u, back, w, args, sys.stdout)
    return 0


def selftest():
    ok = True

    def check(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    ### Old tree: 20 samples.  node_2 = s0-9, node_3 = s0-2, node_4 = s10-14.
    ### New tree: s0 and s1 pruned, node IDs all changed; node_9 = s2-9 (old
    ### node_2 less the pruned), node_8 = s2 (old node_3 less the pruned),
    ### node_7 = s10-13 (old node_4 less one sample), node_6 a new split.
    old = {"node_1": {f"s{i}" for i in range(20)},
           "node_2": {f"s{i}" for i in range(10)},
           "node_3": {f"s{i}" for i in range(3)},
           "node_4": {f"s{i}" for i in range(10, 15)}}
    new = {"node_1": {f"s{i}" for i in range(2, 20)},
           "node_9": {f"s{i}" for i in range(2, 10)},
           "node_8": {"s2"},
           "node_7": {f"s{i}" for i in range(10, 14)},
           "node_6": {f"s{i}" for i in range(15, 20)}}
    rows, unmatched, n_u, back, w = match(old, new)
    by = {r["old"]: r for r in rows}
    check(n_u == 18, f"universe is the samples in both: {n_u}")
    check(by["node_2"]["new"] == "node_9" and by["node_2"]["jaccard"] == 1.0,
          "a pruned sample does not count against its old clade")
    check(by["node_3"]["new"] == "node_8" and by["node_3"]["jaccard"] == 1.0,
          "small clade matched after the prune")
    check(by["node_4"]["new"] == "node_7"
          and abs(by["node_4"]["jaccard"] - 0.8) < 1e-12
          and by["node_4"]["only_old"] == 1, "partial match counted")
    check(all(r["mutual"] for r in rows), "mutual matches")
    args = parse_args(["--selftest"])
    import io
    buf = io.StringIO()
    report(rows, unmatched, n_u, back, w, args, buf)
    text = buf.getvalue()
    check("Jaccard 1.000: 2" in text and "node_6" in text,
          "report counts exact matches and lists the new-only split")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
