#!/usr/bin/env python3
"""Tabulate what qc/prune_mutation_sample_ratio.py removed, by clade.

This builds the tables in README.md step 2 ("Check what it removed, by clade"
and "Where the pruned samples sat") from files that step already produces:

    matUtils summary -i tree.pb.gz -s samples_raw.txt
    matUtils summary -i tree.pb.gz --sample-clades clades_raw.txt
    qc/prune_mutation_sample_ratio.py ... > pruned_report.txt

All three describe the UNPRUNED tree (plus the prune report); the pruned tree
is not needed, because the prune report lists every removed sample.

Tables written (markdown, to stdout):

  1. Removal rate by clade: in tree, removed, rate, and rate relative to the
     overall rate.
  2. Terminal branches vs removals: median terminal branch length (the
     sample's parsimony score from `matUtils summary -s`), how many tips have
     a terminal branch >= threshold, how many were pruned, and how many of
     those went by a terminal-branch crossing vs. an internal-node crossing.
     A clade with pruned > tips >= threshold lost samples through an internal
     node -- the stem-driven behaviour the README describes.
  3. Internal-node prunings: each pruned subtree rooted at an internal node,
     with its subtree mutation total (ratio x tips), the sum of its tips'
     terminal branches, and the remainder (stem + internal branches).  With
     --input_tree the stem is reported exactly.
  4. (needs --input_tree) Where the pruned samples sat: for each surviving
     parent of a pruned subtree, the number of tips under it and its
     commonest clades.
  5. Unlabelled samples: how many parents hold an unlabelled sample, how many
     of those also hold a labelled one, and the contiguous runs of unlabelled
     samples in node-number order with the labelled lineage on either side.
     Usher numbers nodes in traversal order, so a clade occupies a contiguous
     range of parent node numbers.

Clades are collapsed to --depth dot-separated components (default 1, so
lineage4.3.4.2 counts as lineage4).  Use --depth 0 for full names.

Usage:
    prune_clade_check.py --samples samples_raw.txt --clades clades_raw.txt \\
        --prune_report pruned_report.txt --threshold 250
    prune_clade_check.py ... --input_tree tree.pb.gz     # adds tables 3-4 detail
"""

import argparse
import statistics
import sys
from collections import Counter, defaultdict

UNLABELLED = "None"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split("\n", 1)[1])
    p.add_argument("--samples", required=True,
                   help="`matUtils summary -s` output for the unpruned tree")
    p.add_argument("--clades", required=True,
                   help="`matUtils summary --sample-clades` output for the unpruned tree")
    p.add_argument("--prune_report", required=True,
                   help="stdout of qc/prune_mutation_sample_ratio.py")
    p.add_argument("--threshold", type=float, required=True,
                   help="the --threshold given to prune_mutation_sample_ratio.py")
    p.add_argument("--depth", type=int, default=1,
                   help="collapse clade names to this many dot-separated parts "
                        "(0 = full name; default 1)")
    p.add_argument("--min_run", type=int, default=5,
                   help="minimum size of an unlabelled run to list (default 5)")
    p.add_argument("--input_tree",
                   help="unpruned tree; adds exact stems and parent subtree sizes (needs bte)")
    p.add_argument("--exclude", default="ReconstructedAncestor",
                   help="comma-separated sample names to ignore (default: ReconstructedAncestor)")
    return p.parse_args()


def node_number(node_id):
    return int(node_id.rsplit("_", 1)[1]) if node_id.startswith("node_") else -1


def collapse(clade, depth):
    if depth <= 0 or clade == UNLABELLED:
        return clade
    return ".".join(clade.split(".")[:depth])


def read_samples(path, exclude):
    terminal, parent = {}, {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:3] != ["sample", "parsimony", "parent_id"]:
            sys.exit(f"error: {path}: expected 'sample parsimony parent_id' header, got {header}")
        for line in f:
            name, pars, par = line.rstrip("\n").split("\t")[:3]
            if name in exclude:
                continue
            terminal[name] = int(pars)
            parent[name] = par
    return terminal, parent


def read_clades(path, exclude):
    clade = {}
    with open(path) as f:
        f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in exclude:
                continue
            clade[fields[0]] = fields[1] if len(fields) > 1 and fields[1] else UNLABELLED
    return clade


def read_prune_report(path):
    events = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            par, child, tips, ratio = line.rstrip("\n").split("\t")[:4]
            events.append({"parent": par, "child": child, "tips": tips.split(","),
                           "ratio": float(ratio)})
    return events


def pct(x):
    if x == 0:
        return "0"
    if x >= 0.1:
        return f"{100 * x:.0f}%" if x >= 0.995 else f"{100 * x:.1f}%"
    return f"{100 * x:.2g}%"


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    print("\n".join(out))
    print()


def fmt_int(n):
    return f"{n:,}"


def main():
    args = parse_args()
    exclude = set(filter(None, args.exclude.split(",")))
    terminal, parent = read_samples(args.samples, exclude)
    full_clade = read_clades(args.clades, exclude)
    events = read_prune_report(args.prune_report)

    missing = set(terminal) - set(full_clade)
    if missing:
        print(f"warning: {len(missing)} samples have no clade; counted as {UNLABELLED}",
              file=sys.stderr)
    full_clade = {s: full_clade.get(s, UNLABELLED) for s in terminal}
    clade = {s: collapse(c, args.depth) for s, c in full_clade.items()}

    pruned_by = {}  # sample -> "terminal" | "internal"
    for e in events:
        how = "terminal" if e["child"] in terminal else "internal"
        for s in e["tips"]:
            if s not in terminal:
                print(f"warning: pruned sample {s} not in {args.samples}", file=sys.stderr)
                continue
            pruned_by[s] = how

    in_tree = Counter(clade.values())
    removed = Counter(clade[s] for s in pruned_by)
    removed_term = Counter(clade[s] for s, h in pruned_by.items() if h == "terminal")
    removed_int = Counter(clade[s] for s, h in pruned_by.items() if h == "internal")
    by_clade_terms = defaultdict(list)
    for s, t in terminal.items():
        by_clade_terms[clade[s]].append(t)

    n_total = len(terminal)
    n_removed = len(pruned_by)
    overall = n_removed / n_total
    print(f"# Pruning check: {n_removed:,} of {n_total:,} samples removed "
          f"at threshold {args.threshold:g} (overall rate {pct(overall)})\n")

    # ---- 1. removal rate by clade ----------------------------------------
    print("## 1. Removal rate by clade\n")
    order = sorted(in_tree, key=lambda c: (-removed[c] / in_tree[c], -in_tree[c]))
    rows = []
    for c in order:
        rate = removed[c] / in_tree[c]
        rows.append([c, fmt_int(in_tree[c]), removed[c], pct(rate),
                     f"{rate / overall:.2g}x" if removed[c] and overall else "—"])
    table(["clade", "in tree", "removed", "rate", "vs overall"], rows)

    # ---- 2. terminal branches vs removals --------------------------------
    print(f"## 2. Terminal branches >= {args.threshold:g} vs. removals\n")
    print("A clade whose `pruned` exceeds `tips >= threshold` lost samples through an\n"
          "internal node (whole-subtree pruning), not through its own long tips.\n")
    rows = []
    for c in sorted(in_tree, key=lambda c: -in_tree[c]):
        terms = by_clade_terms[c]
        n_long = sum(1 for t in terms if t >= args.threshold)
        flag = " **!**" if removed[c] > n_long else ""
        rows.append([c, fmt_int(in_tree[c]), f"{statistics.median(terms):g}", n_long,
                     f"{removed[c]}{flag}", removed_term[c], removed_int[c]])
    table(["clade", "in tree", "median terminal", f"tips >= {args.threshold:g}",
           "pruned", "via terminal", "via internal node"], rows)

    # ---- optional tree ---------------------------------------------------
    tree = None
    if args.input_tree:
        try:
            import bte
        except ImportError:
            sys.exit("error: --input_tree needs bte; run in the bte environment")
        print(f"Reading {args.input_tree} ...", file=sys.stderr)
        tree = bte.MATree(args.input_tree)

    # ---- 3. internal-node prunings ---------------------------------------
    internal = [e for e in events if e["child"] not in terminal]
    print("## 3. Prunings at internal nodes\n")
    print("subtree total = ratio x tips (stem included); remainder = total minus the\n"
          "tips' own terminal branches, i.e. the stem plus any internal branches.\n")
    rows = []
    for e in internal:
        tips = [s for s in e["tips"] if s in terminal]
        total = round(e["ratio"] * len(e["tips"]))
        tip_terms = [terminal[s] for s in tips]
        rem = total - sum(tip_terms)
        clades = Counter(full_clade[s] for s in tips)
        row = [e["parent"], e["child"], len(e["tips"]), f"{e['ratio']:g}", total,
               ", ".join(map(str, sorted(tip_terms))), rem]
        if tree is not None:
            stem = len(tree.get_node(e["child"]).mutations)
            row += [stem, rem - stem]
        row.append(", ".join(f"{k} ({v})" if v > 1 else k for k, v in clades.most_common()))
        rows.append(row)
    hdr = ["parent", "child", "tips", "ratio", "subtree total", "tip terminals", "remainder"]
    if tree is not None:
        hdr += ["stem", "internal below"]
    table(hdr + ["clades"], rows)

    # ---- 4. where the pruned samples sat ---------------------------------
    if tree is not None:
        print("## 4. Surviving parents of pruned subtrees\n")
        print("Tips are counted in the unpruned tree, so they include the pruned ones.\n")
        by_parent = defaultdict(list)
        for e in events:
            by_parent[e["parent"]].extend(e["tips"])
        rows = []
        for par, ptips in by_parent.items():
            leaves = [s for s in tree.get_leaves_ids(par) if s in full_clade]
            cc = Counter(full_clade[s] for s in leaves)
            top = ", ".join(f"`{k}` {v:,}" for k, v in cc.most_common(3))
            if len(cc) > 3:
                top += f", +{len(cc) - 3} more"
            pcl = Counter(full_clade[s] for s in ptips if s in full_clade)
            rows.append([par, len(ptips),
                         ", ".join(f"{k} ({v})" if v > 1 else k for k, v in pcl.most_common()),
                         fmt_int(len(leaves)), top])
        rows.sort(key=lambda r: node_number(r[0]))
        table(["parent", "pruned", "pruned clades", "tips under parent", "commonest clades"],
              rows)

    # ---- 5. unlabelled samples -------------------------------------------
    print("## 5. Unlabelled samples\n")
    par_labels = defaultdict(set)
    for s, par in parent.items():
        par_labels[par].add(full_clade[s] == UNLABELLED)
    with_unl = [p for p, v in par_labels.items() if True in v]
    mixed = [p for p in with_unl if False in par_labels[p]]

    ordered = sorted(parent, key=lambda s: (node_number(parent[s]), full_clade[s] != UNLABELLED, s))
    runs, i = [], 0
    while i < len(ordered):
        if full_clade[ordered[i]] != UNLABELLED:
            i += 1
            continue
        j = i
        while j < len(ordered) and full_clade[ordered[j]] == UNLABELLED:
            j += 1
        runs.append((i, j))
        i = j
    n_unl = sum(1 for c in full_clade.values() if c == UNLABELLED)
    big = [r for r in runs if r[1] - r[0] >= args.min_run]
    covered = sum(j - i for i, j in big)
    n_pruned_unl = sum(1 for s in pruned_by if full_clade[s] == UNLABELLED)
    table(["", ""], [
        ["unlabelled samples", fmt_int(n_unl)],
        ["...of which pruned", n_pruned_unl],
        ["parents holding at least one unlabelled sample", fmt_int(len(with_unl))],
        ["...that also hold a labelled sample", len(mixed)],
        [f"contiguous unlabelled runs of {args.min_run}+ samples",
         f"{len(big)}, covering {covered:,} of {n_unl:,}"],
    ])
    rows = []
    for i, j in sorted(big, key=lambda r: -(r[1] - r[0])):
        lo, hi = node_number(parent[ordered[i]]), node_number(parent[ordered[j - 1]])
        below = full_clade[ordered[i - 1]] if i > 0 else "(start)"
        above = full_clade[ordered[j]] if j < len(ordered) else "(end)"
        n_pr = sum(1 for s in ordered[i:j] if s in pruned_by)
        rows.append([j - i, f"{lo}..{hi}", f"`{below}`", f"`{above}`", n_pr])
    table(["size", "node range", "below", "above", "pruned"], rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
