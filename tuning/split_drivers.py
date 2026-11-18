#!/usr/bin/env python3
"""Which mutation classes drive each split, and which splits to interpret.

This belongs to step 8 of parameter tuning, once spectrumSplits.py has
finished.  spectra.tsv gives each split's per-class proportions, but reading
sixty rows of eighteen columns by eye is how a one-class artifact gets
mistaken for biology.  This script does the reading:

  1. WHICH CLASSES DRIVE IT.  Each split's residual spectrum (its row in
     spectra.tsv) is set against its parent's residual row -- the parent is
     the subtree it was found in, read from splits.log -- as a 2 x K
     chi-squared table, and each class's share of that chi-squared is
     reported.  A split whose top class carries --one_class_share or more of
     it is driven by one class, which usually means an artifact or a single
     over-mutated site.  The folded classes are used when spectra.tsv has
     them (a --fold run), otherwise the twelve.

     The search scored each split against its parent subtree at the moment
     it was accepted, before later splits were carved out of either side, so
     chi2_residual is not the logged chi-squared.  It answers a different
     question: which classes separate the rows you will actually read.

  2. ONE SITE?  With --input_tree, each split's residual mutations are
     counted by site, and for its --top_sites most recurrent sites the
     chi-squared is recomputed with that site removed from both rows.  If one
     removal takes away --one_site_loss or more of it, one site is the split.
     The residual totals are also checked against spectra.tsv, which catches
     a different tree or --max_branch_length.

  3. BOOTSTRAP SUPPORT.  With --bootstrap, the fraction of replicates that
     split the same node or one within --neighbourhood of its number (as in
     split_confidence.py, whose step 10 ranking this anticipates).

  4. THE DEFINING BRANCH.  With --input_tree, the mutations on each split's
     own branch, as a comma-separated list annotate_sites.py --positions
     takes, so the first look for a mechanism is one command.

A split is listed to interpret when it is multi-class, not one site, and
supported.  Each test is skipped, not failed, when its input is not given.

Usage:
  split_drivers.py --spectra spectra.tsv --splits_log splits.log \\
      --bootstrap bootstrap_*_splits_output.tsv \\
      --input_tree tree.pruned.masked.pb.gz --max_branch_length 250 \\
      --output split_drivers.tsv > split_drivers.txt

  split_drivers.py --selftest
"""

import argparse
import collections
import math
import os
import re
import sys

from split_confidence import parse_splits_log, neighbourhood_support

MUTATION_TYPES = ["AC", "AG", "AT", "CA", "CG", "CT",
                  "GA", "GC", "GT", "TA", "TC", "TG"]
MUT_RE = re.compile(r"^([ACGT])(\d+)([ACGT])$")


def parse_spectra(lines):
    """Return (classes, {node: {class: count}}, {node: tips}, order).

    Counts are rebuilt from proportion x Total_Mutations, so the folded
    columns are used when present: they are what a --fold search tested.
    """
    header = lines[0].rstrip("\n").split("\t")
    folded = [h for h in header if "|" in h]
    classes = folded if folded else [h for h in header if h in MUTATION_TYPES]
    col = {h: i for i, h in enumerate(header)}
    counts, tips, order = {}, {}, []
    for line in lines[1:]:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(classes) + 4 or not f[0]:
            continue
        total = int(f[col["Total_Mutations"]])
        counts[f[0]] = {c: float(f[col[c]]) * total for c in classes}
        tips[f[0]] = int(f[col["Number_Tips"]])
        order.append(f[0])
    return classes, counts, tips, order


def class_of(mut_type, classes):
    """The column a 12-type mutation counts toward."""
    if mut_type in classes:
        return mut_type
    for c in classes:
        if mut_type in c.split("|"):
            return c
    return None


def chi2_contributions(a, b, classes):
    """Pearson chi-squared of a 2 x K table, split into per-class shares."""
    ta, tb = sum(a[c] for c in classes), sum(b[c] for c in classes)
    n = ta + tb
    contrib = {}
    for c in classes:
        col = a[c] + b[c]
        if col <= 0 or ta <= 0 or tb <= 0:
            contrib[c] = 0.0
            continue
        x = 0.0
        for obs, tot in ((a[c], ta), (b[c], tb)):
            exp = tot * col / n
            x += (obs - exp) ** 2 / exp
        contrib[c] = x
    return contrib


def drivers(a, b, classes):
    """(chi2, ranked [(class, share, direction)], classes to reach 80%)."""
    contrib = chi2_contributions(a, b, classes)
    total = sum(contrib.values())
    ta, tb = sum(a.values()), sum(b.values())
    ranked = []
    for c in sorted(classes, key=lambda c: -contrib[c]):
        share = contrib[c] / total if total else 0.0
        up = ta and tb and a[c] / ta > b[c] / tb
        ranked.append((c, share, "+" if up else "-"))
    run, n80 = 0.0, 0
    for _, share, _ in ranked:
        run += share
        n80 += 1
        if run >= 0.8:
            break
    return total, ranked, n80


def residual_sites(tree, node_id, split_ids, max_branch_length):
    """Per-site class counts of a split's residual, as spectrumSplits counts.

    Walks the subtree, stopping at nested splits, and skips any branch with
    more than max_branch_length mutations, as compute_mutation_spectrum does.
    """
    sites = collections.defaultdict(collections.Counter)
    stack = [tree.get_node(node_id)]
    while stack:
        node = stack.pop()
        muts = node.mutations
        if len(muts) <= max_branch_length:
            for m in muts:
                mm = MUT_RE.match(m)
                if mm:
                    sites[int(mm.group(2))][mm.group(1) + mm.group(3)] += 1
        for child in node.children:
            if child.id not in split_ids:
                stack.append(child)
    return sites


def fold_sites(sites, classes):
    out = {}
    for pos, by_type in sites.items():
        c = collections.Counter()
        for t, k in by_type.items():
            cls = class_of(t, classes)
            if cls:
                c[cls] += k
        out[pos] = c
    return out


def site_loss(a, b, a_sites, b_sites, classes, top):
    """Largest fraction of chi-squared lost by removing one recurrent site.

    The site is removed from both rows, as masking it would.  Returns
    (loss, position, mutations at that site in the split).
    """
    base = sum(chi2_contributions(a, b, classes).values())
    if base <= 0:
        return None, None, None
    best = (None, None, None)
    ranked = sorted(a_sites, key=lambda p: -sum(a_sites[p].values()))[:top]
    for pos in ranked:
        a2 = {c: a[c] - a_sites[pos].get(c, 0) for c in classes}
        bs = b_sites.get(pos, {})
        b2 = {c: b[c] - bs.get(c, 0) for c in classes}
        loss = 1 - sum(chi2_contributions(a2, b2, classes).values()) / base
        if best[0] is None or loss > best[0]:
            best = (loss, pos, sum(a_sites[pos].values()))
    return best


def read_replicates(paths):
    """replicate name -> set of node IDs it split, one file per replicate."""
    reps = {}
    for path in paths:
        m = re.search(r"bootstrap_(\d+)_", os.path.basename(path))
        name = m.group(1) if m else path
        with open(path) as fh:
            reps[name] = {line.split("\t", 1)[0] for line in fh
                          if line.startswith("node_")}
    return reps


def verdict(top_share, loss, support, args):
    why = []
    if top_share >= args.one_class_share:
        why.append("one class")
    if loss is not None and loss >= args.one_site_loss:
        why.append("one site")
    if support is not None and support < args.min_support:
        why.append("weak support")
    return ", ".join(why) if why else "interpret"


def build_rows(args, spectra_lines, log_lines, reps=None, tree=None):
    classes, counts, tips, order = parse_spectra(spectra_lines)
    chi, parent, log_order = parse_splits_log(log_lines)
    split_ids = set(counts)
    sites = {}
    if tree is not None:
        for node in order:
            raw = residual_sites(tree, node, split_ids, args.max_branch_length)
            sites[node] = fold_sites(raw, classes)
            got = sum(sum(c.values()) for c in sites[node].values())
            want = sum(counts[node].values())
            if abs(got - want) > 0.5:
                print(f"warning: {node} residual has {got} mutations in the "
                      f"tree but {want:.0f} in spectra; is this the tree and "
                      f"--max_branch_length the search used?", file=sys.stderr)
    rows = []
    ordered = [n for n, _ in log_order if n in counts]
    ordered += [n for n in order if n not in set(ordered)]
    for node in ordered:
        par = parent.get(node)
        row = {"node": node, "parent": par or "", "chi_logged": chi.get(node),
               "mutations": round(sum(counts[node].values())),
               "tips": tips[node]}
        if par is None or par not in counts:
            row.update(chi2_residual=None, drivers="", top_share=None,
                       classes_80=None, site_loss=None, site=None,
                       site_mutations=None, support=None, verdict="root")
        else:
            total, ranked, n80 = drivers(counts[node], counts[par], classes)
            loss = pos = nsite = None
            if tree is not None:
                loss, pos, nsite = site_loss(counts[node], counts[par],
                                             sites[node], sites[par],
                                             classes, args.top_sites)
            support = None
            if reps:
                support = neighbourhood_support(
                    reps, node, args.neighbourhood) / len(reps)
            row.update(
                chi2_residual=total,
                drivers=",".join(f"{c}{d}:{s:.2f}" for c, s, d in ranked[:3]),
                top_share=ranked[0][1], classes_80=n80, site_loss=loss,
                site=pos, site_mutations=nsite, support=support,
                verdict=verdict(ranked[0][1], loss, support, args))
        row["branch_mutations"] = ""
        if tree is not None:
            row["branch_mutations"] = ",".join(tree.get_node(node).mutations)
        rows.append(row)
    return classes, rows


COLUMNS = ["node", "parent", "chi_logged", "chi2_residual", "mutations",
           "tips", "drivers", "top_share", "classes_80", "site", "site_mutations",
           "site_loss", "support", "verdict", "branch_mutations"]


def fmt(v, digits=2):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.{digits}f}"
    return str(v)


def write_tsv(rows, handle):
    handle.write("\t".join(COLUMNS) + "\n")
    for r in rows:
        handle.write("\t".join(fmt(r[c], 3) for c in COLUMNS) + "\n")


def report(classes, rows, args, out):
    have_tree = any(r["branch_mutations"] for r in rows)
    have_boot = any(r["support"] is not None for r in rows)
    print(f"SPLIT DRIVERS ({len(classes)} classes"
          f"{', folded' if '|' in classes[0] else ''}; each split's residual "
          f"against its parent's residual)", file=out)
    print(f"  one class: top class >= {args.one_class_share:.0%} of "
          f"chi2_residual", file=out)
    if have_tree:
        print(f"  one site:  removing one of the top {args.top_sites} sites "
              f"loses >= {args.one_site_loss:.0%} of it", file=out)
    else:
        print("  one site:  not tested (needs --input_tree)", file=out)
    if have_boot:
        print(f"  weak support: < {args.min_support:.0%} of replicates split "
              f"within +/-{args.neighbourhood} node numbers", file=out)
    else:
        print("  support:   not tested (needs --bootstrap)", file=out)
    print(file=out)
    head = ["node", "parent", "chi_logged", "mutations", "tips",
            "top 3 classes (share of chi2)", "site", "site_loss", "support",
            "verdict"]
    table = [head]
    for r in rows:
        table.append([r["node"], r["parent"], fmt(r["chi_logged"], 1),
                      str(r["mutations"]), str(r["tips"]), r["drivers"],
                      fmt(r["site"]), fmt(r["site_loss"]),
                      fmt(r["support"]), r["verdict"]])
    widths = [max(len(t[i]) for t in table) for i in range(len(head))]
    for t in table:
        print("  " + "  ".join(v.ljust(w) for v, w in zip(t, widths)).rstrip(),
              file=out)
    keep = [r for r in rows if r["verdict"] == "interpret"]
    counts = collections.Counter(r["verdict"] for r in rows)
    print(file=out)
    print("  " + ", ".join(f"{v}: {k}" for v, k in counts.most_common()),
          file=out)
    print(file=out)
    print(f"SPLITS TO INTERPRET ({len(keep)})", file=out)
    if not keep:
        return
    print("  " + " ".join(r["node"] for r in keep), file=out)
    if not have_tree:
        print("  (give --input_tree to list each one's defining-branch "
              "mutations)", file=out)
        return
    tsv = args.output or "split_drivers.tsv"
    v = COLUMNS.index("verdict") + 1
    b = COLUMNS.index("branch_mutations") + 1
    print(file=out)
    print("  Branch mutations per split (a branch longer than "
          "--max_branch_length was left", file=out)
    print("  out of the spectrum, so it labels the lineage and is skipped "
          "below):", file=out)
    print(file=out)
    for r in keep:
        n = len(r["branch_mutations"].split(",")) if r["branch_mutations"] \
            else 0
        flag = "  (longer than --max_branch_length: skipped)" \
            if n > args.max_branch_length else ""
        print(f"  {r['node']}: {n}{flag}", file=out)
    print(file=out)
    print("  To annotate every remaining split's defining branch:", file=out)
    print(file=out)
    loop = f"""\
  awk -F'\\t' 'NR>1 && ${v}=="interpret" && split(${b}, m, ",") <= {args.max_branch_length} {{print $1}}' {tsv} |
  while read n; do
    echo "== $n" >&2
    annotate_sites.py --annotation {args.annotation} --reference {args.reference} \\
        --highlight {args.highlight} --output ${{n}}_branch.tsv \\
        --positions "$(awk -F'\\t' -v n=$n '$1==n {{print $NF}}' {tsv})"
  done 2> branch_annotation.log
  awk '/^== /{{n=$2}} /HIGHLIGHTED/{{print n":"$0}}' branch_annotation.log"""
    print(loop, file=out)
    print(file=out)
    print("  Then test the hits against chance (the mask BEDs are the "
          "--mask_bed outputs of step 3):", file=out)
    print(file=out)
    print(f"  branch_hits.py --branches node_*_branch.tsv --annotation "
          f"{args.annotation} \\", file=out)
    print(f"      --highlight {args.highlight} --mask <mask BEDs> "
          f"--reference {args.reference}", file=out)

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Which classes drive each split, and which to interpret.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spectra", help="spectrumSplits.py --output_spectrum")
    p.add_argument("--splits_log", help="spectrumSplits.py stderr; gives each "
                   "split's parent (bootstrap lines are ignored)")
    p.add_argument("--bootstrap", nargs="+", default=None,
                   help="bootstrap_*_splits_output.tsv files")
    p.add_argument("--input_tree", default=None,
                   help="The tree the search ran on, for the one-site test "
                        "and the defining-branch mutations")
    p.add_argument("--max_branch_length", type=int, default=100000,
                   help="As given to spectrumSplits.py (its default 100000)")
    p.add_argument("--one_class_share", type=float, default=0.8)
    p.add_argument("--one_site_loss", type=float, default=0.5)
    p.add_argument("--top_sites", type=int, default=5)
    p.add_argument("--min_support", type=float, default=0.5)
    p.add_argument("--neighbourhood", type=int, default=10)
    p.add_argument("--annotation", default="NC_012345.1.gff3.gz",
                   help="Only used to print the annotate_sites.py command")
    p.add_argument("--reference", default="NC_012345.1.fna",
                   help="Only used to print the annotate_sites.py command")
    p.add_argument("--highlight", default="repair_genes.txt",
                   help="Only used to print the annotate_sites.py command")
    p.add_argument("--output", default=None, help="Per-split TSV")
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args(argv)
    if not args.selftest and not (args.spectra and args.splits_log):
        p.error("--spectra and --splits_log are required")
    return args


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    with open(args.spectra) as fh:
        spectra_lines = fh.readlines()
    with open(args.splits_log) as fh:
        log_lines = fh.readlines()
    reps = read_replicates(args.bootstrap) if args.bootstrap else None
    tree = None
    if args.input_tree:
        import bte
        print(f"Loading {args.input_tree}", file=sys.stderr)
        tree = bte.MATree(args.input_tree)
    classes, rows = build_rows(args, spectra_lines, log_lines, reps, tree)
    if args.output:
        with open(args.output, "w") as fh:
            write_tsv(rows, fh)
    report(classes, rows, args, sys.stdout)
    return 0


class FakeNode:
    def __init__(self, nid, muts, children=()):
        self.id, self.mutations, self.children = nid, list(muts), list(children)


class FakeTree:
    def __init__(self, root):
        self.nodes = {}
        stack = [root]
        while stack:
            n = stack.pop()
            self.nodes[n.id] = n
            stack.extend(n.children)

    def get_node(self, nid):
        return self.nodes[nid]


def selftest():
    ok = True

    def check(cond, label):
        nonlocal ok
        print(("ok   " if cond else "FAIL ") + label)
        ok = ok and cond

    ### Contributions sum to the textbook statistic.
    a = {"x": 30, "y": 70}
    b = {"x": 60, "y": 40}
    chi = sum(chi2_contributions(a, b, ["x", "y"]).values())
    check(abs(chi - 18.18) < 0.01, f"2x2 chi-squared {chi:.2f} == 18.18")
    _, ranked, _ = drivers(a, b, ["x", "y"])
    check(ranked[0][2] in "+-" and abs(ranked[0][1] + ranked[1][1] - 1) < 1e-9,
          "shares sum to 1")

    ### Folded mapping.
    classes = ["CA|GT", "CT|GA"]
    check(class_of("GT", classes) == "CA|GT", "GT folds into CA|GT")

    ### A tree where node_70's difference is one site's 40 C>T mutations.
    leaves_a = [FakeNode(f"a{i}", ["C10A", "G20T"]) for i in range(20)]
    hot = [FakeNode(f"h{i}", ["C500T"]) for i in range(40)]
    leaves_b = [FakeNode(f"b{i}", ["C11A", "G21T"]) for i in range(20)]
    node_b = FakeNode("node_70", ["A7G", "T9C"], leaves_b + hot)
    root = FakeNode("node_1", [], leaves_a + [node_b])
    tree = FakeTree(root)
    sites = residual_sites(tree, "node_1", {"node_1", "node_70"}, 100000)
    check(sum(sum(c.values()) for c in sites.values()) == 40,
          "residual stops at a nested split")
    sb = fold_sites(residual_sites(tree, "node_70", {"node_1", "node_70"},
                                   100000), ["CA|GT", "CT|GA", "TC|AG"])
    sr = fold_sites(sites, ["CA|GT", "CT|GA", "TC|AG"])
    cls = ["CA|GT", "CT|GA", "TC|AG"]
    cb = {c: sum(s.get(c, 0) for s in sb.values()) for c in cls}
    cr = {c: sum(s.get(c, 0) for s in sr.values()) for c in cls}
    loss, pos, n = site_loss(cb, cr, sb, sr, cls, 5)
    check(pos == 500 and n == 40 and loss > 0.9,
          f"one hot site found: pos {pos}, loss {loss:.2f}")
    check(residual_sites(tree, "node_70", {"node_70"}, 1)[500]["CT"] == 40,
          "max_branch_length drops long branches only")

    ### End to end on text inputs.
    spectra = [
        "Node_ID\tTotal_Mutations\tNumber_Tips\tMutations:Tips\t"
        "CA|GT\tCT|GA\tTC|AG\n",
        "node_1\t40\t20\t2\t1.0\t0\t0\n",
        "node_70\t82\t60\t1.4\t0.4878\t0.4878\t0.0244\n",
    ]
    log = ["Computing distances between splits in (sub)tree node_1\n",
           "New split found at node_70 with x2 99.0\n",
           "Bootstrapping splits with 2 replicates using 1 threads.\n",
           "New split found at node_70 with x2 1.0\n"]
    args = parse_args(["--selftest"])
    reps = {"1": {"node_1", "node_70"}, "2": {"node_1"}}
    _, rows = build_rows(args, spectra, log, reps, tree)
    by = {r["node"]: r for r in rows}
    check(by["node_1"]["verdict"] == "root", "root has no parent")
    rb = by["node_70"]
    check(rb["parent"] == "node_1" and rb["chi_logged"] == 99.0,
          "parent and chi read from before the bootstrap marker")
    check(abs(rb["support"] - 0.5) < 1e-9, f"support {rb['support']}")
    check("one site" in rb["verdict"] and "weak support" not in rb["verdict"],
          f"verdict '{rb['verdict']}'")
    check(rb["branch_mutations"] == "A7G,T9C", "defining-branch mutations")
    import io
    buf = io.StringIO()
    report(["CA|GT", "CT|GA", "TC|AG"], rows, args, buf)
    check("SPLITS TO INTERPRET (0)" in buf.getvalue(), "report renders")
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
