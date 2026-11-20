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
     takes.

A split is listed to interpret when it is multi-class, not one site, and
supported.  Each test is skipped, not failed, when its input is not given.

With --annotation, --reference and --highlight too, it goes on to look for
a mechanism, as step 8's manual loop and branch_hits.py used to:

  5. GENE NAME CHECK.  check_gene_names.py on the --highlight list: a name
     the annotation spells differently would drop out of the highlight
     without anything looking wrong.

  6. DEFINING-BRANCH ANNOTATION.  Each interpret split's defining branch,
     no longer than --max_branch_length (a longer one was left out of the
     spectrum, so it labels the lineage), annotated by annotate_sites.py's
     own code and written as node_N_branch.tsv in --branch_dir; the
     highlighted-gene hits on each are listed.

  7. BRANCH HITS.  Those hits against chance.  A branch with hundreds of
     mutations hits a gene list covering 1-2% of the genome whatever the
     mechanism.  The rate per mutation is the genes' length outside the
     --mask BEDs over the genome's, times the branches' nonsynonymous share
     of coding mutations; only protein-changing hits count unless --count
     all.  Per branch, per gene on one branch, and per gene across
     branches, each with a "by chance" count: how many of the tests made
     would reach that P by chance.  Read that first.

Outputs: PREFIX.tsv (--prefix, default split_drivers) per split, and
PREFIX_hits.tsv per annotated branch.

Usage:
  split_drivers.py --spectra spectra.tsv --splits_log splits.log \\
      --bootstrap bootstrap_*_splits_output.tsv \\
      --input_tree tree.pruned.masked.pb.gz --max_branch_length 250 \\
      --annotation genes.gff3.gz --reference ref.fa \\
      --highlight repair_genes.txt --mask selection_loci.bed,alignment_loci.bed \\
      > split_drivers.txt

  split_drivers.py --selftest
"""

import argparse
import collections
import math
import os
import re
import sys
from collections import defaultdict

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
    print(file=out)
    print("  Branch mutations per split (a branch longer than "
          "--max_branch_length was left", file=out)
    print("  out of the spectrum, so it labels the lineage and is not "
          "annotated):", file=out)
    print(file=out)
    for r in keep:
        n = len(r["branch_mutations"].split(",")) if r["branch_mutations"] \
            else 0
        flag = "  (longer than --max_branch_length: not annotated)" \
            if n > args.max_branch_length else ""
        print(f"  {r['node']}: {n}{flag}", file=out)


# --------------------------------------------------------------------------
# defining-branch annotation and hits against chance
# (formerly the annotate_sites.py loop and branch_hits.py)
# --------------------------------------------------------------------------

def poisson_tail(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    term = math.exp(-lam)
    below = term
    for i in range(1, k):
        term *= lam / i
        below += term
    return max(0.0, 1.0 - below)


def chance_count(p, lams):
    """Expected number of tests, of Poisson means lams, reaching P <= p.

    Each test's counts are discrete, so its chance of reaching P <= p is the
    tail at the smallest count whose tail is <= p, not p itself.
    """
    total = 0.0
    for lam in lams:
        k = 0
        while poisson_tail(k, lam) > p:
            k += 1
            if k > 10000:
                break
        total += poisson_tail(k, lam)
    return total


def overlap(intervals, mask):
    """Bases of merged 1-based intervals that fall inside merged mask."""
    n = 0
    for s, e in intervals:
        for ms, me in mask:
            if me < s:
                continue
            if ms > e:
                break
            n += min(e, me) - max(s, ms) + 1
    return n


def length(intervals):
    return sum(e - s + 1 for s, e in intervals)


def gene_intervals(features, names):
    """name -> merged intervals, plus the names the annotation lacks."""
    from annotate_sites import (LANDMARK_TYPES, merge_intervals,
                                resolve_region_names)
    features = [f for f in features if f.ftype not in LANDMARK_TYPES]
    matched, missing = resolve_region_names(features, names)
    out = {}
    for name, feats in matched.items():
        out[name] = merge_intervals([(f.start, f.end) for f in feats])
    return out, missing


def summarise_branch(rows):
    """One branch's annotate_site() rows as branch-hit counts.

    hits: (gene, mutation, aa_change, nonsynonymous) per site INSIDE a
    highlighted gene; a site only near one (name@+N) is not a hit.
    """
    n = coding = nonsyn = 0
    hits = []
    for row in rows:
        n += 1
        aa = str(row["aa_change"])
        is_nonsyn = "(nonsyn)" in aa
        if is_nonsyn or "(syn)" in aa:
            coding += 1
            nonsyn += is_nonsyn
        hl = row.get("highlight", "")
        if hl and "@" not in hl:
            mut = row["mutation"] or str(row["position"])
            hits.append((hl, mut, aa.split(" ")[0], is_nonsyn))
    return {"mutations": n, "coding": coding, "nonsyn": nonsyn, "hits": hits}


def analyse(branches, gene_ivs, mask, genome_size, count="nonsyn"):
    """All the numbers the report prints, as a dict."""
    from annotate_sites import merge_intervals
    masked = length([(s, min(e, genome_size)) for s, e in mask
                     if s <= genome_size])
    searchable = genome_size - masked
    union = merge_intervals([iv for ivs in gene_ivs.values() for iv in ivs])
    hl_bp = length(union)
    hl_free = hl_bp - overlap(union, mask)
    gene_free = {g: length(ivs) - overlap(ivs, mask)
                 for g, ivs in gene_ivs.items()}

    all_hits = sum(len(d["hits"]) for d in branches.values())
    hit_nonsyn = sum(h[3] for d in branches.values() for h in d["hits"])
    coding = sum(d["coding"] for d in branches.values())
    nonsyn = sum(d["nonsyn"] for d in branches.values())
    ### P(a mutation is counted | it lands in a highlighted gene).
    frac = 1.0 if count == "all" else (nonsyn / coding if coding else 0.0)
    rate = hl_free / searchable * frac

    counted = {b: [h for h in d["hits"] if count == "all" or h[3]]
               for b, d in branches.items()}
    lams = {b: d["mutations"] * rate for b, d in branches.items()}
    rows = []
    for b, d in branches.items():
        k = len(counted[b])
        genes = defaultdict(list)
        for h in counted[b]:
            genes[h[0]].append(h)
        rows.append({"branch": b, "mutations": d["mutations"], "hits": k,
                     "expected": lams[b], "p": poisson_tail(k, lams[b]),
                     "genes": dict(genes)})
    for r in rows:
        r["chance"] = chance_count(r["p"], lams.values()) if r["hits"] \
            else None
    rows.sort(key=lambda r: (r["p"], -r["hits"]))

    n_mut = sum(d["mutations"] for d in branches.values())
    n_counted = sum(len(v) for v in counted.values())
    tot = {"mutations": n_mut, "hits": n_counted, "expected": n_mut * rate,
           "all_hits": all_hits, "hit_nonsyn": hit_nonsyn, "coding": coding,
           "nonsyn": nonsyn}
    tot["p"] = poisson_tail(n_counted, tot["expected"])

    ### One gene, one branch.
    pair_lams = {(b, g): d["mutations"] * frac * bp / searchable
                 for b, d in branches.items() for g, bp in gene_free.items()}
    cluster = []
    for r in rows:
        for g, hs in r["genes"].items():
            if len(hs) < 2 or g not in gene_free:
                continue
            lam = pair_lams[(r["branch"], g)]
            p = poisson_tail(len(hs), lam)
            cluster.append({"branch": r["branch"], "gene": g, "hits": hs,
                            "expected": lam, "p": p,
                            "chance": chance_count(p, pair_lams.values())})
    cluster.sort(key=lambda r: r["p"])

    ### One gene, several branches.
    gene_lams = {g: sum(1 - math.exp(-d["mutations"] * frac * bp / searchable)
                        for d in branches.values())
                 for g, bp in gene_free.items()}
    recur = []
    for g in gene_free:
        hit_by = [(r["branch"], r["genes"][g]) for r in rows if g in r["genes"]]
        if len(hit_by) < 2:
            continue
        p = poisson_tail(len(hit_by), gene_lams[g])
        recur.append({"gene": g, "bp": gene_free[g], "hit_by": hit_by,
                      "expected": gene_lams[g], "p": p,
                      "chance": chance_count(p, gene_lams.values())})
    recur.sort(key=lambda r: r["p"])
    return {"genome": genome_size, "masked": masked, "searchable": searchable,
            "hl_bp": hl_bp, "hl_free": hl_free, "frac": frac, "rate": rate,
            "count": count, "rows": rows, "total": tot, "cluster": cluster,
            "recur": recur, "n_genes": len(gene_free)}


def pct(a, b):
    return f"{100 * a / b:.0f}%" if b else "n/a"


def hit_label(h):
    return f"{h[0]} {h[2]}" if h[2] else f"{h[0]} {h[1]}"


def print_table(table, out):
    widths = [max(len(t[i]) for t in table) for i in range(len(table[0]))]
    for t in table:
        print("  " + "  ".join(v.ljust(w) for v, w in zip(t, widths)).rstrip(),
              file=out)


def hits_report(res, names_file, n_names, missing, out):
    t = res["total"]
    what = "nonsynonymous hits" if res["count"] == "nonsyn" else "hits"
    print(f"BRANCH HITS: genes from {names_file}; counting {what}", file=out)
    print(f"  {res['n_genes']} of {n_names} names found in the annotation"
          + (f"; not found: {', '.join(missing)}" if missing else ""), file=out)
    print(f"  highlighted genes: {res['hl_bp']:,} bp, {res['hl_free']:,} "
          f"outside the mask", file=out)
    print(f"  genome: {res['genome']:,} bp, {res['masked']:,} masked, "
          f"{res['searchable']:,} searchable", file=out)
    line = (f"  chance rate per mutation: {res['hl_free']:,} / "
            f"{res['searchable']:,}")
    if res["count"] == "nonsyn":
        line += (f" x {t['nonsyn']:,}/{t['coding']:,} nonsynonymous "
                 f"coding mutations")
    print(f"{line} = {100 * res['rate']:.3f}%", file=out)
    print(file=out)
    print("  'by chance' is how many of the tests made would reach that P by "
          "chance. Read it", file=out)
    print("  first: near or above 1, the result is what the best of that many "
          "draws looks", file=out)
    print("  like, and its rank means nothing.", file=out)
    print(file=out)

    print("ALL BRANCHES", file=out)
    print(f"  {len(res['rows'])} branches, {t['mutations']:,} mutations, "
          f"{t['hits']} {what}, {t['expected']:.1f} expected, "
          f"P(>={t['hits']}) = {t['p']:.3g}", file=out)
    print(f"  nonsynonymous share: {t['hit_nonsyn']} of {t['all_hits']} hits "
          f"({pct(t['hit_nonsyn'], t['all_hits'])}), against {t['nonsyn']:,} "
          f"of {t['coding']:,} ({pct(t['nonsyn'], t['coding'])}) coding "
          f"mutations on these branches", file=out)
    print(file=out)

    n = len(res["rows"])
    print(f"PER BRANCH ('by chance' out of {n} branches)", file=out)
    table = [["branch", "by chance", "P", what.split()[0] if
              res["count"] == "nonsyn" else "hits", "expected", "mutations",
              "genes"]]
    for r in res["rows"]:
        if not r["hits"]:
            continue
        genes = ", ".join(hit_label(h) for hs in r["genes"].values()
                          for h in hs)
        table.append([r["branch"], f"{r['chance']:.2f}", f"{r['p']:.3g}",
                      str(r["hits"]), f"{r['expected']:.2f}",
                      str(r["mutations"]), genes])
    if len(table) > 1:
        print_table(table, out)
    else:
        print("  none", file=out)
    print(file=out)

    pairs = n * res["n_genes"]
    print(f"ONE GENE, ONE BRANCH ('by chance' out of {pairs:,} branch-gene "
          f"pairs)", file=out)
    if res["cluster"]:
        table = [["branch", "gene", "by chance", "P", "hits", "expected",
                  "mutations"]]
        for r in res["cluster"]:
            table.append([r["branch"], r["gene"], f"{r['chance']:.3g}",
                          f"{r['p']:.3g}", str(len(r["hits"])),
                          f"{r['expected']:.3f}",
                          ", ".join(f"{h[1]} {h[2]}" for h in r["hits"])])
        print_table(table, out)
    else:
        print("  none", file=out)
    print(file=out)

    print(f"ONE GENE, SEVERAL BRANCHES ('by chance' out of {res['n_genes']} "
          f"genes)", file=out)
    if res["recur"]:
        table = [["gene", "by chance", "P", "branches", "expected", "bp",
                  "hits"]]
        for r in res["recur"]:
            hits = "; ".join(f"{b} " + ", ".join(h[2] or h[1] for h in hs)
                             for b, hs in r["hit_by"])
            table.append([r["gene"], f"{r['chance']:.2f}", f"{r['p']:.3g}",
                          str(len(r["hit_by"])), f"{r['expected']:.2f}",
                          f"{r['bp']:,}", hits])
        print_table(table, out)
    else:
        print("  none", file=out)


def write_hits_tsv(res, path):
    with open(path, "w") as fh:
        fh.write("branch\tmutations\thits\texpected\tp\tby_chance\tgenes\n")
        for r in res["rows"]:
            vals = [r["branch"], str(r["mutations"]), str(r["hits"]),
                    f"{r['expected']:.4g}", f"{r['p']:.4g}",
                    "" if r["chance"] is None else f"{r['chance']:.4g}",
                    ",".join(f"{h[0]}:{h[1]}" for hs in r["genes"].values()
                             for h in hs)]
            fh.write("\t".join(vals) + "\n")



BRANCH_HEADER = ["position", "count", "mutation", "context", "feature_type",
                 "gene", "locus_tag", "strand", "codon", "codon_pos",
                 "ref_codon", "ref_aa", "aa_change", "distance", "product",
                 "highlight"]


def load_reference(args):
    """(features without landmarks, sequence) from --annotation and
    --reference, choosing the sequence as annotate_sites.py does."""
    from annotate_sites import LANDMARK_TYPES, load_annotation, load_fasta
    features = [f for f in load_annotation(args.annotation, args.seqid)
                if f.ftype not in LANDMARK_TYPES]
    seqs = load_fasta(args.reference)
    if args.seqid and args.seqid in seqs:
        seq = seqs[args.seqid]
    else:
        seq = max(seqs.values(), key=len) if seqs else ""
    return features, seq


def annotate_branches(rows, args, features, seq, highlight_names):
    """annotate_site() rows for every interpret split's defining branch no
    longer than --max_branch_length, as {node: [row, ...]}, written to
    --branch_dir/node_N_branch.tsv when that is set.  Also returns the
    splits skipped for length."""
    from annotate_sites import FeatureIndex, annotate_site
    index = FeatureIndex(features)
    done, skipped = {}, []
    for r in rows:
        if r["verdict"] != "interpret" or not r["branch_mutations"]:
            continue
        muts = r["branch_mutations"].split(",")
        if len(muts) > args.max_branch_length:
            skipped.append((r["node"], len(muts)))
            continue
        out = []
        for m in muts:
            mm = MUT_RE.match(m)
            if not mm:
                continue
            row, _kind, _label, _hl = annotate_site(
                int(mm.group(2)), None, m.upper(), index, seq, BRANCH_HEADER,
                args.flank, highlight_names, 0)
            out.append(row)
        done[r["node"]] = out
        if args.branch_dir:
            path = os.path.join(args.branch_dir, f"{r['node']}_branch.tsv")
            with open(path, "w") as fh:
                fh.write("\t".join(BRANCH_HEADER) + "\n")
                for row in out:
                    fh.write("\t".join(str(row[k]) for k in BRANCH_HEADER)
                             + "\n")
    return done, skipped


def name_check(names, features, args, out):
    """check_gene_names.py on the --highlight list: a name the annotation
    spells differently drops out of the highlight without anything looking
    wrong.  Returns True if every entry resolves to one locus."""
    from check_gene_names import build_loci, check
    from check_gene_names import report as names_report
    loci, by_name = build_loci(features)
    rows = check(names, loci, by_name, 2)
    shim = argparse.Namespace(names=args.highlight, annotation=args.annotation,
                              max_candidates=5)
    names_report(rows, shim, out)
    return not any(r["status"] in ("several", "not found") for r in rows)


def branch_report(done, skipped, out):
    print(f"DEFINING-BRANCH ANNOTATION ({len(done)} interpret splits)",
          file=out)
    if skipped:
        print("  not annotated, longer than --max_branch_length: "
              + ", ".join(f"{n} ({k})" for n, k in skipped), file=out)
    hits = [(n, [r for r in rows if r["highlight"]
                 and "@" not in r["highlight"]])
            for n, rows in done.items()]
    hits = [(n, h) for n, h in hits if h]
    print(f"  {len(hits)} branches hit a highlighted gene:", file=out)
    for n, h in hits:
        print(f"    {n}: " + ", ".join(
            f"{r['highlight']} {str(r['aa_change']) or r['mutation']}"
            for r in h), file=out)
    print(file=out)


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
    g = p.add_argument_group(
        "defining branches (needs --input_tree, --annotation, --reference "
        "and --highlight)")
    g.add_argument("--annotation", default=None,
                   help="GFF3 or GenBank annotation, as for annotate_sites.py")
    g.add_argument("--reference", default=None,
                   help="Reference FASTA: codons, and the genome length")
    g.add_argument("--highlight", default=None,
                   help="Gene list whose hits on the interpret splits' "
                        "defining branches are tested against chance")
    g.add_argument("--mask", default="",
                   help="Comma-separated BED files of the masked regions "
                        "(step 3's --mask_bed outputs): masked sites cannot "
                        "carry a mutation, so they leave the chance rate")
    g.add_argument("--seqid", default=None)
    g.add_argument("--count", choices=("nonsyn", "all"), default="nonsyn",
                   help="Count only protein-changing hits (default), or all")
    g.add_argument("--flank", type=int, default=200,
                   help="Promoter window, as in annotate_sites.py "
                        "(default 200)")
    g.add_argument("--branch_dir", default=".",
                   help="Where to write node_N_branch.tsv per annotated "
                        "branch (default the current directory; '' skips)")
    p.add_argument("--prefix", default="split_drivers",
                   help="Prefix of the TSVs: PREFIX.tsv per split, and "
                        "PREFIX_hits.tsv per annotated branch (default "
                        "split_drivers)")
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
    with open(f"{args.prefix}.tsv", "w") as fh:
        write_tsv(rows, fh)
    out = sys.stdout
    report(classes, rows, args, out)
    print(file=out)

    need = [o for o in ("input_tree", "annotation", "reference", "highlight")
            if not getattr(args, o)]
    if need:
        print("  To annotate the interpret splits' defining branches and test "
              "their highlighted-", file=out)
        print("  gene hits against chance, also give "
              + ", ".join(f"--{o}" for o in need) + ".", file=out)
        return 0
    from check_gene_names import read_names
    features, seq = load_reference(args)
    names = read_names(args.highlight)
    highlight_names = {}
    for n in names:
        highlight_names.setdefault(n.lower(), n)
    name_check(names, features, args, out)
    print(file=out)
    done, skipped = annotate_branches(rows, args, features, seq,
                                      highlight_names)
    branch_report(done, skipped, out)
    from annotate_sites import read_bed_intervals
    gene_ivs, missing = gene_intervals(features, names)
    mask = read_bed_intervals([m for m in args.mask.split(",") if m]) \
        if args.mask else []
    ### By name, as the shell glob over node_N_branch.tsv gave them, so
    ### ties sort the same way every run.
    branches = {n: summarise_branch(done[n]) for n in sorted(done)}
    res = analyse(branches, gene_ivs, mask, len(seq), args.count)
    write_hits_tsv(res, f"{args.prefix}_hits.tsv")
    hits_report(res, args.highlight, len(names), missing, out)
    print(f"Wrote {args.prefix}.tsv, {args.prefix}_hits.tsv"
          + (f", and {len(done)} node_N_branch.tsv in {args.branch_dir}"
             if args.branch_dir else ""), file=sys.stderr)
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

    _selftest_hits(check)
    _selftest_branches(check)
    print("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


def _selftest_branches(check):
    """Annotate two constructed branches from a small GFF3 and FASTA."""
    import io
    import tempfile
    d = tempfile.mkdtemp()
    gff = os.path.join(d, "g.gff3")
    with open(gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for gid, s0, e0, name in (("g1", 1, 30, "recA"),
                                  ("g2", 41, 70, "other")):
            fh.write(f"chr\tx\tgene\t{s0}\t{e0}\t.\t+\t.\tID={gid};"
                     f"gene={name};locus_tag=L{gid}\n")
            fh.write(f"chr\tx\tCDS\t{s0}\t{e0}\t.\t+\t0\tID=c{gid};"
                     f"Parent={gid};gene={name};locus_tag=L{gid}\n")
    fa = os.path.join(d, "r.fa")
    with open(fa, "w") as fh:
        fh.write(">chr\n" + "ATG" + "AAA" * 9 + "C" * 10 + "ATG" + "AAA" * 9
                 + "C" * 30 + "\n")
    args = parse_args(["--selftest"])
    args.annotation, args.reference, args.highlight = gff, fa, "genes.txt"
    args.max_branch_length, args.branch_dir = 5, d
    features, seq = load_reference(args)
    check(len(seq) == 100, f"reference read: {len(seq)} bp")
    rows = [{"node": "node_70", "verdict": "interpret",
             "branch_mutations": "A1G,A44C"},
            {"node": "node_9", "verdict": "interpret",
             "branch_mutations": ",".join(f"A{p}G" for p in range(1, 9))},
            {"node": "node_5", "verdict": "one class",
             "branch_mutations": "A1G"}]
    done, skipped = annotate_branches(rows, args, features, seq,
                                      {"reca": "recA"})
    check(list(done) == ["node_70"] and skipped == [("node_9", 8)],
          "only interpret splits within --max_branch_length are annotated")
    b = summarise_branch(done["node_70"])
    check(b["hits"] == [("recA", "A1G", "M1V", True)] and b["coding"] == 2
          and b["nonsyn"] == 2, f"hit in recA, M1V: {b}")
    check(os.path.exists(os.path.join(d, "node_70_branch.tsv")),
          "node_70_branch.tsv written")
    buf = io.StringIO()
    ok_names = name_check(["recA", "nucS"], features, args, buf)
    check(not ok_names and "nucS" in buf.getvalue(),
          "a highlight name the annotation lacks is reported")


def _selftest_hits(check):
    check(abs(poisson_tail(8, 3.16) - 0.0158) < 0.001,
          f"P(X>=8 | 3.16) = {poisson_tail(8, 3.16):.4f}")
    check(poisson_tail(0, 2.0) == 1.0, "P(X>=0) is 1")
    check(chance_count(0.01, [0.0001]) < 0.01, "a short branch adds ~0")
    check(abs(chance_count(0.5, [5.0]) - poisson_tail(6, 5.0)) < 1e-12,
          "discrete tail used, not p")
    check(overlap([(10, 20)], [(15, 30)]) == 6, "overlap is inclusive")

    import io
    hdr = ("position\tcount\tmutation\tcontext\tfeature_type\tgene\t"
           "locus_tag\tstrand\tcodon\tcodon_pos\tref_codon\tref_aa\t"
           "aa_change\tdistance\tproduct\thighlight\n")
    rows = ["1\t\tA1G\tCDS\tCDS\tx\tR1\t+\t1\t1\tATG\tM\tM1V (nonsyn)\t0\tp\trecA\n",
            "2\t\tA2G\tCDS\tCDS\tx\tR1\t+\t1\t2\tATG\tM\tM1I (nonsyn)\t0\tp\trecA\n",
            "3\t\tA3G\tCDS\tCDS\tx\tR1\t+\t1\t3\tATG\tM\tM1M (syn)\t0\tp\trecA\n",
            "4\t\tA4G\tCDS\tCDS\ty\tR2\t+\t2\t1\tATG\tM\tM2V (nonsyn)\t0\tp\t\n",
            "5\t\tA5G\tpromoter\t\t\t\t\t\t\t\t\t\t-5\t\tmutY@-5\n"]
    cols = hdr.rstrip("\n").split("\t")
    b = summarise_branch([dict(zip(cols, r.rstrip("\n").split("\t")))
                          for r in rows])
    check(b["mutations"] == 5 and len(b["hits"]) == 3 and b["coding"] == 4
          and b["nonsyn"] == 3, "hits inside the gene only; flank not counted")

    gene_ivs = {"recA": [(1, 100)], "mutY": [(201, 300)]}
    syn_only = {"mutations": 10, "coding": 10, "nonsyn": 5,
                "hits": [("mutY", "C250T", "L5L", False)]}
    branches = {"node_7": b, "node_9": {"mutations": 10, "coding": 10,
                                        "nonsyn": 5,
                                        "hits": [("recA", "A9G", "K3E", True),
                                                 ("mutY", "C251T", "L6L",
                                                  False)]},
                "node_80": syn_only}
    res = analyse(branches, gene_ivs, [(51, 100), (901, 1000)], 1000)
    check(res["searchable"] == 850 and res["hl_free"] == 150,
          f"mask removed from both: {res['hl_free']}/{res['searchable']}")
    check(abs(res["frac"] - 13 / 24) < 1e-12, "nonsyn fraction pooled")
    check(abs(res["rate"] - 150 / 850 * 13 / 24) < 1e-12, "rate scaled")
    check(res["total"]["hits"] == 3, "synonymous hits not counted")
    check([r["gene"] for r in res["recur"]] == ["recA"],
          "a synonymous-only gene does not recur")
    check([(r["branch"], r["gene"], len(r["hits"])) for r in res["cluster"]]
          == [("node_7", "recA", 2)], "clustering counts nonsyn only")
    res_all = analyse(branches, gene_ivs, [(51, 100), (901, 1000)], 1000,
                      "all")
    check(res_all["total"]["hits"] == 6 and res_all["frac"] == 1.0
          and {r["gene"] for r in res_all["recur"]} == {"recA", "mutY"},
          "--count all counts every hit, unscaled")
    buf = io.StringIO()
    hits_report(res, "genes.txt", 2, [], buf)
    check("ONE GENE, SEVERAL BRANCHES" in buf.getvalue(), "report renders")



if __name__ == "__main__":
    sys.exit(main())
