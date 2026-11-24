#!/usr/bin/env python3
"""Find sites and samples where an UShER tree has too many reversions to the
reference allele on its tips: the signature of reference backfilling.

When a pipeline writes the reference base where a sample has no confident
call (low depth, a homopolymer, a deletion) instead of N, the sample looks as
if it has reverted to the reference there.  Placed on a tree whose internal
nodes carry the non-reference allele, that "reversion" becomes a mutation on
the sample's own branch.  Backfilled calls recur in sample after sample at
the same hard-to-call sites, so they show up as highly recurrent mutations
that point back at the reference, almost all of them on tips.

The test: a reversion is a mutation whose to-allele is the reference base;
every other mutation is a forward mutation.  For real mutations, reversions
and forward mutations should fall on tips and internal branches in about the
same proportion.  Backfilled calls are private to one sample, so they land on
tips.  The tree's own forward mutations give the expected tip:internal ratio;
reversions on tips beyond that ratio are the excess.  The same comparison is
made per site (tip vs internal reversions at that site) and per sample (the
sample's share of reversions vs the internal-branch share).  Internal
branches can carry backfilled calls too, when a clade of samples shares the
no-call, so the excess is a lower bound.

Output:
  stdout            a report: overall verdict, how concentrated the excess
                    is, the top sites and samples with a verdict for each,
                    and what to do.
  --output          per-site TSV, every site with a tip reversion, ranked by
                    excess tip reversions.
  --sample_output   per-sample TSV, every sample with a reversion, ranked
                    the same way.
  --pairs_output    sample, position and mutation for every tip reversion at
                    a flagged site: the input for branch-specific masking.

The --reference FASTA must be the sequence the samples' calls were made
against (and backfilled from), not the root of a rerooted tree.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

import bte
from scipy.stats import binom

MUT_RE = re.compile(r'^([ACGT])(\d+)([ACGT])$')


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input_tree', required=True,
                   help='UShER MAT protobuf (.pb or .pb.gz)')
    p.add_argument('--reference', required=True,
                   help='FASTA of the reference the sample calls were made '
                        'against (plain or .gz)')
    p.add_argument('--chrom', default=None,
                   help='sequence name to use if --reference holds more than '
                        'one sequence')
    p.add_argument('--top', type=int, default=20,
                   help='sites and samples to show in the report (default 20)')
    p.add_argument('--alpha', type=float, default=0.01,
                   help='family-wise significance level for flagging a site '
                        'or sample, Bonferroni-corrected over all tested '
                        '(default 0.01)')
    p.add_argument('--min_reversions', type=int, default=3,
                   help='minimum tip reversions for a site or sample to be '
                        'flagged (default 3)')
    p.add_argument('--homopolymer_window', type=int, default=15,
                   help='report the longest homopolymer within this many '
                        'bases of each site (default 15)')
    p.add_argument('--output', default=None, help='per-site TSV')
    p.add_argument('--sample_output', default=None, help='per-sample TSV')
    p.add_argument('--pairs_output', default=None,
                   help='TSV of sample, position, mutation for each tip '
                        'reversion at a flagged site')
    return p.parse_args()


def read_reference(path, chrom):
    opener = gzip.open if path.endswith('.gz') else open
    seqs, name, buf = {}, None, []
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    if not seqs:
        sys.exit(f'No sequences in {path}')
    if chrom:
        if chrom not in seqs:
            sys.exit(f'--chrom {chrom} not in {path} (has: {", ".join(seqs)})')
        return chrom, seqs[chrom].upper()
    if len(seqs) > 1:
        sys.exit(f'{path} has {len(seqs)} sequences; choose one with --chrom')
    name = next(iter(seqs))
    return name, seqs[name].upper()


def longest_homopolymer(ref, pos, window):
    """Longest single-base run within [pos-window, pos+window] (1-based)."""
    lo, hi = max(0, pos - 1 - window), min(len(ref), pos + window)
    best = run = 0
    prev = ''
    for b in ref[lo:hi]:
        run = run + 1 if b == prev else 1
        prev = b
        best = max(best, run)
    return best


def homopolymer_background(ref, min_run, window):
    """Fraction of positions within `window` of a run of >= min_run."""
    cover = bytearray(len(ref))
    i, n = 0, len(ref)
    while i < n:
        j = i
        while j < n and ref[j] == ref[i]:
            j += 1
        if j - i >= min_run:
            cover[max(0, i - window):min(n, j + window)] = \
                b'\x01' * (min(n, j + window) - max(0, i - window))
        i = j
    return sum(cover) / n


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    if n <= 0 or k > n:
        return 0.0
    return float(binom.sf(k - 1, n, p))


def pct(a, b):
    return f'{100.0 * a / b:.1f}%' if b else 'n/a'


class Counts:
    __slots__ = ('tip_rev', 'tip_fwd', 'int_rev', 'int_fwd', 'types')

    def __init__(self):
        self.tip_rev = self.tip_fwd = self.int_rev = self.int_fwd = 0
        self.types = defaultdict(int)

    @property
    def tip(self):
        return self.tip_rev + self.tip_fwd

    @property
    def total(self):
        return self.tip_rev + self.tip_fwd + self.int_rev + self.int_fwd


def main():
    args = parse_args()
    chrom, ref = read_reference(args.reference, args.chrom)
    tree = bte.MATree(args.input_tree)

    sites = defaultdict(Counts)
    samples = defaultdict(Counts)
    tip_rev_events = defaultdict(list)   # pos -> [(sample, mutation)]
    tot = Counts()
    skipped = out_of_range = 0

    for node in tree.depth_first_expansion():
        leaf = node.is_leaf()
        for m in node.mutations:
            mm = MUT_RE.match(m)
            if not mm:
                skipped += 1
                continue
            pos = int(mm.group(2))
            if pos < 1 or pos > len(ref):
                out_of_range += 1
                continue
            r = ref[pos - 1]
            if r not in 'ACGT':
                skipped += 1
                continue
            rev = mm.group(3) == r
            s = sites[pos]
            s.types[m[0] + '>' + m[-1]] += 1
            targets = (s, tot, samples[node.id]) if leaf else (s, tot)
            for c in targets:
                if leaf:
                    if rev:
                        c.tip_rev += 1
                    else:
                        c.tip_fwd += 1
                elif rev:
                    c.int_rev += 1
                else:
                    c.int_fwd += 1
            if leaf and rev:
                tip_rev_events[pos].append((node.id, m))

    if tot.int_fwd == 0 or tot.int_rev + tot.int_fwd == 0 or tot.tip == 0:
        sys.exit('Too few mutations to estimate a baseline.')

    # Expected tip:internal ratio, from forward mutations.
    ratio = tot.tip_fwd / tot.int_fwd
    p_tip = ratio / (1 + ratio)
    exp_tip_rev = tot.int_rev * ratio
    excess_total = tot.tip_rev - exp_tip_rev
    # Expected share of a sample's mutations that are reversions: the
    # internal-branch share.
    int_rev_share = tot.int_rev / (tot.int_rev + tot.int_fwd)

    site_rows = []
    tested = sum(1 for c in sites.values() if c.tip_rev + c.int_rev > 0)
    for pos, c in sites.items():
        if c.tip_rev == 0:
            continue
        exp = c.int_rev * ratio
        exc = c.tip_rev - exp
        pval = binom_sf(c.tip_rev, c.tip_rev + c.int_rev, p_tip)
        flag = (c.tip_rev >= args.min_reversions
                and pval * max(tested, 1) < args.alpha)
        site_rows.append((pos, c, exc, exp, pval, flag))
    site_rows.sort(key=lambda r: (-r[2], r[0]))
    flagged_sites = [r for r in site_rows if r[5]]

    sample_rows = []
    tested_s = sum(1 for c in samples.values() if c.tip > 0)
    for sid, c in samples.items():
        if c.tip_rev == 0:
            continue
        exp = c.tip * int_rev_share
        exc = c.tip_rev - exp
        pval = binom_sf(c.tip_rev, c.tip, int_rev_share)
        flag = (c.tip_rev >= args.min_reversions
                and pval * max(tested_s, 1) < args.alpha)
        sample_rows.append((sid, c, exc, exp, pval, flag))
    sample_rows.sort(key=lambda r: (-r[2], r[0]))
    flagged_samples = [r for r in sample_rows if r[5]]

    hp = {r[0]: longest_homopolymer(ref, r[0], args.homopolymer_window)
          for r in site_rows}

    # ---------------- report ----------------
    w = sys.stdout.write
    w('REFERENCE-REVERSION CHECK\n')
    w(f'  tree: {args.input_tree}\n  reference: {args.reference} ({chrom}, '
      f'{len(ref):,} bp)\n')
    w(f'  {tree.count_leaves():,} tips; {tot.tip:,} tip mutations, '
      f'{tot.int_rev + tot.int_fwd:,} internal-branch mutations\n')
    if skipped or out_of_range:
        w(f'  skipped {skipped:,} non-ACGT mutations and {out_of_range:,} '
          f'outside the reference\n')
    w('\n')

    w('OVERALL\n')
    w('  A reversion is a mutation to the reference base; every other\n'
      '  mutation is forward.  Real reversions should be split between tips\n'
      '  and internal branches like forward mutations are.\n\n')
    w(f'  forward mutations:  {tot.tip_fwd:>10,} on tips, {tot.int_fwd:>10,} '
      f'internal  (tip:internal = {ratio:.2f})\n')
    w(f'  reversions:         {tot.tip_rev:>10,} on tips, {tot.int_rev:>10,} '
      f'internal  (tip:internal = '
      f'{tot.tip_rev / tot.int_rev if tot.int_rev else float("inf"):.2f})\n')
    w(f'  expected tip reversions at the forward ratio: {exp_tip_rev:,.0f}\n')
    w(f'  excess tip reversions: {excess_total:,.0f} '
      f'({pct(max(excess_total, 0), tot.tip)} of all tip mutations)\n\n')

    exc_share = max(excess_total, 0) / tot.tip
    if exc_share >= 0.10:
        verdict = ('STRONG SIGN OF REFERENCE BACKFILLING.  More than a tenth '
                   'of all tip mutations are excess reversions to the '
                   'reference.  Check how no-calls are encoded upstream (they '
                   'should be N, not the reference base).')
    elif exc_share >= 0.02:
        verdict = ('SOME SIGN OF REFERENCE BACKFILLING.  A few percent of tip '
                   'mutations are excess reversions.  Check the pipeline at a '
                   'few flagged sites, and mask or prune what is flagged.')
    else:
        verdict = ('NO CLEAR SIGN OF REFERENCE BACKFILLING.  Excess tip '
                   'reversions are under 2% of tip mutations.  Any flagged '
                   'sites below are individual problem sites.')
    w('  VERDICT: ' + verdict + '\n')
    w('  (Strong: excess >= 10% of tip mutations; some: >= 2%.  The excess is\n'
      '  a lower bound, since internal branches can carry shared backfill.)\n\n')

    if excess_total > 0 and site_rows:
        w('HOW CONCENTRATED IS IT?\n')
        pos_exc = [max(r[2], 0) for r in site_rows]
        tot_pos = sum(pos_exc) or 1
        for k in (20, 100, 1000, 10000):
            if len(pos_exc) >= k:
                w(f'  top {k:>6,} sites hold {pct(sum(pos_exc[:k]), tot_pos)} '
                  f'of the excess\n')
        w(f'  flagged sites: {len(flagged_sites):,}, holding '
          f'{pct(sum(max(r[2], 0) for r in flagged_sites), tot_pos)} of the '
          f'excess\n')
        spos = [max(r[2], 0) for r in sample_rows]
        stot = sum(spos) or 1
        for k in (10, 100, 1000):
            if len(spos) >= k:
                w(f'  worst {k:>5,} samples hold {pct(sum(spos[:k]), stot)} '
                  f'of the excess\n')
        w(f'  flagged samples: {len(flagged_samples):,} of '
          f'{tree.count_leaves():,}\n')
        if exc_share >= 0.02:
            if sum(pos_exc[:20]) < 0.5 * tot_pos:
                w('  -> Spread over many sites: masking the top sites will not '
                  'remove most\n     of it.\n')
            top_tenth = max(1, len(spos) // 10)
            if sum(spos[:top_tenth]) < 0.5 * stot:
                w('  -> Spread over most samples, not a few bad ones: pruning '
                  'samples will\n     not remove most of it.  This points at '
                  'how the calls were made.\n')
        if flagged_sites:
            near = sum(1 for r in flagged_sites if hp[r[0]] >= 8)
            bg = homopolymer_background(ref, 8, args.homopolymer_window)
            w(f'  flagged sites within {args.homopolymer_window} bp of a run '
              f'of 8+ identical bases: {near:,} of {len(flagged_sites):,} '
              f'({pct(near, len(flagged_sites))}), against {100 * bg:.1f}% '
              f'of all positions\n')
            if len(flagged_sites) >= 10 and near / len(flagged_sites) > 2 * bg:
                w('  -> Homopolymers are over-represented: the no-calls are '
                  'likely where\n     reads lose confidence.\n')
        w('\n')

    w(f'TOP {min(args.top, len(site_rows))} SITES BY EXCESS TIP REVERSIONS\n')
    w('  pos         ref  tip_rev  int_rev  expected  rev_share  homopol  '
      'main_changes              verdict\n')
    for pos, c, exc, exp, pval, flag in site_rows[:args.top]:
        main = ', '.join(f'{k} {v}' for k, v in
                         sorted(c.types.items(), key=lambda kv: -kv[1])[:2])
        share = c.tip_rev / c.total
        if flag and share >= 0.5:
            v = 'BACKFILL-LIKE'
        elif flag:
            v = 'excess, mixed'
        else:
            v = 'not significant'
        w(f'  {pos:<11} {ref[pos - 1]:<4} {c.tip_rev:>7,}  {c.int_rev:>7,}  '
          f'{exp:>8.1f}  {100 * share:>8.1f}%  {hp[pos]:>7}  {main:<24}  {v}\n')
    w('\n  tip_rev, int_rev: reversions on tips and internal branches.  '
      'expected: tip\n  reversions at the forward tip:internal ratio.  '
      'rev_share: tip reversions as\n  a share of all mutations at the site.  '
      'homopol: longest single-base run\n  within '
      f'{args.homopolymer_window} bp.  BACKFILL-LIKE: significant excess, and '
      'most changes at the\n  site are tip reversions.  excess, mixed: '
      'significant, but other changes\n  recur there too.\n\n')

    w(f'TOP {min(args.top, len(sample_rows))} SAMPLES BY EXCESS TIP REVERSIONS\n')
    w(f'  (expected share of reversions: {100 * int_rev_share:.1f}%, the '
      f'internal-branch share)\n')
    w('  sample                      mutations  reversions  rev_share  verdict\n')
    for sid, c, exc, exp, pval, flag in sample_rows[:args.top]:
        share = c.tip_rev / c.tip
        if flag and share >= 0.8:
            v = 'MOSTLY REVERSIONS: likely low coverage; consider pruning'
        elif flag:
            v = 'excess reversions'
        else:
            v = 'not significant'
        w(f'  {sid:<27} {c.tip:>9,}  {c.tip_rev:>10,}  {100 * share:>8.1f}%  '
          f'{v}\n')
    w('\n')

    w('WHAT TO DO\n')
    w('  - Look at a few top sites in the per-sample calls (VCF/GVCF) of '
      'samples\n    that revert there (--pairs_output lists them).  If those '
      'samples have\n    no confident call at the site, the reference base was '
      'backfilled.\n'
      '  - Fix upstream if you can: encode no-calls as N so placement imputes\n'
      '    them, and rebuild.\n'
      '  - Otherwise mask branch-specifically: the --pairs_output rows are the\n'
      '    sample-site pairs to mask, which keeps real mutations at those '
      'sites on\n    other branches.  Masking whole sites removes those too.\n'
      '  - Samples flagged MOSTLY REVERSIONS are candidates for pruning.\n')

    # ---------------- files ----------------
    if args.output:
        with open(args.output, 'w') as f:
            f.write('position\tref\ttip_reversions\tinternal_reversions\t'
                    'expected_tip_reversions\texcess\tp_value\tflagged\t'
                    'tip_forward\tinternal_forward\treversion_share\t'
                    'longest_homopolymer\tchanges\n')
            for pos, c, exc, exp, pval, flag in site_rows:
                ch = ','.join(f'{k}:{v}' for k, v in
                              sorted(c.types.items(), key=lambda kv: -kv[1]))
                f.write(f'{pos}\t{ref[pos - 1]}\t{c.tip_rev}\t{c.int_rev}\t'
                        f'{exp:.2f}\t{exc:.2f}\t{pval:.3g}\t{int(flag)}\t'
                        f'{c.tip_fwd}\t{c.int_fwd}\t'
                        f'{c.tip_rev / c.total:.3f}\t{hp[pos]}\t{ch}\n')
    if args.sample_output:
        with open(args.sample_output, 'w') as f:
            f.write('sample\ttip_mutations\treversions\texpected\texcess\t'
                    'p_value\tflagged\treversion_share\n')
            for sid, c, exc, exp, pval, flag in sample_rows:
                f.write(f'{sid}\t{c.tip}\t{c.tip_rev}\t{exp:.2f}\t{exc:.2f}\t'
                        f'{pval:.3g}\t{int(flag)}\t{c.tip_rev / c.tip:.3f}\n')
    if args.pairs_output:
        with open(args.pairs_output, 'w') as f:
            f.write('sample\tposition\tmutation\n')
            for pos, *_ in flagged_sites:
                for sid, m in sorted(tip_rev_events[pos]):
                    f.write(f'{sid}\t{pos}\t{m}\n')


if __name__ == '__main__':
    main()
