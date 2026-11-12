# Parameter tuning helpers

Two read-only scripts for choosing `spectrumSplits.py` parameters on a tree
other than the one the defaults were calibrated on. Neither modifies the tree.

Both need only `bte` and the Python standard library.

## Why the defaults do not transfer

`spectrumSplits.py` scores each candidate node with a 2 x 12 chi-squared test:
the candidate subtree's mutation spectrum against the rest of the subtree it
sits in. That statistic grows roughly as

```
chi-squared  ~  (size of the spectral difference)  x  (effective mutation count)
```

where the effective count is dominated by the smaller of the two sides. So the
same biological effect gives a chi-squared ten times larger on a tree with ten
times the mutations. `--min_chi 500` is a working point on the SARS-CoV-2 tree,
not a significance cutoff — with 11 degrees of freedom it corresponds to
p ~ 1e-100. It is high because the algorithm takes the maximum over every node
in a subtree (massive multiple testing) and because recurrent sites violate the
test's independence assumption.

That means `--min_chi` has to be calibrated, not copied.

## 1. `profile_tree.py` — measure the tree first

```
profile_tree.py --input_tree tree.pb.gz
profile_tree.py --input_tree tree.pb.gz --node_table nodes.tsv
```

Reports:

| Section | What to do with it |
| --- | --- |
| Tree size, total usable SNVs | Sets the ceiling on how much splitting is possible |
| Global 12-type and strand-folded 6-type spectrum | Shows whether folding to 6 categories is worth it, and gives the rarest-class floor on `--min_mutations` (5 / p_rarest, so expected cell counts stay at or above 5) |
| Per-branch mutation counts, and what each `--max_branch_length` cap would exclude | Pick a cap only if there is a clear outlier tail; on a diverged, sparsely sampled tree long branches are usually real |
| Subtree mutation counts, and eligible candidate nodes at each `--min_mutations` | If only a handful of nodes clear your threshold, the analysis will be coarse whatever `--min_chi` is |
| Mutations:tips ratio, with the `prune_mutation_sample_ratio.py` auto-threshold | Shows what that QC script would prune before you run it |
| Site recurrence (mutations per site, top sites' share) | The independence violation that forces `--min_chi` above any nominal cutoff; `qc/mask_site_splits.py` targets it |

`--node_table` writes one row per node (branch and subtree counts, tips,
mutations:tips, per-branch 12-type spectrum) if you want to plot distributions.

## 2. `permutation_null.py` — calibrate `--min_chi`

Destroys any association between mutation type and position in the tree, then
runs the same search and the same statistic, and reports the chi-squared values
it still finds. Set `--min_chi` comfortably above the maximum.

```
# headline calibration
permutation_null.py --input_tree tree.pb.gz --min_mutations 500 \
    --replicates 20 --mode sites --include_observed
```

Two null models:

- `--mode types` shuffles mutation-type labels across all mutations, keeping
  topology and per-branch counts fixed. The strict independence null: it also
  removes the type clumping that recurrent sites cause, so it gives the lower,
  more permissive threshold.
- `--mode sites` (default) gives each site a single type (its most common one)
  and shuffles that assignment across sites. Which branches carry mutations at
  which site is untouched, so a site mutating 500 times inside one clade still
  dumps 500 counts of one randomized type into that clade. This keeps the
  independence violation and gives the higher, more conservative threshold.
  Prefer it, especially with convergent selection at resistance loci.

Key options:

- `--min_mutations` / `--max_branch_length` should match what you intend for the
  real run; they change which nodes are eligible and so change the null.
- `--max_splits` defaults to `1`: one pass over the whole tree per replicate,
  which is exactly the quantity the first real split has to beat. Raise it to
  follow the greedy cascade into subtrees, where chi-squared is necessarily
  smaller and not comparable to the first pass.
- `--fold` scores the strand-folded 6-type spectrum (6 categories, 5 degrees of
  freedom) instead of 12. Worth considering for a double-stranded genome with
  no strand-asymmetric replication process.
- `--include_observed` scores the unpermuted tree through the same code path, so
  the real maximum and the null maximum are directly comparable. If the ratio is
  under about 2x, the strongest real split is close to what chance produces.
- `--nthreads` runs replicates concurrently where `fork` is available, and falls
  back to serial otherwise; running several instances with different `--seed`
  also works.
- `--selftest` checks the chi-squared implementation against a naive reference
  and against `scipy.stats.chi2_contingency`, and exits. Worth running once in
  your environment to confirm the statistic matches `spectrumSplits.py`.

The statistic is reimplemented in pure Python rather than called through scipy,
so all candidate nodes can be scored in one vectorless pass without scipy's
per-call overhead. It matches `chi2_contingency` exactly for tables larger than
2 x 2 (no continuity correction); columns whose total is zero are skipped, where
scipy would raise.

## Suggested order

1. Run `profile_tree.py` on the raw tree.
2. Run the `qc/` scripts, inspecting what each one removes.
3. Run `profile_tree.py` again on the cleaned tree; pick `--min_mutations` from
   the rarest-class floor and the precision you need, sanity-checked against how
   many nodes stay eligible.
4. Run `permutation_null.py` with that `--min_mutations`; set `--min_chi` above
   the worst null replicate, with a margin.
5. Run `spectrumSplits.py` once with a permissive `--min_chi`, logging every
   split's chi-squared from stderr, and check where the splits-versus-threshold
   curve bends.
6. Confirm the surviving splits with `--bootstrap_splits`.
7. Inspect each split's per-class drivers and effect size in proportion terms
   before believing it. A split driven almost entirely by one mutation class is
   usually an artifact or a single over-mutated site.
