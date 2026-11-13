# Spectrum splits

```
usage: spectrumSplits.py [-h] [--input_tree INPUT_TREE] [--output_spectrum OUTPUT_SPECTRUM] [--min_chi MIN_CHI] [--min_mutations MIN_MUTATIONS] [--ntips NTIPS]
                         [--bootstrap_splits BOOTSTRAP_SPLITS] [--bootstrap_spectra BOOTSTRAP_SPECTRA] [--nthreads NTHREADS] [--max_branch_length MAX_BRANCH_LENGTH]
                         [--fold]

Process a phylogenetic tree to find splits, compute spectra, and get representative tips.

options:
  -h, --help            show this help message and exit
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
  --output_spectrum OUTPUT_SPECTRUM
                        Output TSV file for spectra
  --min_chi MIN_CHI     Minimum Chi-square value to accept a split
  --min_mutations MIN_MUTATIONS
                        Minimum number of mutations required for a split
  --ntips NTIPS         Number of tips to retrieve for each split
  --bootstrap_splits BOOTSTRAP_SPLITS
                        Number of bootstrap replicates to attempt in defining splits
  --bootstrap_spectra BOOTSTRAP_SPECTRA
                        Number of bootstrap replicates to attempt in defining spectra
  --nthreads NTHREADS   Number of threads for concurrent bootstrapping
  --max_branch_length MAX_BRANCH_LENGTH
                        Maximum branch length to include in spectrum calculations
  --fold                Collapse the twelve mutation classes onto the six strand-complementary pairs (C>T with G>A, and so on) for the chi-squared test.
                        Appropriate for a double-stranded genome with no strand-asymmetric replication process; gives 5 degrees of freedom instead of 11 and
                        doubles the count in every cell. Changes the scale of the statistic, so --min_chi must be recalibrated.
```

## Strand folding (`--fold`)

By default each candidate split is scored on a 2 x 12 contingency table: the
twelve substitution classes of the candidate subtree against the same twelve
classes for the rest of the current split. That is the right table when the two
strands can differ, as they do for SARS-CoV-2.

For a double-stranded genome with no strand-asymmetric replication process, a
C>T on the reference strand and a G>A on the other strand are the same event,
and the twelve classes carry only six classes' worth of information. `--fold`
collapses them onto the six strand-complementary pairs before the test:

```
TG|AC   TC|AG   TA|AT   CA|GT   CG|GC   CT|GA
```

What this buys, confirmed by simulation under the null:

* **Degrees of freedom drop from 11 to 5.** A chi-squared statistic averages its
  degrees of freedom under the null, so the noise floor falls with it
  (simulated null means: 11.29 unfolded, 4.99 folded).
* **Every cell count doubles**, which matters for rare classes. The expected
  cell count floor on `--min_mutations` falls accordingly.
* **A strand-symmetric difference is preserved in full** -- folding loses none
  of that signal while halving the degrees of freedom, so the same difference
  becomes more significant.
* **Zero-count columns become far less likely**, which matters because the test
  cannot be computed for a class absent from the whole subtree.

What it costs: a genuinely strand-*asymmetric* difference loses half its
statistic, since the two halves of a pair now cancel. If you expect or want to
detect strand asymmetry, do not fold.

Decide empirically rather than by assumption. `tuning/profile_tree.py` prints
the global spectrum as both twelve classes and six folded pairs; if every
complementary pair agrees closely, the twelve-class test is spending degrees of
freedom on nothing.

`--fold` affects the chi-squared test and adds the six folded proportions to the
output TSV. The twelve unfolded columns are kept, so strand symmetry stays
visible per split -- a clade with a real strand-asymmetric process would show up
there and nowhere else.

**`--min_chi` does not carry over between folded and unfolded runs.** The
statistic changes scale and the null changes shape, so a threshold calibrated
one way is meaningless the other. Recalibrate with
`tuning/permutation_null.py --fold`.
