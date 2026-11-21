# Tuning spectrumSplits on a new tree

`spectrumSplits.py` finds clades whose mutation spectrum differs from the
rest of a phylogeny. Its parameters were set on the SARS-CoV-2 tree and do
not transfer to another tree. Its splits also need checking before they can
be read as biology: against chance, against sequencing quality, against
what each split contains, and against lineage. This directory holds the
scripts that do that tuning and checking, and the documentation for using
them. None of the scripts modifies a tree.

The workflow was developed on a *Mycobacterium tuberculosis* complex tree of
129,314 genomes, and every decision in it is illustrated with that tree's
numbers. Applying it to another species means reading the same reports on
your own tree.

## Where to start

| Document | Read it to |
| --- | --- |
| [HOWTO.md](HOWTO.md) | Run the workflow on a new tree, step by step: the commands, which part of each report to read, and what to decide from it |
| [METHODS.md](METHODS.md) | Understand why each step works as it does: the statistics, the caveats, the full MTBC worked examples, the record of how the MTBC tree was cleaned in four rounds, and the tests of the workflow itself |
| [SCRIPTS.md](SCRIPTS.md) | Look up a script's options, inputs and output columns, or the gene lists |

Each step in HOWTO.md links to its section in METHODS.md, and back.

## The workflow

| Step | What it does | Main script |
| --- | --- | --- |
| 1 | Profile the raw tree, and choose a divergence threshold | [`profile_tree.py`](SCRIPTS.md#profile_treepy) |
| 2 | Prune genomes with implausible divergence, and laboratory strains | [`prune_clade_check.py`](SCRIPTS.md#prune_clade_checkpy) |
| 3 | Mask regions under selection or prone to misalignment | [`annotate_sites.py`](SCRIPTS.md#annotate_sitespy), [`check_gene_names.py`](SCRIPTS.md#check_gene_namespy) |
| 4 | Re-profile, and decide whether the mask needs another round | [`profile_tree.py`](SCRIPTS.md#profile_treepy), [`annotate_sites.py`](SCRIPTS.md#annotate_sitespy) |
| 5 | Decide whether to fold strand-complementary classes, then profile trinucleotide contexts | [`profile_tree.py`](SCRIPTS.md#profile_treepy) |
| 6 | Choose `--min_mutations` and `--max_branch_length` | [`profile_tree.py`](SCRIPTS.md#profile_treepy) |
| 7 | Calibrate `--min_chi` against permutation nulls, by subtree size, and check for a hotspot | [`permutation_null.py`](SCRIPTS.md#permutation_nullpy) |
| 8 | Run the search with bootstrap; screen splits for one-class and one-site drivers; test hits on their defining branches against chance | [`split_drivers.py`](SCRIPTS.md#split_driverspy) |
| 9 | Rule out a sequencing-quality gradient | [`branch_class_spectra.py`](SCRIPTS.md#branch_class_spectrapy) |
| 10 | Rank splits by size-matched margin and bootstrap support | [`split_confidence.py`](SCRIPTS.md#split_confidencepy) |
| 11 | Check what each split contains, and test it for a trinucleotide-context shift | [`split_composition.py`](SCRIPTS.md#split_compositionpy), [`split_contexts.py`](SCRIPTS.md#split_contextspy), [`plot_split_contexts.py`](SCRIPTS.md#plot_split_contextspy) |
| 12 | Test candidate genes genotype-first, against their siblings and a control list | [`gene_hit_contrasts.py`](SCRIPTS.md#gene_hit_contrastspy), [`control_genes.py`](SCRIPTS.md#control_genespy), [`plot_gene_hits.py`](SCRIPTS.md#plot_gene_hitspy) |
| 13 | Ask whether the splits' spectra change how often, or with which variants, they acquire resistance mutations | [`resistance_spectra.py`](SCRIPTS.md#resistance_spectrapy) |

The order matters more than any single setting. Masking changes the
spectrum, the spectrum sets the floor under `--min_mutations`, and both
reshape the null that `--min_chi` is calibrated against. After a change to
the tree, [`match_splits.py`](SCRIPTS.md#match_splitspy) matches the new
splits to the old by their samples
([HOWTO.md](HOWTO.md#after-a-change-to-the-tree-replicate-and-match)).

## Requirements

`matUtils` and `spectrumSplits.py`, and Python with the standard library and
`bte`. Some scripts need less:

- `prune_clade_check.py` and `split_drivers.py` use `bte` only with
  `--input_tree`.
- `check_gene_names.py`, `control_genes.py`, `split_confidence.py`,
  `split_composition.py` and `match_splits.py` read text files only.

Some need more:

- `split_contexts.py` and `gene_hit_contrasts.py --context_genes` need numpy.
- `resistance_spectra.py` needs numpy, and `openpyxl` to read the WHO
  catalogue's `.xlsx`.
- `plot_gene_hits.py` and `plot_split_contexts.py` draw figures, and need
  matplotlib rather than `bte`.

All but `profile_tree.py`, `prune_clade_check.py` and `plot_gene_hits.py`
have `--selftest`.
