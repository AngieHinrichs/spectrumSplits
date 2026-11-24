# Automated qc tools
These scripts provide automated QC for reducing the impacts of potentially spurious mutations and samples. 

## prune_mutation_sample_ratio.py usage:
```usage: prune_mutation_sample_ratio.py [-h] [--input_tree INPUT_TREE] [--output_tree OUTPUT_TREE] [--threshold THRESHOLD]

Process a phylogenetic tree for changepoint detection in mutation/descendant ratios.

options:
  -h, --help            show this help message and exit
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
  --output_tree OUTPUT_TREE
                        Output tree file (protobuf format)
  --threshold THRESHOLD
                        Mutation:Leaf Ratio to prune
```

## mask_site_splits.py usage:
```usage: mask_site_splits.py [-h] [--input_tree INPUT_TREE] [--output_tree OUTPUT_TREE] [--min_total MIN_TOTAL] [--min_count MIN_COUNT] [--nthreads NTHREADS] [--mask_chi MASK_CHI]

Process a phylogenetic tree to mask highly imbalanced mutations.

options:
  -h, --help            show this help message and exit
  --input_tree INPUT_TREE
                        Input tree file (protobuf format)
  --output_tree OUTPUT_TREE
                        Output tree file (protobuf format)
  --min_total MIN_TOTAL
                        Minimum mutation count to accept a split
  --min_count MIN_COUNT
                        Minimum mutation count to accept a split
  --nthreads NTHREADS   Number of concurrent threads for processing
  --mask_chi MASK_CHI   Minimum chi2 value for masking mutation below a node (defaults to off)
```

## find_reverted_sites.py usage:
Finds the signature of reference backfilling (no-calls written as the
reference base instead of N): reversions to the reference that fall on tips
far more often than forward mutations do.  Prints a report with an overall
verdict, how concentrated the excess is across sites and samples, the top
sites and samples with a verdict for each, and what to do.  Needs `bte` and
`scipy`.  `--reference` must be the sequence the sample calls were made
against, not the root of a rerooted tree.

```
find_reverted_sites.py --input_tree tree.pb.gz --reference ref.fa \
    --output reverted_sites.tsv --sample_output reverted_samples.tsv \
    --pairs_output reverted_pairs.tsv > reverted_report.txt
```

`--output` ranks every site with a tip reversion by excess tip reversions;
`--sample_output` does the same per sample; `--pairs_output` lists the
sample, position and mutation of every tip reversion at a flagged site, for
branch-specific masking.  Options: `--top` (rows in the report, default 20),
`--alpha` (Bonferroni family-wise level, default 0.01), `--min_reversions`
(default 3), `--homopolymer_window` (default 15), `--chrom` (for a
multi-sequence FASTA).
