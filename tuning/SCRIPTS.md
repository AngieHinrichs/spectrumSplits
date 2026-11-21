# Scripts and gene lists

Usage, inputs, outputs and options for each script, and the gene lists they read. [HOWTO.md](HOWTO.md) says when to run each one; [METHODS.md](METHODS.md) says why.

## Script reference

| Script | Step | What it answers |
| --- | --- | --- |
| [`profile_tree.py`](#profile_treepy) | 1, 4, 5 | What is in this tree, and what do its numbers imply for `--min_mutations`, `--max_branch_length` and whether to fold? |
| [`prune_clade_check.py`](#prune_clade_checkpy) | 2 | What did `qc/prune_mutation_sample_ratio.py` remove, by clade, and was it a divergent genome or a clade on a long stem? |
| [`annotate_sites.py`](#annotate_sitespy) | 3, 4 | What are the recurrent sites, and which regions should be masked before anything else is trusted? |
| [`check_gene_names.py`](#check_gene_namespy) | 3, 8 | Does every name in a `--regions` or `--highlight` list resolve to exactly one gene in this annotation, and if not, what is it called there? |
| [`permutation_null.py`](#permutation_nullpy) | 7 | How large a chi-squared does chance alone reach on this tree, and so where should `--min_chi` sit? How does that ceiling fall with subtree size, which splits does a fixed `--min_chi` wrongly reject, and is one hotspot setting the threshold? |
| [`split_drivers.py`](#split_driverspy) | 8 | Which mutation classes drive each split, does one site or weak bootstrap support explain it, which splits are worth interpreting, and are the highlighted-gene hits on their defining branches more than chance? |
| [`branch_class_spectra.py`](#branch_class_spectrapy) | 9 | Is a split a real mutational difference, or a sequencing-quality gradient on terminal branches? |
| [`split_confidence.py`](#split_confidencepy) | 10 | Which accepted splits are worth interpreting? Scores each against a size-matched null and the bootstrap replicates. |
| [`split_composition.py`](#split_compositionpy) | 11 | What is each split made of: one lineage, or several averaged together? |
| [`split_contexts.py`](#split_contextspy) | 11 | Does a split also differ in trinucleotide context inside its classes, by more than any clade of its size would, on all branches and on internal branches alone? |
| [`plot_split_contexts.py`](#plot_split_contextspy) | 11 | Figures of `split_contexts.py`'s results: which splits stand out from random clades, the context profiles of those that pass, and the context shifts of every split side by side. |
| [`gene_hit_contrasts.py`](#gene_hit_contrastspy) | 12 | Does a candidate gene's mutation shift the spectrum below it? Tests every independent hit against its siblings, so lineage cannot stand in for mechanism, with every hit and with internal hits only, against a control list, and reports the control, the excess, the genes that clear the bar and the leads. |
| [`control_genes.py`](#control_genespy) | 12 | Is a gene list's excess its own, or what any genes of that length would show? Draws a length-matched random control list. |
| [`plot_gene_hits.py`](#plot_gene_hitspy) | 12 | A figure: one gene's shift against every gene in its list, and against its own permutation null. |
| [`resistance_spectra.py`](#resistance_spectrapy) | 13 | Does a split's spectrum change how often, or with which variants, it acquires resistance mutations, against the rest of the tree, by a graded catalogue? |
| [`plot_resistance_spectra.py`](#plot_resistance_spectrapy) | 13 | Figures of `resistance_spectra.py`'s results: predicted against observed resistance supply per split, llr against what the spectrum alone predicts, and the shares behind one split's llr. |
| [`match_splits.py`](#match_splitspy) | [testing](METHODS.md#testing-the-method) | Which splits of a second search on a changed tree are the same clades as the first's? |
| [`trinucleotide.py`](#trinucleotidepy) | 5, 11, 12 | Not run directly: calls each mutation's trinucleotide context from its parent's genome, for `profile_tree.py`, `split_contexts.py` and `gene_hit_contrasts.py --context`. |

The step numbers are [HOWTO.md](HOWTO.md)'s. The entries below are in
the order the steps use them.

### `profile_tree.py`

Used in [step 1](HOWTO.md#step-1-profile-the-raw-tree), again in
[step 4](HOWTO.md#step-4-re-profile-and-check-what-the-mask-bought), and for the context
spectrum at the end of [step 5](HOWTO.md#step-5-decide-whether-to-fold).

```
profile_tree.py --input_tree tree.pb.gz
profile_tree.py --input_tree tree.pb.gz --node_table nodes.tsv
profile_tree.py --input_tree tree.pb.gz --reference ref.fna --annotation genes.gff3.gz
profile_tree.py --input_tree tree.pb.gz --reference ref.fna --fold \
    --mask masked.bed --context_table context.tsv
```

Reports:

| Section | What to do with it |
| --- | --- |
| Tree size, total usable SNVs | Sets the ceiling on how much splitting is possible |
| Global 12-type and strand-folded 6-type spectrum | Gives the rarest-class floor on `--min_mutations` (5 / p_rarest, so expected cell counts stay at or above 5). Don't use it to decide whether to fold: strand asymmetries cancel in it (see [step 5](HOWTO.md#step-5-decide-whether-to-fold)) |
| Per-branch mutation counts, and what each `--max_branch_length` cap would exclude | Pick a cap only if there is a clear outlier tail; on a diverged, sparsely sampled tree long branches are usually real |
| Subtree mutation counts, and eligible candidate nodes at each `--min_mutations` | If only a handful of nodes clear your threshold, the analysis will be coarse whatever `--min_chi` is |
| Mutations:tips ratio, with the `prune_mutation_sample_ratio.py` auto-threshold | Shows what that QC script would prune before you run it |
| Site recurrence (mutations per site, top sites' share) | The independence violation that forces `--min_chi` above any nominal cutoff; `qc/mask_site_splits.py` targets it |
| Context spectrum (with `--reference`) | The trinucleotide context classes, 96 folded with `--fold` and 192 forward-strand without, with flanks from each mutation's parent genome: how many could be called, a check that the FASTA is the root's genome, count, share, opportunity and relative rate per context, any context dominated by one site, and the count floors a context test needs. [Step 5](HOWTO.md#step-5-decide-whether-to-fold) reads it, once the fold decision is made, and a note follows the FOLD TEST if `--fold` disagrees with its verdict |
| Strand asymmetry (with `--reference`, or `--terminus`) | Complementary-pair log2 ratios on the leading and lagging arms, in + and − strand genes (with `--annotation`), and per annotated clade, with a heterogeneity chi-squared for each pair. It ends with a FOLD TEST and a verdict. This is what step 5 reads |

The strand asymmetry section needs the origin and terminus of replication.
`--reference` estimates them from cumulative GC skew: the minimum is the origin
and the maximum the terminus, which assumes the leading strand is G-rich, as
it is in most bacteria. Without a reference, give `--terminus`; `--origin`
defaults to 1. Either can override the estimate. In the per-clade table each
mutation goes to its nearest annotated ancestor, so a parent clade's row
leaves out its annotated subclades. A clade is listed only if its rarest
pair's log2 ratio would have a standard error of at most `--strand_max_se`
(default 0.15); the report prints the clade size that implies for the tree's
spectrum, and `--strand_min_clade` replaces it with a fixed size. The FOLD
TEST uses the same clades. Its limits are `--fold_max_dispersion` (default 3)
and `--fold_max_z` (default 5); step 5 explains both.

The context spectrum needs the reference to be the root's genome (the
rerooted reference on a rerooted tree). `--mask` takes comma-separated BED
files: a mutation with a masked flank gets an unknown context, and masked
trinucleotides are left out of the opportunity counts. Pass the BEDs the tree
was masked with. `relative rate` divides each context's count by how often its
trinucleotide occurs, relative to its class's average (1 is average).
`--fold` picks the class set, as in the other scripts: the 96 folded classes
with it, the 192 forward-strand ones without, whose opportunity counts are
over the 64 forward-strand trinucleotides. It changes nothing else in the
report. `--context_table` writes the rows as a TSV, and `--no_context` skips
the section, as step 4 does before the fold decision.

`--node_table` writes one row per node (branch and subtree counts, tips,
mutations:tips, per-branch 12-type spectrum) if you want to plot distributions.

### `prune_clade_check.py`

Used in [step 2](HOWTO.md#step-2-prune-samples-with-implausible-divergence). Tabulates
what `qc/prune_mutation_sample_ratio.py` removed, by clade, and whether each
removal was a single divergent genome or a whole subtree. It reads only text
files; `bte` is needed only for `--input_tree`.

```
prune_clade_check.py --samples samples_raw.txt --clades clades_raw.txt \
    --prune_report pruned_report.txt --threshold 250 \
    [--exclude ReconstructedAncestor] [--input_tree tree.pb.gz] > prune_check.md
```

Inputs, all describing the **unpruned** tree:

- `--samples`: `matUtils summary -s` output. Its `parsimony` column is each
  sample's terminal branch length, and `parent_id` places it in the tree.
- `--clades`: `matUtils summary --sample-clades` output.
- `--prune_report`: the prune script's stdout, which lists every removed
  sample, so the pruned tree is not needed.
- `--threshold`: the value given to the prune script. It sets only the
  `tips >=` column.

Output is markdown, in five sections: removal rate by clade; terminal branches
at or above the threshold against removals, split into terminal-branch and
internal-node crossings; each internal-node pruning with its subtree total and
the tips' terminal branches; the surviving parent of each pruned subtree, with
the tips under it and their commonest clades (needs `--input_tree`); and the
unlabelled samples, as contiguous runs in node-number order with the labelled
lineage on either side.

- `--depth` collapses clade names to that many dot-separated parts (default 1,
  so `lineage4.3.4.2` counts as `lineage4`; 0 keeps full names).
- `--input_tree` reads the unpruned tree with `bte` to report each
  internal-node pruning's exact stem, and to fill in the parent section.
  Without it, the stem and the internal branches below it are reported
  together as one remainder.
- `--min_run` is the smallest unlabelled run listed (default 5).
- `--exclude` names samples to ignore, comma-separated (default none). On
  the MTBC tree that is `ReconstructedAncestor`, a placeholder sample put on
  the root when the tree was rerooted (`matUtils mask --rename-samples`).
  Other trees will have their own such names, or none.

### `annotate_sites.py`

Used in [step 3](HOWTO.md#step-3-build-and-apply-a-region-mask) to build the mask, and in
[step 4](HOWTO.md#step-4-re-profile-and-check-what-the-mask-bought) to identify what survived it.

`profile_tree.py` tells you *which* sites mutate repeatedly; this tells you
*what they are*, so you can decide what to mask. It does two jobs.

**Identify sites.** Reads positions straight out of a profile report:

```
annotate_sites.py --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --from_profile profile.txt --output annotated_sites.tsv
```

Reports feature type, gene, locus tag, strand, codon number, codon position,
reference codon and amino acid for each site. Given a bare position it lists all
three possible amino-acid changes; given a mutation string such as `C1057G` it
reports the specific change and whether it is synonymous. Reverse-strand genes
are handled (codon numbering runs from the feature's `end` downwards and bases
are reverse-complemented) — `--selftest` checks that arithmetic on both strands,
and it is worth running once.

**Build mask regions.** Masking a top-N list of individual sites does not
converge: mask `rpoB` codon 445 position 1 and position 2 of the same codon
surfaces in the next profile. Mask the locus instead:

```
annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions tb_resistance_genes.txt --promoter_flank 200 \
    --mask_bed resistance_loci.bed
```

`--promoter_flank` extends each feature upstream in a strand-aware way, which
matters because several important mutations are promoter mutations
(`inhA`/`fabG1` −15, `embA` −16, `pncA` −11). `--mask_pattern '^(PE|PPE)'` adds
features by name, for the repetitive families where alignment artifacts
accumulate. Output is a merged BED, or `--mask_positions` for a flat 1-based
position list (expect tens of thousands of lines).

`--mask_matutils` writes those same positions as `N<position>N`, one per line,
which is the format `matUtils mask -m` reads, so the mask can go straight into
`matUtils mask` with no conversion step.

`--mask_product_pattern` masks by **function** instead of by name, matching
each feature's product description: `'(?i)(transposase|integrase|insertion
element)'` finds mobile elements whatever they happen to be called. This is not
a convenience — `--mask_pattern` matches gene names and locus tags only, so a
transposase annotated as `Rv2648` is unreachable by name, and the genes most
worth masking for alignment reasons are exactly the ones with uninformative
names. Matching is case sensitive, like `--mask_pattern`, so prefix `(?i)`;
only CDS and RNA features carry products, so bare `gene` features never match
(their CDS does, over essentially the same interval). The run prints the loci
it matched.

`--near_mask` takes the BED files a mask was built from and reports, for
each site within `--near_mask_within` bp (default 200) of a masked interval,
its distance outside it (`mask_gap`: positive past the interval's right end,
negative before its left end) and the masked genes (`mask_interval`). The
columns are added only when it is given. The stderr summary groups the sites
by interval, with their mutation count and the widening that would cover
them, and warns if any site lies *inside* the mask, which means the profile
was not made from the tree that mask was applied to. Step 4 has the worked
case.

`--mask_pad` widens every masked feature on **both** sides, which
`--promoter_flank` does not: a promoter flank is upstream-only by definition,
so it reaches nothing downstream of a feature. Use it for alignment artifacts,
which spill past a feature's annotated edges in either direction. On a
129k-sample MTBC tree, nine of the twenty most recurrent sites left after
masking the resistance loci sat 1-56 bp *outside* a feature, seven of them
downstream; `--mask_pad 100` covered all nine. The two options are additive, so
a feature given both grows by their sum on its upstream side.

Names that don't resolve against your annotation are reported as NOT FOUND
rather than silently dropped. Matching is case-insensitive and covers every
name a feature carries — `gene`, `Name`, `locus_tag`, `old_locus_tag`,
`gene_synonym` and `Alias` — not just the one displayed in the output. That
matters more than it sounds: NCBI GFF3 often puts the locus tag in `Name` and
the name people actually look a gene up by in `gene_synonym`, so `mabA` and
`gidB` resolve to `Rv1482c` and `Rv3919c` even though neither is the displayed
name.

A NOT FOUND line is still worth reading rather than ignoring, because it has
three quite different causes: the annotation genuinely lacks that name (pair
the entry with its locus tag), the name is wrong (`clpC` where the genome has
`clpC1` and `clpC2`), or the gene is absent. Only the first is harmless, and
only if the locus is covered by another entry.

A position in no feature is reported as `intergenic`, with the nearest gene and
a signed distance to it (negative = upstream). Whole-sequence landmark records
are ignored for this: an NCBI GFF3 opens with a `region` line spanning the
whole chromosome, and GenBank uses `source` the same way, so every position
falls inside one. Counted as a hit, they would turn every intergenic site into
a nameless feature at distance 0 — and let `--mask_pattern` mask the entire
genome. They are still listed in the feature counts, then dropped with a note.

### `check_gene_names.py`

Used in steps [3](HOWTO.md#step-3-build-and-apply-a-region-mask) and
[8](HOWTO.md#step-8-run-the-search-with-bootstrap-and-check-where-it-saturates), before any
list is given to `--regions` or `--highlight`.

```
check_gene_names.py --names repair_genes.txt --annotation NC_012345.1.gff3.gz
```

Resolves each entry exactly as `annotate_sites.py` does, against the
annotation's gene, Name, gene_synonym, Alias, locus_tag and old_locus_tag
fields, with a gene and its CDS merged into one locus. The report has three
parts:

- **Shown under another name**: entries that matched through a synonym or
  locus tag, with the name output tables will use instead (`dinB1` appears as
  `dinX`).
- **Same locus as an earlier entry**: harmless in a `--regions` list, and
  often deliberate there, since a name plus its locus tag keeps a mask
  working on annotations that name the gene differently. In a `--highlight`
  list keep one: `split_drivers.py` counts each name as a gene of its own,
  which inflates its "by chance" counts.
- **To fix**: entries matching several loci (use the locus tag of the one
  meant) and entries not found, each with candidates. A candidate is a name
  within `--near` characters (default 2) that contains the entry or is
  contained in it, or a locus whose product mentions the entry as a word,
  proposed by locus tag when that is the annotation's only name for it. A
  candidate already in the list under another entry is marked. An entry with
  no candidate needs a search of the product lines by function.

Candidates are suggestions: check each one's product before putting it in the
list, because a locus tag resolves exactly and a wrong one misleads quietly.
The script exits 1 when anything is left to fix. `--seqid` restricts the
annotation to one sequence, `--max_candidates` (default 5) caps the
candidates shown per entry, and `--selftest` checks every case above on a
constructed annotation.

### `permutation_null.py`

Used in [step 7](HOWTO.md#step-7-calibrate---min_chi).

Destroys any association between mutation type and position in the tree, then
runs the same search and the same statistic, and reports the chi-squared values
it still finds. Set `--min_chi` comfortably above them. One run does all of
step 7: both nulls, the size-matched null and the hotspot check.

```
permutation_null.py --input_tree tree.pb.gz --fold --min_mutations 500 \
    --max_branch_length 250 --replicates 100 --seed 1 --nthreads 8 \
    --annotation genes.gff3.gz --reference ref.fna > null_report.txt
```

Two null models, both run by default (`--mode both`; `--mode sites` or
`--mode types` runs one):

- `types` shuffles mutation-type labels across all mutations, keeping
  topology and per-branch counts fixed. The strict independence null: it also
  removes the type clumping that recurrent sites cause, so it gives the lower,
  more permissive threshold. It scores the whole tree only.
- `sites` gives each site a single type (its most common one) and shuffles
  that assignment across sites. Which branches carry mutations at which site
  is untouched, so a site mutating 500 times inside one clade still dumps 500
  counts of one randomized type into that clade. This keeps the independence
  violation and gives the higher, more conservative threshold. Calibrate
  against it, especially with convergent selection at resistance loci. It
  follows the greedy cascade `--max_splits` (default 25) deep, to measure
  the null for smaller subtrees.

The unpermuted tree is always scored too, through the same code path.

The report has four sections:

| Section | What to do with it |
| --- | --- |
| Permutation null summary | Each null's replicate maxima (min, median, `--quantile`, max), how far `sites` sits above `types`, and the recommended `--min_chi`: `--margin` (default 2) times the `--quantile` (default 0.95) of the `sites` maxima. Bare, 1.5x and 2x lines are printed too. With fewer than about 5 / (1 − quantile) replicates it warns that the quantile is not resolved. The observed maximum's ratio to the null: under about 2x, the strongest real split is close to what chance produces |
| Null ceiling by subtree size | How far the chance maximum falls below the whole-tree figure as subtrees shrink. A `--min_chi` set from the root is that much too strict for deep splits. Below the table: the lowest `--min_chi` that clears the whole-tree null, with its margin over each bin's ceiling in the last column |
| Observed splits against a size-matched null | Each observed cascade split's margin over the `--null_quantile` (default 0.99) of the nulls from subtrees within `--null_window` (default 1.5x) of its size, as in `split_confidence.py`. It flags splits the recommended `--min_chi` (or `--flag_min_chi`) rejects despite a margin of at least `--margin`, and accepted splits with a margin under 1. `--min_nulls` (default 10) marks a ceiling drawn from too few null splits |
| Hotspot check | The largest whole-tree `sites` maxima and the clade behind each. It flags a clade that is the maximum in more than `--hotspot_fraction` (default half) of the replicates, or that holds the `--top_maxima` (default 2) largest maxima. A node and a child whose `mutations_below` agree within `--same_clade` (default 0.95) count as one clade. A flagged clade's `--hotspot_sites` (default 10) most recurrent sites are listed, with an `annotate_sites.py` command for those holding `--hotspot_min` (default 20) or more of its mutations; `--annotation` and `--reference` go into that command. Annotate them, and mask them if they look like an artifact |

Other options:

- `--min_mutations` / `--max_branch_length` should match what you intend for the
  real run; they change which nodes are eligible and so change the null.
- `--fold` scores the strand-folded 6-type spectrum (6 categories, 5 degrees of
  freedom) instead of 12. This must match how you intend to run
  `spectrumSplits.py`, which takes the same `--fold` flag: the statistic changes
  scale and the null changes shape, so a threshold calibrated one way is
  meaningless the other. See [step 5](HOWTO.md#step-5-decide-whether-to-fold).
- `--min_chi` (default 1) is where the null searches stop, not the threshold
  being calibrated. Keep it at 1, so the cascade runs its full depth.
- `--seed` (default 1): replicate r uses seed + r in both nulls, so they are
  paired, and a deeper `--max_splits` run's iteration-1 rows are identical to
  a shallower one's.
- `--nthreads` runs replicates concurrently where `fork` is available, and falls
  back to serial otherwise.
- `--prefix` (default `null`) names the per-split TSVs: `null_types.tsv` and
  `null_sites_cascade.tsv` (`null_sites.tsv` with `--max_splits 1`), one row
  per split scored per replicate. Step 10's `split_confidence.py --null`
  reads the cascade.
- `--selftest` checks the chi-squared implementation against a naive reference
  and against `scipy.stats.chi2_contingency`, the summary, and the size bins,
  ceilings, flags and both hotspot rules on a constructed cascade, including a
  clade under two IDs and the same clade read from an extracted tree. Run it
  once in the environment where `spectrumSplits.py` runs: where scipy is not
  importable the scipy comparison is SKIPPED with a note, so a pass there
  confirms only the naive reference. Confirmed on the MTBC server, where both
  agree on 200 random tables to 1.8e-12.

The observed cascade accepts every split at `--min_chi 1`, and so carves
subtrees differently from a production run at a higher threshold. Treat its
deep rows as approximate. For the production run's own splits,
`split_confidence.py` models the carving exactly. Its `ceiling` and window
are used here too, so the two scripts give the same margin for the same
size.

The statistic is reimplemented in pure Python rather than called through scipy,
so all candidate nodes can be scored in one vectorless pass without scipy's
per-call overhead. It matches `chi2_contingency` exactly for tables larger than
2 x 2 (no continuity correction); columns whose total is zero are skipped, where
scipy would raise.

### `split_drivers.py`

Used in [step 8](HOWTO.md#step-8-run-the-search-with-bootstrap-and-check-where-it-saturates).

```
split_drivers.py --spectra spectra.tsv --splits_log splits.log \
    --bootstrap bootstrap_*_splits_output.tsv \
    --input_tree tree.pruned.masked.pb.gz --max_branch_length 250 \
    --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --highlight repair_genes.txt --mask selection_loci.bed,alignment_loci.bed \
    > split_drivers.txt
```

Reads the search's results and prints, in order:

| Section | What to do with it |
| --- | --- |
| Split drivers | Every split's residual against its parent's, its top classes, one-site loss, bootstrap support and verdict |
| Splits to interpret | The `interpret` splits, and each one's defining-branch length; those longer than `--max_branch_length` are not annotated |
| Gene name check | `check_gene_names.py` on the `--highlight` list: entries shown under another name, duplicates, and entries matching several loci or none, with candidates. Fix those and rerun |
| Defining-branch annotation | The `interpret` splits' defining branches annotated as `annotate_sites.py` would, and the highlighted-gene hits on each |
| Branch hits | The hits against chance, below |

The last three need `--input_tree`, `--annotation`, `--reference` and
`--highlight`; without them the report says what to add.

`PREFIX.tsv` (`--prefix`, default `split_drivers`) has one row per split:

| Column | Meaning |
| --- | --- |
| `parent`, `chi_logged` | The subtree the split was found in, and its chi-squared, from the lines of `splits.log` above the bootstrap marker |
| `chi2_residual`, `drivers`, `top_share`, `classes_80` | The split's `spectra.tsv` row against its parent's, as a 2 x K table (K = 6 folded classes when `spectra.tsv` has them, else 12). `drivers` lists the top three classes with `+` or `-` for the direction in the split and each class's share of the chi-squared. `classes_80` is how many classes it takes to reach 80% |
| `site`, `site_mutations`, `site_loss` | With `--input_tree`: among the split's `--top_sites` (default 5) most recurrent sites, the one whose removal from both rows loses the most chi-squared, and that fraction |
| `support` | With `--bootstrap`: fraction of replicates splitting the node or one within `--neighbourhood` (default 10) node numbers |
| `verdict` | `interpret`, or the flags it failed: `one class` (`top_share` >= `--one_class_share`, default 0.8), `one site` (`site_loss` >= `--one_site_loss`, default 0.5), `weak support` (`support` < `--min_support`, default 0.5) |
| `branch_mutations` | With `--input_tree`: the mutations on the split's own branch, comma-separated for `annotate_sites.py --positions` |

The rows are residuals: each split minus the splits nested in it, which is
what `spectra.tsv` reports. The search scored each split before later splits
were carved out of either side, so `chi2_residual` is not `chi_logged`. It
answers a different question: which classes separate the rows you will read.
With `--input_tree`, each residual total is rebuilt from the tree and checked
against `spectra.tsv`; a mismatch means a different tree or
`--max_branch_length`, and is reported.

**Branch hits.** A hit is a site inside a highlighted gene. Only
nonsynonymous hits, stop gains included, count unless `--count all` is
given, and every result's "by chance" column, not its rank, says whether it
stands out:

| Section | What to do with it |
| --- | --- |
| Chance rate | Highlighted genes' length outside the `--mask` BEDs over the genome's, times the branches' nonsynonymous share of coding mutations. Check the names reported as not found |
| All branches | Counted hits against expectation summed over every branch, and the nonsynonymous share of all hits beside the share for every coding mutation on the branches |
| Per branch | Each branch with counted hits, by P. "by chance" is how many of the branches would reach that P by chance, from each one's own length, so short branches that cannot reach it do not inflate it |
| One gene, one branch | Genes with several counted hits on one branch, with the mutations. Hits a few bases apart may be one event; check them in the reads before believing them |
| One gene, several branches | Genes hit on two or more branches, against the branches expected to hit them from the gene's length and each branch's; "by chance" corrects across the gene list |

`PREFIX_hits.tsv` has the per-branch numbers. Each annotated branch is also
written as `node_N_branch.tsv`, with `annotate_sites.py`'s columns, in
`--branch_dir` (default the current directory; `''` skips them). `--flank`
(default 200) is the promoter window, and `--seqid` picks the sequence, as
in `annotate_sites.py`.

`--selftest` checks the chi-squared arithmetic, the residual walk and the
one-site test on a constructed tree; the Poisson tails, the chance counts,
the hit counting and the mask arithmetic; and the annotation of two
constructed branches from a small GFF3 and FASTA, including the name check.

### `branch_class_spectra.py`

Used in [step 9](HOWTO.md#step-9-check-the-splits-against-sequencing-quality). Splits
each accepted split's mutations into those on internal branches and those on
terminal branches, and compares the two, to tell a lineage signal from a
per-genome quality gradient. Needs `bte`.

```
branch_class_spectra.py --input_tree tree.pruned.masked.pb.gz \
    --spectra spectra.tsv --fold --max_branch_length 250 \
    --output branch_class.tsv > branch_class_summary.txt
```

- `--input_tree` must be the same file the search ran on, because node IDs are
  assigned when a tree is read (see
  [Node IDs are not stable across trees](METHODS.md#node-ids-are-not-stable-across-trees)).
  Split roots missing from the tree are reported with a warning.
- `--spectra` is the search's `--output_spectrum` TSV; split roots are read
  from its `Node_ID` column. `--splits` takes a plain list of node IDs instead.
- `--fold` and `--max_branch_length` must match the search. A branch over the
  cap contributes nothing, as in `spectrumSplits.py`.

Each node's mutations go to the nearest split root at or above it, which is
how `spectrumSplits.py` reports spectra. The report on stdout gives the
whole-tree spectrum by branch class, the correlation of each class with
divergence per tip computed separately on internal and terminal branches (the
test), and one line per split. `--output` writes two rows per split
(`internal`, `terminal`) with `Node_ID`, `Branch_Class`, `Tips`, `Mutations`
and the proportion of each class. `--selftest` checks the group assignment and
partitioning, and exits.

### `split_confidence.py`

Used in [step 10](HOWTO.md#step-10-rank-the-splits-before-interpreting-them), and again in
[Testing the method](METHODS.md#testing-the-method). Ranks the accepted splits by two
independent measures, margin over a size-matched null and bootstrap support. It reads only files earlier steps
produce, runs nothing on the tree, and needs nothing beyond the standard
library.

```
split_confidence.py --splits_log splits.log --spectra spectra.tsv \
    --null null_sites_cascade.tsv \
    --bootstrap_replicates bootstrap_by_replicate.tsv \
    --clade_samples split_clade_samples.tsv \
    --output split_confidence.tsv > split_confidence.txt
```

Inputs:

- `--splits_log`: the search's stderr. Bootstrap replicates print to the same
  stream; only the lines before the bootstrap marker are read.
- `--spectra`: the `--output_spectrum` TSV from the same run.
- `--null`: the cascade null from `permutation_null.py --max_splits`, as in
  [step 7](HOWTO.md#step-7-calibrate---min_chi).
- `--bootstrap_replicates`: `replicate<TAB>node` per line, made from the
  replicate outputs with

  ```
  for f in bootstrap_*_splits_output.tsv; do n=${f#bootstrap_}; n=${n%%_*};
    cut -f1 "$f" | grep ^node_ | sed "s/^/$n\t/"; done > bootstrap_by_replicate.tsv
  ```

  `--bootstrap_support` takes aggregate counts instead, but then support can
  only be counted at the exact node.
- `--clade_samples`: `node<TAB>sample` per line, from `matUtils extract -I
  <node> -u` over the accepted splits. It establishes which split is nested in
  which; without it the tested-against sizes are upper bounds and late-found
  splits' margins are understated.

Options:

- `--null_window` (default 1.5): how far, as a fold range, a null replicate's
  subtree size may be from the split's and still count as size-matched.
- `--null_quantile` (default 0.99): the quantile of those size-matched null
  values taken as the ceiling; 1 takes their maximum, the old behaviour,
  which rises with the replicate count. Step 7 says why 0.99.
- `--neighbourhood` (default 10): the node-number window within which a
  bootstrap replicate's split counts as the same split (0 disables).
- `--strong_margin` (default 5) and `--strong_support` (default 70%): tier A
  needs both, tier B one, tier C neither. An input that is missing counts as a
  failure.

Output columns: `node`, `search_root`, `clade_parent`, `chi2`,
`tested_in_mutations`, `null_ceiling`, `null_replicates`, `margin`,
`support_exact`, `support_neighbourhood`, `tier`. With `--output`, a short
summary (tier counts, the weakest margin, and any warnings) goes to stdout;
without it, the table goes to stdout and the summary to stderr. `--selftest`
runs the self-tests and exits.

### `split_composition.py`

Used in [step 11](HOWTO.md#step-11-check-what-each-split-contains-before-claiming-a-mechanism).

```
split_composition.py --clade_samples split_clade_samples.tsv \
    --clades clades.txt --confidence split_confidence.tsv \
    --output split_composition.tsv > split_composition.txt
```

Inputs: `--clade_samples` is the `node<TAB>sample` table from step 10,
`--clades` the `sample<TAB>clade` table from `matUtils summary
--sample-clades` (its header and `None` labels are skipped), and
`--confidence` the `split_confidence.tsv` that supplies tiers and order.
Nesting comes from sample membership, as in `split_confidence.py`, and each
split is reduced to its residual samples before counting.

Per split it prints the tier, all samples, residual samples, the number of
distinct clades, `other` (the share outside the top clade's own lineage, by
dotted names), the `--top` (default 3) clades with shares, and the flag
`mixed` (`other` >= `--mixed_share`, default 0.2). `--tiers` (default
`A,B`) picks which tiers to list; `--tiers all` also lists tier C and the
root, and works without `--confidence`. `--output` writes the same per
split as a TSV, with the top clades as `clade:count`. `--selftest` checks
the residual, the lineage naming and the flag on a constructed tree.

### `split_contexts.py`

Used in [step 11](HOWTO.md#step-11-check-what-each-split-contains-before-claiming-a-mechanism).
Tests every split for a trinucleotide-context shift inside its classes, on
all branches and on internal branches only, and summarises the splits that
pass on both. Needs `bte` and numpy.

```
split_contexts.py --input_tree tree.pruned.masked.pb.gz --spectra spectra.tsv \
    --fold --max_branch_length 250 --reference root_genome.fa \
    --mask selection_loci.bed,alignment_loci.bed \
    --composition split_composition.tsv > split_contexts.txt
```

Inputs:

- `--input_tree`, `--spectra`, `--fold` and `--max_branch_length` must be
  the search's. Residual totals are rebuilt from the tree and checked against
  `spectra.tsv`'s `Total_Mutations`. A `spectra.tsv` with folded columns and
  no `--fold`, or the reverse, is warned about.
- `--reference` must be the root's genome. More than 1% of mutations
  disagreeing with it is reported as a warning.
- `--mask` takes BED files. A context with a masked flank is not called.
- `--composition` (optional) is `split_composition.py --output`, for the
  summarised splits' clades.

Options:

- `--fold`: the 96 folded context classes in 6 classes. Without it, the
  192 forward-strand classes in the 12 types, for an unfolded search.
- `--against`. `clade` (the default) compares against the nearest
  enclosing split's residual. `search` compares against the subtree the
  search found the split in, read from `--splits_log`, as `split_drivers.py`
  does.
- `--permutations` (default 1,000) sets the branch shuffles per split.
- `--controls` (default 30, 0 to skip) sets the random non-split clades drawn
  per split. They are drawn within `--control_window` (default 2) fold of the
  split's size, each with `--control_permutations` (default 200) shuffles.
  A split whose window holds few clades gets a coarse `ctrl_p`, and the
  report says whether more `--controls` or a wider `--control_window` would
  help.
- `--max_q` (default 0.05) and `--control_alpha` (default 0.05): a split
  passes both nulls with q below the first and `ctrl_q`, the
  Benjamini–Hochberg-corrected `ctrl_p`, at or below the second.
- `--top` (default 3) sets how many contexts each run's detailed listing
  gives per split.
- `--summary_top` (default 0, all) trims the summary to the N splits passing
  both runs with the lowest all-branch `ctrl_p`, and `--summary_contexts` (default 2) sets
  the contexts per split in the summary.
- `--one_site` (default 0.2) is the top-site share at which a context is
  flagged as one site.
- `--min_class` (default 80) is the count below which a class is marked
  `few`.
- `--seed` (default 1). Each run starts its own random stream from it, so a
  run's numbers don't depend on the other.
- `--prefix` (default `split_contexts`) names the TSVs.

The report has one section per run. Each gives the calling statistics, the
share of controls that pass the shuffle, and one row per split sorted by
`p`. It then lists the top contexts, by contribution to G(within), for each
split that passes both nulls, each with the most frequent site's share of
it. A last section summarises the splits that pass both nulls on both runs,
names those that pass on only one, and gives the one-site check.

The TSVs, with the default prefix:

- `split_contexts.tsv` and `split_contexts_internal.tsv`, one row per split:
  `g_between`, `g_within`, `within_share`, `null_mean`, `excess`, `p` and
  `q`; `control_n`, `control_median`, `control_p` and `control_q`;
  `control_floor`, the smallest `ctrl_p` its control count allows, and
  `control_reachable`, whether that floor could reach `--control_alpha`;
  G(within) and p per class; and `top_contexts`.
- `split_contexts_96.tsv` and `split_contexts_internal_96.tsv` (`_192`
  without `--fold`), the 2 x 96 (or 2 x 192) table: counts, within-class
  shares, each context's G contribution, and the split's most frequent site
  in that context with its share of the count (`split_top_site`,
  `split_top_site_share`).
- `split_contexts_controls.tsv`, the random control clades of both runs
  (`branches`, `node`, `group`, `mutations`, `excess`, `p`), for
  `plot_split_contexts.py`.
- `split_contexts_table.tsv`, the summary: `split`, `comparator`, `clade`
  (top clades of the residual samples and flags, from `--composition`),
  excess, `ctrl_p`, `ctrl_q` and q on each run, `contexts` (the contexts carrying most
  of the all-branch G(within), with their shares of class, split vs
  comparator), `top_site_shares` and `max_site_share`.

`--selftest` checks the following, then exits:

- that G(between) and G(within) add up to the 2 x 96 and 2 x 192 G exactly;
- the shuffle's calibration, under the null and under a proportional class
  shift;
- its power for a one-context shift;
- the pipeline on a constructed tree, in both class sets;
- that a difference only in strand is invisible folded and a class
  difference unfolded;
- that control clades stop at nested splits, and that the internal-only run
  counts no terminal branch;
- `ctrl_q`, and the flag on a split whose control count can't reach
  `--control_alpha`;
- the summary's selection, ordering, contexts, site shares and
  `--summary_top`.

### `plot_split_contexts.py`

Used in [step 11](HOWTO.md#step-11-check-what-each-split-contains-before-claiming-a-mechanism).
Draws three figures from one `split_contexts.py` run. Needs only matplotlib.

```
plot_split_contexts.py --prefix split_contexts
plot_split_contexts.py --prefix split_contexts --format png --output figs/sc
```

It reads `PREFIX.tsv`, `PREFIX_internal.tsv`, the two context tables
(`_96` or `_192`), `PREFIX_table.tsv` and `PREFIX_controls.tsv`. Without the
controls file, from a `split_contexts.py` older than the figures, the
overview draws without the cloud and says so.

- `PREFIX_overview`: excess against split size, all branches and internal
  branches only, over the random control clades. Summary splits are filled
  and labelled. Splits that pass the shuffle but whose control count can't
  reach `--control_alpha` are crosses.
- `PREFIX_profiles`: one row per summary split. Bars are the change in each
  context's share of its class, split minus comparator, on all branches;
  dots, internal branches only. The summary's contexts are named. A
  triangle marks a context where one site carries `--one_site` (default
  0.2) or more of a count of at least `--site_min_count` (default 20).
- `PREFIX_heatmap`: every split tested, by all-branch excess, against every
  context, coloured by sign(share change) × √(G contribution). Cells below
  `--min_cell_g` (default 3.84, nominal p 0.05 at one degree of freedom) are
  drawn neutral.

Colours follow the context classes: one per folded class, with the class
names printed over the blocks. In the 192 classes, each folded pair of types
shares a colour and the purine-reference type is hatched. `--output` sets
the figures' prefix (default `--prefix`), and `--format` sets pdf, png or
svg. `--selftest` checks the share change, the signed deviance, the
colours, and that all three figures draw from constructed 96- and 192-class
files, and without controls.

### `gene_hit_contrasts.py`

Used in [step 12](HOWTO.md#step-12-test-a-candidate-mechanism-genotype-first). Finds every
branch with a coding SNV in a listed gene and tests whether the spectrum below
it differs from its siblings', pooled across all the independent hits. One
run tests one list, read three ways, and reports on them together. Needs
`bte`; it reuses `annotate_sites.py`'s annotation and codon code and
`permutation_null.py`'s spectrum code, so keep the three scripts together.

```
gene_hit_contrasts.py --input_tree tree.pruned.masked.pb.gz \
    --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --genes tb_repair_genes.txt --control control_genes.txt \
    --fold --max_branch_length 250 --permutations 10000 > hit_summary.txt
```

Inputs:

- `--annotation` and `--reference`: as for `annotate_sites.py`, and in the
  tree's coordinates. `--seqid` picks one sequence if they hold several.
- `--genes`: gene names or locus tags, one per line, `#` comments allowed;
  [`tb_repair_genes.txt`](#tb_repair_genestxt) works as it is.
- `--control` (optional): a control list from `control_genes.py`, tested
  the same way with every hit counted, as the baseline the list's excess is
  read against.
- `--fold` and `--max_branch_length`: match the search. A branch over the cap
  contributes no mutations but still counts for topology.

Options:

- `--effects` (default `nonsynonymous,synonymous`): comma-separated effect
  sets to test, from `nonsynonymous`, `lof` (nonsense and start-loss),
  `missense` and `synonymous`. The synonymous set is the built-in control.
- `--top_branch both` (the default) runs two readings. `include` counts each
  child's own branch as well as what is below it. It is the only way terminal
  hits contribute anything, and most branches are terminal. Its cost is
  dilution: mutations that predate the hit on the same branch are counted as
  if after it. That weakens a real effect but cannot manufacture one, since
  siblings are treated identically. `exclude` is the cleaner check, usable
  only on internal hits. `--top_branch include` or `exclude` runs one.
- `--min_focal` (default 1): the fewest counted mutations in the hit child
  for an event to be used. Siblings below it are not drawn as pseudo-hits.
- `--min_comparator` (default 10): the fewest mutations in the pooled
  siblings for an event to be used.
- `--permutations` (default 1000) and `--seed` (default 1). Use 10,000 when
  the list is long enough that its bar is below 0.001 (step 12). Each
  reading starts its draws from the seed, so it gives what a run of it
  alone would.
- `--null_draws_genes g1,g2` (or `ALL`, the pooled row): write every
  permutation's class shifts for those genes, in both readings, for
  `plot_gene_hits.py`.
- `--context_genes g1,g2`: test those genes, which must be in `--genes`,
  again in trinucleotide context classes, in both readings: the 96 folded
  ones with `--fold`, the 192 forward-strand ones without (needs numpy; see
  step 12). `--mask` gives the BEDs the tree was masked with, so masked
  flanks are unknown.
- `--gene_sets FILE`: group the genes into pathways, one `gene<TAB>set` per
  line ([`tb_repair_pathways.tsv`](#tb_repair_pathwaystsv) is the MTBC
  draft). Every gene and every set is then permuted together, in both
  readings, and the run adds the joint tests (step 12; needs numpy). Names
  must be ones `--genes` tests; others are reported and left out.
- `--joint`: the joint count and family-wise bar without gene sets.
- `--calibrate_context N`: instead of testing a tree, simulate N small trees
  in each of three scenarios and report how often the context tests reach
  p < 0.05: no shift, C>T raised 3× in proportion across its contexts, and
  one C>T context raised (step 12). It honours `--fold`, `--permutations`
  and `--seed`, and needs numpy but not `bte`.
- `--calibrate_joint N`: instead of testing a tree, simulate N small trees
  in each of three scenarios and report how often the joint tests (the gene
  count, and a two-gene set's pooled and any-direction tests) reach
  p < 0.05: no effect, both genes shifted the same way, and shifted opposite
  ways (step 12). It honours `--permutations` and `--seed`, and needs numpy
  but not `bte`.
- `--prefix` (default `hit`) names the outputs.

Outputs, with the default prefix (`_internal` for the internal-only
reading, `_control` for the control list):

- `hit_events.tsv`, `hit_events_internal.tsv`, `hit_events_control.tsv`: one
  row per hit event and effect set, with the gene, hit and parent node, codon
  effect and amino-acid changes, `nested_under`, how many siblings were used,
  the tips and mutations on each side, the event's chi-squared, and the class
  proportions on both sides (`focal_<class>`, `comparator_<class>`).
- `hit_summary.tsv`, `hit_summary_internal.tsv`, `hit_summary_control.tsv`:
  one row per gene and effect set, plus a pooled `ALL`, with event counts,
  `effective_events`, `diff_<class>` and `p_<class>`, and `omnibus` and
  `omnibus_p`. Step 12 says how to read them. The last columns are the
  heaviest-events check: `top_event_share` and `top2_event_share`, the
  share of the test's weight on its heaviest one and two events;
  `diff_without_top2_<class>`, the shift with those two dropped;
  `retained_without_top2`, how much of the shift that keeps (1 unchanged,
  0 or less gone, in the omnibus metric); and `heavy`, `yes` when under 50%
  is kept.
- `hit_null_draws.tsv` and `hit_null_draws_internal.tsv`, with
  `--null_draws_genes`.
- With `--context_genes`: `hit_summary_context.tsv` and
  `hit_events_context.tsv` (and `_internal`), whose summary adds
  `shape_p_<class>` and `shape_omnibus_p`, and `hit_context.tsv` (and
  `_internal`), one row per test and context. The events file still reports
  the 6 classes, or the 12 types.
- With `--gene_sets` or `--joint`: `hit_joint.tsv` (and `_internal`), one
  row per gene and effect set from the joint shuffle, with its omnibus p and
  `familywise_p`, the min-p-adjusted p; and with `--gene_sets`,
  `hit_sets.tsv` (and `_internal`), one row per set and effect set, with
  `genes`, `genes_with_events`, `events`, `informative`,
  `effective_events`, the pooled `diff_<class>` and `p_<class>`,
  `omnibus_p` (pooled, directional), `any_direction_p`, the heaviest-events
  check (`top2_event_share`, `retained_without_top2`, `heavy`) and
  `members`. The report marks a set HEAVY only where its pooled p is below
  0.05; the `heavy` column is filled for every set.

The report on stdout, in order:

1. Whether the synonymous control looks like chance: how many genes reach
   omnibus p < 0.05, and how many class tests reach p < 0.05, against what
   chance predicts, with a binomial tail p for the excess.
2. Whether the nonsynonymous hits exceed chance, and with `--control` the
   control genes: a one-sided Fisher test of the list's rate against the
   control's, and the control's pooled shift per class, the baseline a single
   gene's shift should be read against.
3. Each gene that clears the bar, 0.05 over the number of genes tested, with
   every class at p < 0.05 (or the best one) for all hits, internal hits and
   its own synonymous hits, with events and effective events.
4. The leads: genes reaching p < 0.05 without clearing the bar. Chance alone
   gives about one in twenty.
5. With `--gene_sets` or `--joint`, GENE SETS AND JOINT TESTS, per reading
   and effect set, synonymous first: the count of genes at p < 0.05 against
   the count expected with every gene permuted together, the family-wise
   bar against 0.05 over the number of genes and the genes clearing it, and
   per set the pooled p, its largest class shift and the any-direction p.
6. With `--context_genes`, the shape p values per gene and reading.

It also says when the smallest p is the permutation floor (more permutations
could only lower it), when the bar is below that floor (no gene could clear
it), and when a test rests on few effective events. Each gene that clears the
bar gets its shift without its two heaviest branches, and any gene or set
whose shift doesn't survive that is marked HEAVY. The weights grow with
subtree size, so a hit on a deep, lineage-defining branch compared with a
single sibling can carry a whole test; it is then one lineage against
another. The number of
permutations is read from the p-values, which are (k + 1) / (permutations +
1). Chance counts use a binomial tail at 5%; neither they nor the bar correct
for the six classes being correlated, so read them as guides.

`--selftest` checks codon effects on both strands (including ancestral codon
state and two changes in one codon), checks the cut-subtree bookkeeping
against brute force, and simulates star-of-polytomy trees for null
calibration and power. In a larger run of the same simulation, 1,000 null
trees of 100 hit subtrees each gave omnibus p < 0.05 at 4.9%. It checks the
report on a constructed summary, and that each reading writes its files and
reproduces whatever ran before it. With numpy it checks the joint tests: a
one-gene set is permuted with exactly its gene's draws; under no effect the
count and both set tests stay near 5% (`--calibrate_joint` measures this
on more trees); two
genes shifting the same way are found by the pooled test; and two shifting
opposite ways cancel there but are found by the any-direction test. With numpy it also checks the context test: the
summed contexts equal the 6-class statistic for the same events; the shape
test stays near 5% with no effect and with a whole-class proportional shift;
and it finds a shift in one context.

The simulated trees draw each hit child with the size-biased weight the test
assumes, so the hit child is exchangeable with its siblings under no effect.
(Before November 2026 they drew it uniformly; the calibration figures in
step 12 were remeasured after the change.)

### `control_genes.py`

Used in [step 12](HOWTO.md#step-12-test-a-candidate-mechanism-genotype-first), which has
the command. Draws a length-matched random control list for
`gene_hit_contrasts.py`.

For each gene of `--match`, `--per_gene` (default 3) protein-coding genes
whose CDS length is within `--tolerance` (default 0.2) of it are drawn at
random, or the closest ones if too few qualify, never the same gene twice. A
gene is not eligible if it is named in `--match` or an `--exclude` list, if
more than `--max_masked` (default 0.1) of its CDS lies in the `--mask` BED
intervals, or if its product matches `--exclude_product`. The output is a
gene list of locus tags, each with a comment naming the gene it was matched
to. `--seed` (default 1) makes the draw reproducible. `--selftest` checks the
exclusions, the length window, the fallback and the seed.

### `plot_gene_hits.py`

Used in [step 12](HOWTO.md#step-12-test-a-candidate-mechanism-genotype-first), which has
the command. Needs only matplotlib, so it can run outside the `bte`
environment.

- `--prefix`: the `gene_hit_contrasts.py` run's prefix. It reads
  `PREFIX_summary.tsv` and `PREFIX_null_draws.tsv`, and the internal-only
  reading's `_internal` files, for a third row in panel B, when they exist.
  `--summary`, `--null`, `--internal_summary` and `--internal_null` name the
  files one by one instead.
- `--gene`, `--class`: the gene, and the classes to show as the summary
  columns name them, comma-separated (`'CG|GC,CT|GA'`), one column each.
  `gene_hit_contrasts.py`'s report lists the classes worth showing.
- `--label`: other genes to name in panel A. `--title` replaces the default
  title. `--output`: the format follows the extension (`.png`, `.pdf`,
  `.svg`).

### `resistance_spectra.py`

Used in [step 13](HOWTO.md#step-13-ask-whether-the-splits-spectra-show-in-their-resistance-mutations), which has
the command. Finds every resistance substitution on the unmasked tree, from
a catalogue of graded variants, and asks of each split and drug whether its
spectrum changes how many it acquires (rate) or which (composition). Needs
`bte`, numpy and `openpyxl`; it reuses `annotate_sites.py`'s annotation
code, `trinucleotide.py`'s contexts and `split_contexts.py`'s helpers.

Inputs:

- `--input_tree`: the tree from before the mask. On the masked tree these
  sites have no mutations.
- `--clade_samples`: every sample of every split, as assembled in step 10.
  Each split is placed at its samples' most recent common ancestor. A split
  whose ancestor's clade has extra samples, or some of whose samples are
  missing, is reported.
- `--catalogue`: the WHO catalogue's Web Annex B (`.xlsx`), or a catalogue
  in its layout:
  - `Graded_Mutation_List`, with `drug`, `Gene`, `variant`, `effect`,
    `SUPPLEMENTARY GRADING CONSIDERATIONS` and `FINAL CONFIDENCE GRADING`;
  - `Gene_Coordinates`, with `Gene`, `Rv ID`, `Left` and `Right`.
- `--groups` (default `1,2`): the confidence groups that count as
  resistance.
- `--annotation`: each catalogue gene's strand and type, found by the
  catalogue's locus tag. Where its coordinates differ from the catalogue's,
  the catalogue's are used and the difference is reported.
- `--reference`: the root genome.
- `--mask`: the masks the search's tree was built with. Spectra and
  opportunity are counted outside them; resistance sites are read
  regardless.
- `--fold`, `--max_branch_length`: match the search.

Options:

- `--prior` (default 100): the pseudo-mutations shrinking each group's
  spectrum towards the tree's.
- `--min_events` (default 5): the events a split needs for its composition
  to be tested.
- `--permutations` (default 10000): label shuffles for G, and spectrum-ratio
  permutations for llr, per drug and reading.
- `--max_q` (default 0.05): the q below which the report lists a split under
  a drug.
- `--seed` and `--prefix` (default `resistance_spectra`).

Outputs:

- **The report**, on stdout:
  - the catalogue table (graded rows, those one substitution makes, those
    it can't reach, the rules, and what is invisible);
  - events per drug, with back-mutations;
  - then, per reading, the drug table, the `any drug` table per split, and
    the per-drug hits.
- **`PREFIX_events.tsv`:** one row per event and drug. Columns: `node`,
  `terminal`, `split`, `pos`, `parent_base`, `alt`, `drug`, `gene`,
  `variant` (named against the parent), `effect`, `group`, and `route` (the
  catalogue's supplementary grading, `listed`, or `Rule: LoF (not listed)`
  / `Rule: RRDR (not listed)` for a rule reaching an unlisted variant).
- **`PREFIX.tsv`:** one row per reading, drug and split. Columns:
  - `events` and `mutations`;
  - `e_global`, `e_spectrum`, `spectrum_factor` and `obs_over_spectrum`;
  - `llr`, `llr_expected`, `llr_p` and `llr_q`;
  - `g`, `g_p` and `g_q`;
  - `top_variant`, the variant most over-represented against the rest, as
    observed / expected.
- **`PREFIX_drugs.tsv`:** one row per reading and drug. Columns: `events`;
  `variants_seen`; `substitutions` (the site-base pairs that make a
  resistance variant from the reference); `splits_tested`; `llr_total` and
  its p; `g_all` and its p.
- **`PREFIX_variants.tsv`:** the shares behind each llr, for
  `plot_resistance_spectra.py`. One row per reading, drug, tested split and
  variant: each variant the split has, or that makes up 0.5% or more of q
  or q'. Columns: `observed`, `split_events`, `observed_share`,
  `rest_share` (q) and `spectrum_share` (q').

`--selftest` checks:

- naming on both strands, including upstream, rRNA, start-lost, stop-lost
  and a background codon change;
- the catalogue's precedence: an explicit grade over a rule, the LoF and
  RRDR rules, and a group outside `--groups`;
- events on a small tree, with a back-mutation that would otherwise be named
  as resistance;
- split placement;
- the G statistic against the textbook formula, and llr_expected against
  simulated draws;
- that a spectrum-driven composition is found, and that the variants
  table's shares reproduce its llr;
- that a composition shifted for other reasons is not called a spectrum
  effect by the ratio permutation (2 to 4 in 100), where label shuffling
  would (27 to 34 in 100).

### `plot_resistance_spectra.py`

Used in [step 13](HOWTO.md#step-13-ask-whether-the-splits-spectra-show-in-their-resistance-mutations).
Draws three figures from one `resistance_spectra.py` run. Needs only
matplotlib.

```
plot_resistance_spectra.py --prefix resistance_spectra \
    --composition split_composition.tsv
plot_resistance_spectra.py --prefix resistance_spectra --format png \
    --example_drug Isoniazid --example_split node_36254 --example_reading all
```

It reads `PREFIX.tsv` and `PREFIX_variants.tsv`. `--composition`
(`split_composition.py`'s output) names each split's main lineage and
colours it:

- L1 to L4 each get a colour;
- *M. africanum* (L5, L6 and La) shares one;
- anything else is gray.

`--drug` (default `any drug`, every drug pooled) picks the drug of the
first two figures.

Captions:

- **`PREFIX_rate`.** *Resistance supply: what each split's spectrum
  predicts, and what it got.* Each row is a split, on one log axis of
  resistance events relative to the tree-wide expectation, `E_global` (the
  dashed line at 1).
  - The hollow circle is `spectrum_factor`: how many resistance mutations
    the split's own spectrum would make available, relative to the tree's.
  - The filled circle, coloured by lineage, is what the split acquired.
  - Rows are sorted by the all-branch observed rate. The right panel, for
    internal branches only, keeps that order.
  - A split with no events is a left-pointing triangle at the axis's left
    edge.

  How far the hollow circles stray from 1 is how much mutation supply
  could matter. How far the filled ones stray is how much treatment and
  transmission do.

- **`PREFIX_composition`.** *Resistance variants: does a split's spectrum
  predict which it acquires?* One point per split tested
  (`--min_events`), on all branches and on internal branches only.
  - x is `llr_expected`, the lean towards the spectrum's favoured variants
    that the spectrum alone would produce. y is the observed `llr`.
  - On the dashed diagonal, a split leans its spectrum's way by as much as
    its spectrum predicts.
  - Above the diagonal, it leans further than the spectrum explains. Below
    zero, it leans against its spectrum.
  - Filled, labelled points pass at llr q below `--max_q`; the p value
    comes from permuting the spectrum ratios across variants.
  - Both axes are symmetric-log: linear within 1 of zero, logarithmic
    beyond.

- **`PREFIX_example`.** *What an llr is made of.* For one split and drug
  (`--example_split`, `--example_drug`, `--example_reading`; default
  `node_14336`, Streptomycin, internal), its `--top_variants` (default 8)
  commonest variants, each with three bars:
  - gray: the variant's share of the rest of the tree's events (q);
  - orange: the share the split's spectrum predicts (q');
  - blue: the split's own share, with its count.

  The title gives the split's llr, llr_expected and q. llr sums log(q'/q)
  over the split's events, so it is positive where blue follows orange
  away from gray. Blue going well past orange is a lean the spectrum
  alone doesn't explain.

`--output` sets the figures' prefix (default `--prefix`), and `--format`
sets pdf, png or svg. `--selftest` checks the lineage names and colour
groups, and that all three figures draw from constructed files, with and
without lineages, and without a variants table.

### `match_splits.py`

Used in [Testing the method](METHODS.md#testing-the-method).

```
match_splits.py --old split_clade_samples.tsv \
    --new split_clade_samples_nobcg.tsv --output split_matches.tsv > split_matches.txt
```

`--old` and `--new` are `node<TAB>sample` tables for two searches, each
made by running `matUtils extract -I` on every node in that run's
`spectra.tsv` (the loop in step 10). Both are reduced to the samples present
in both trees, taken from their `node_1` rows, before anything is compared,
and the root itself is not matched. Each old split is paired with the new
split of highest Jaccard index, and each new split with its best old one;
a match is mutual when the two agree. The summary counts exact matches,
near ones (`--near`, default 0.95) and old splits with no counterpart
(`--lost`, default 0.5). `--output` writes the per-split rows as a TSV.
`--selftest` checks the shared-sample reduction, partial matches and the
report on a constructed pair of trees.

### `trinucleotide.py`

Used by `profile_tree.py` ([step 5](HOWTO.md#step-5-decide-whether-to-fold)),
`split_contexts.py` ([step 11](HOWTO.md#step-11-check-what-each-split-contains-before-claiming-a-mechanism))
and `gene_hit_contrasts.py --context` ([step 12](HOWTO.md#step-12-test-a-candidate-mechanism-genotype-first));
not run on its own except as `trinucleotide.py --selftest`.

It calls every SNV's context in one depth-first pass over the tree. The
positions changed on the current root-to-node path are kept in a dictionary,
and each node's changes are applied on the way down and undone on the way back
up. So each flank is the base in the mutation's parent genome, not the
reference's, at a cost of one step per mutation. The classes are the 96
pyrimidine-centred ones, in the usual signature order: C>A, C>G, C>T, T>A,
T>C, T>G, matching the folded classes' order, then 5' and 3' base. A purine
mutation is reverse-complemented along with its trinucleotide.

A mutation gets no context when a flank is N, masked or off the end of the
sequence (ends are not wrapped, even for a circular genome), or when the
adjacent position changed on the same branch. Double substitutions are a
signature of their own and are counted separately rather than called twice.
Folding to 96 has the same blind spot as `--fold`: a leading- or
transcribed-strand process cancels between the replication arms.

For an unfolded tree there are also 192 forward-strand classes: the 12 types
in the order of the 12-type spectrum, each with its 16 flank pairs as read on
the reference's forward strand, with no reverse complement. They need no
notion of which strand leads or is transcribed. `scheme(fold)` bundles
either set's size, labels and opportunity counts, so `profile_tree.py`,
`split_contexts.py` and `gene_hit_contrasts.py` choose in one place. Oriented
classes, by leading or coding strand, would need the `--orient` recoding that
is parked for `spectrumSplits.py`.

`--selftest` checks the class order and folding against `folded_index()`,
flanks that change on an ancestral branch, reversions, sibling isolation,
doublets, masks and a wrong reference. It also checks 30 random trees against
a brute-force walk up the tree, and the opportunity counts. For the 192
classes, it checks that each strand-specific type gets its own class, that
the 192 calls fold onto the 96 on random trees, and that the 64
forward-strand opportunity counts fold onto the 32.

## Gene lists

Four lists of *M. tuberculosis* loci ship with these scripts, for
`annotate_sites.py --regions` and `--highlight` and `gene_hit_contrasts.py
--genes`. Two are masked and two are deliberately not. A fifth file groups
the unmasked genes into pathways for `gene_hit_contrasts.py --gene_sets`. Run [`check_gene_names.py`](#check_gene_namespy) on any of
them, or on a list of your own, against the annotation you will use.

### `tb_resistance_genes.txt`, and a caution

`tb_resistance_genes.txt` holds gene **names and locus tags only — no
coordinates**. Every interval is resolved from the annotation *you* supply, so
the coordinates come from your reference rather than from this repository. That
is deliberate: a name you can check beats a coordinate you cannot.

The list is a starting point and is not authoritative. The authoritative source
is the WHO *Catalogue of mutations in Mycobacterium tuberculosis complex and
their association with drug resistance*, 2nd edition (2023), which publishes
genomic coordinates and is what you should cite. Cross-check against
TB-Profiler's database too. This script does not predict drug resistance; for
resistance calls use TB-Profiler or Mykrobe.

### `tb_repair_genes.txt`

`tb_repair_genes.txt` lists DNA repair, replication-fidelity and dNTP-pool
genes, and exists for the opposite purpose. **Do not mask these.** A
repair-gene mutation is the *explanation* you want when a clade turns out to
have a distinct spectrum: a lineage that lost or altered a repair function is
what a real mutational-spectrum split looks like. Masking them destroys the
signal.

Use it with `--highlight`, which flags matching sites in a `highlight` column
and counts them, changing nothing else:

```
annotate_sites.py --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --positions <split_drivers.py's branch_mutations for the clade> \
    --highlight tb_repair_genes.txt
```

This list carries gene names, and locus tags only where they were checked. A
name missing from your annotation is reported as NOT FOUND, which is harmless,
whereas a misremembered locus tag would resolve to the wrong gene and mislead
you silently. So most Rv numbers are in comments, as hints to check. The five
locus-tag entries (Rv1021 mazG, Rv1321 nucS, Rv0944 fpg2, Rv2464c nei1, Rv2228c
RNase HI) are genes that NCBI's NC_000962.3 annotation names only by locus tag,
each checked against that annotation's product line; without them, nucS, the
mycobacterial mismatch-repair surrogate, went unmatched. With them, every entry
resolves against that annotation. Verify any you rely on against Mycobrowser or
the NCBI H37Rv annotation.

### `tb_replication_genes.txt`

`tb_replication_genes.txt` lists replisome and dNTP-supply genes that
`tb_repair_genes.txt` does not: `dnaA`, `dnaB`, `dnaG`, `dnaZX`, `ssb`,
`priA`, `topA`, the ribonucleotide reductases and their partners, `ndkA`,
`dcd` and `tmk`. Like the repair list, it is never masked. It is kept
separate so that each list is its own family in step 12: adding its 17 genes
to the repair list would raise every repair gene's bar. `gyrA`, `gyrB`,
`thyA` and `thyX` are in the resistance list, so they are masked and carry no
hits. Every name resolves against NC_000962.3 (`check_gene_names.py`).

### `tb_repair_pathways.tsv`

`tb_repair_pathways.tsv` assigns the 73 genes of `tb_repair_genes.txt` and
`tb_replication_genes.txt`, plus four that NCBI's NC_000962.3 names only by
locus tag, to 13 pathways. The four are `Rv3395c` (`imuA'`) and `Rv3394c`
(`imuB`), the rest of `dnaE2`'s mutasome, and `Rv3202c` (`adnA`) and
`Rv3201c` (`adnB`), the AdnAB end-resection helicase-nuclease. The pathways
are:

- nucleotide-pool sanitation, ribonucleotide excision, alkylation, base
  excision, the mismatch surrogate NucS and nucleotide excision;
- the `dnaE2`/`imuA'`/`imuB` mutasome and the Y-family translesion
  polymerases;
- the replisome;
- RecA-dependent recombination and SOS, with AdnAB;
- single-strand annealing, with RecBCD. In mycobacteria RecBCD works in
  this RecA-independent pathway, not in recombination (Gupta et al. 2011,
  *Mol Microbiol* 79:316, in *M. smegmatis*).
- end joining and dNTP supply.

It is a draft for review: some genes could defensibly sit in another set,
and the groupings should be settled before looking at per-gene results.
None of the four locus-tag genes is in a list; add them to `--genes` as the
file's header shows. The header also rebuilds the first draft, with RecBCD
in recombination and no AdnAB, which METHODS.md step 12 reports beside the
final grouping.

### `tb_alignment_loci.txt`

`tb_alignment_loci.txt` is also masked, but for the opposite reason to
`tb_resistance_genes.txt`. Resistance loci recur because the same mutation is
selected repeatedly — the variants are real. These loci recur because short
reads cannot be placed unambiguously: paralogous families (PE/PPE, `esx`, the
`Rv1318c`-`Rv1320c` adenylate cyclases, `lldD2`), repeats (the 13E12 family,
and transposases by product) and prophages (phiRv1, phiRv2). There the
variants are largely not real at all. NC_000962.3 annotates its insertion
sequences, IS6110 among them, as mobile elements with no gene, so nothing in
the list reaches them; the list explains why that left no hole on this tree.

Keeping the two lists apart keeps each auditable, and they want different
widening, so [step 3](HOWTO.md#step-3-build-and-apply-a-region-mask) resolves them in two
separate passes and merges the results.

The strongest signal that a site belongs in this list rather than the
resistance list is in `annotate_sites.py`'s `aa_change` column: if every
possible substitution at a position is synonymous, or a third codon position
recurs hundreds of times, selection cannot be the cause.

Two cautions. The list carries a handful of genes whose recurrence has **no
established mechanism** — `ponA1`, `hsdM` and `Rv2024c` are nonsynonymous at
every alternate, so they may well be convergent selection rather than
miscalling. They are listed because they break the independence assumption
either way, but if you are looking for the biology, look at them before you
mask them. And prefer `--mask_pattern` for the large families: it matches names,
so it picks up members the file does not list.
