# Parameter tuning helpers

Read-only scripts for choosing `spectrumSplits.py` parameters on a tree other
than the one the defaults were calibrated on, and for checking the splits it
then finds. None of them modifies the tree. They need only the Python standard
library and `bte`, and six need less: `prune_clade_check.py` and
`split_drivers.py` use `bte` only with `--input_tree`, and `null_by_size.py`,
`branch_hits.py`, `split_confidence.py` and `split_composition.py` read text
files only.

## Contents

- [Background](#background): why the defaults do not transfer between trees, and two properties of a tree to know before reading any result from it.
  - [Why the spectrumSplits.py defaults need calibration](#why-the-spectrumsplitspy-defaults-need-calibration)
  - [Know how your tree was rooted](#know-how-your-tree-was-rooted)
  - [Node IDs are not stable across trees](#node-ids-are-not-stable-across-trees)
- [Suggested order](#suggested-order): the walkthrough, start to finish, with worked numbers from an *M. tuberculosis* tree. This is the part to work from; the step numbers in the script table refer to it.
  - [1. Profile the raw tree](#1-profile-the-raw-tree)
  - [2. Prune samples with implausible divergence](#2-prune-samples-with-implausible-divergence)
  - [3. Build and apply a region mask](#3-build-and-apply-a-region-mask)
  - [4. Re-profile, and check what the mask bought](#4-re-profile-and-check-what-the-mask-bought)
  - [5. Decide whether to fold](#5-decide-whether-to-fold)
  - [6. Pick `--min_mutations` and `--max_branch_length`](#6-pick---min_mutations-and---max_branch_length)
  - [7. Calibrate `--min_chi`](#7-calibrate---min_chi)
  - [8. Run the search with bootstrap, and check where it saturates](#8-run-the-search-with-bootstrap-and-check-where-it-saturates)
  - [9. Check the splits against sequencing quality](#9-check-the-splits-against-sequencing-quality)
  - [10. Rank the splits before interpreting them](#10-rank-the-splits-before-interpreting-them)
  - [11. Check what each split *contains* before claiming a mechanism](#11-check-what-each-split-contains-before-claiming-a-mechanism)
  - [12. Replicate the search on a changed tree, and match by membership](#12-replicate-the-search-on-a-changed-tree-and-match-by-membership)
- [Script reference](#script-reference): usage, inputs, outputs and options for each script.
  - [`profile_tree.py`](#profile_treepy)
  - [`prune_clade_check.py`](#prune_clade_checkpy)
  - [`annotate_sites.py`](#annotate_sitespy)
  - [`permutation_null.py`](#permutation_nullpy)
  - [`null_by_size.py`](#null_by_sizepy)
  - [`split_drivers.py`](#split_driverspy)
  - [`branch_hits.py`](#branch_hitspy)
  - [`branch_class_spectra.py`](#branch_class_spectrapy)
  - [`split_confidence.py`](#split_confidencepy)
  - [`split_composition.py`](#split_compositionpy)
- [Gene lists](#gene-lists): the lists `annotate_sites.py` reads, and what each is for.
  - [`tb_resistance_genes.txt`, and a caution](#tb_resistance_genestxt-and-a-caution)
  - [`tb_repair_genes.txt`](#tb_repair_genestxt)
  - [`tb_alignment_loci.txt`](#tb_alignment_locitxt)

| Script | Step | What it answers |
| --- | --- | --- |
| [`profile_tree.py`](#profile_treepy) | 1, 4 | What is in this tree, and what do its numbers imply for `--min_mutations`, `--max_branch_length` and whether to fold? |
| [`prune_clade_check.py`](#prune_clade_checkpy) | 2 | What did `qc/prune_mutation_sample_ratio.py` remove, by clade, and was it a divergent genome or a clade on a long stem? |
| [`annotate_sites.py`](#annotate_sitespy) | 3, 4 | What are the recurrent sites, and which regions should be masked before anything else is trusted? |
| [`permutation_null.py`](#permutation_nullpy) | 7 | How large a chi-squared does chance alone reach on this tree, and so where should `--min_chi` sit? |
| [`null_by_size.py`](#null_by_sizepy) | 7 | How does that chance ceiling fall with subtree size, which splits does a fixed `--min_chi` wrongly reject, and is one hotspot setting the threshold? |
| [`split_drivers.py`](#split_driverspy) | 8 | Which mutation classes drive each split, does one site or weak bootstrap support explain it, and what mutations define the ones worth interpreting? |
| [`branch_hits.py`](#branch_hitspy) | 8 | Are the highlighted-gene hits on the split branches more than chance, per branch, per gene on one branch, and per gene across branches? |
| [`branch_class_spectra.py`](#branch_class_spectrapy) | 9 | Is a split a real mutational difference, or a sequencing-quality gradient on terminal branches? |
| [`split_confidence.py`](#split_confidencepy) | 10, 12 | Which accepted splits are worth interpreting? Scores each against a size-matched null and the bootstrap replicates. |
| [`split_composition.py`](#split_compositionpy) | 11 | What is each split made of: one lineage, several, or partly a laboratory strain? |

The steps depend on one another, so the order matters more than any single
setting: masking changes the spectrum, the spectrum sets the floor under
`--min_mutations`, and both reshape the null that `--min_chi` is calibrated
against.

## Background

### Why the spectrumSplits.py defaults need calibration

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

### Know how your tree was rooted

Ask before you interpret anything near the root. The MTBC tree used here was
built against NC_000962.3 as both reference and root, then rerooted to a node
picked as basal-looking, with a new reference made by applying that node's
substitutions to NC_000962.3. Coordinates stay H37Rv throughout, so nothing
about masking or annotation changes.

What it does change is that branch lengths around the root are a convention of
the tree build rather than a measurement. Three consequences:

- The root residual (`node_1` in the spectra output) is not a baseline for
  anything. On this tree there is a second reason as well -- once the search
  accepted its splits, `node_1`'s residual was down to 260 tips.
- Any criterion that fires on a clade's stem can be decided by the rooting for
  clades near the root. Step 2 of the walkthrough has the worked case.
- Claims about the deepest branches need the rooting's provenance attached.

Splits well inside named lineages are unaffected, which on this tree is all of
them that survived ranking.

### Node IDs are not stable across trees

`matUtils` assigns node IDs when it reads a tree, so `node_5289` in one
`.pb.gz` has nothing to do with `node_5289` in another. Masking makes this
bite: a node whose every mutation is masked away is collapsed, which changes
the topology, so IDs shift. Pruning does the same.

A split you want to compare across trees therefore has to be identified by
**its samples**, and the sample list has to be captured *before* the old tree
is overwritten:

```
matUtils extract -i old.pb.gz -I node_5289 -u node_5289_samples.txt
```

To find the same clade in the new tree, select that sample list and let
`-U`/`--from-mrca` expand the selection to everything descended from their
MRCA. `-U` takes no argument of its own; it modifies the `-s` selection:

```
matUtils extract -i new.pb.gz -s node_5289_samples.txt -U -u found_samples.txt
```

Then compare `found_samples.txt` against the input list. Identical means the
clade survived masking intact; larger means it was fragmented and the MRCA has
reached up past it, and the extra samples tell you what it reached over. If you
need the new node's ID itself, add `-A`/`--all-paths`, or extract the subtree
with `-o` so the MRCA becomes the new root.

Writing the masked tree to a new filename instead of reusing the old one keeps
both trees available, which is the easier fix.

Collapsing has a second consequence worth knowing before you re-read the
per-branch distribution: a collapsed node's mutations do not disappear, they
consolidate onto fewer, longer branches. The branch-length tail can shift
upwards even though the total mutation count went down, so any
`--max_branch_length` cap chosen from the old profile has to be re-derived
rather than carried over.

## Suggested order

A new tree, start to finish. File names are generic: `tree.pb.gz` is the
starting tree, `NC_012345.1.gff3.gz` and `NC_012345.1.fna` the reference
annotation and sequence, `resistance_genes.txt` a list of loci to mask.

The worked examples are measurements from a real *M. tuberculosis* tree
(129,314 samples, 4.4 Mb reference), not invented figures, so the reasoning
can be checked against the reports it came from. They are not
recommendations: another species, or another sampling scheme, will give
different numbers and sometimes different decisions. Step 8's figures are
marked hypothetical because that run has not been done.

### 1. Profile the raw tree

```
profile_tree.py --input_tree tree.pb.gz > profile_raw.txt
```

Every decision below is read off this report, so keep it.

### 2. Prune samples with implausible divergence

```
qc/prune_mutation_sample_ratio.py --input_tree tree.pb.gz \
    --threshold 250 --output_tree tree.pruned.pb.gz > pruned_report.txt
```

Choose `--threshold` from the MUTATIONS:TIPS section rather than accepting the
auto-threshold, and keep stdout: the script prints one line per pruned node
with the tips it took, which is the only record of what was removed. Note that
it prunes the whole subtree at the first node that crosses the threshold and
does not recurse into it, so one crossing can remove several samples.

*Worked example.* The report shows mean 18.6, median 7, p90 54, p99 114,
p99.9 192, max 909, and an auto-threshold (mean + 2 SD) of 71.7 with 10,930
nodes at or above it. Taking the auto-threshold would prune 10,930 nodes — 5%
of the tree — on a species where a high ratio is often real divergence. A
threshold of 250, well above p99.9, removed only 60 samples and still cut the
worst terminal branch from 909 mutations to 234.

Read the distribution before choosing. A median of 7 with a p90 of 54 is not a
tight distribution with outliers, it is a broad one, and mean + 2 SD on a broad
distribution cuts into the bulk.

#### Check what it removed, by clade

Run this first, so you have your own tree's numbers beside the worked example.
The removed samples' clades come from the unpruned tree, since the pruned one
no longer contains them, and the prune report already lists every removed
sample:

```
matUtils summary -i tree.pb.gz -s samples_raw.txt
matUtils summary -i tree.pb.gz --sample-clades clades_raw.txt
prune_clade_check.py --samples samples_raw.txt --clades clades_raw.txt \
    --prune_report pruned_report.txt --threshold 250 \
    --exclude ReconstructedAncestor --input_tree tree.pb.gz > prune_check.md
```

`prune_clade_check.py` writes, as markdown, every table in this step and the
next one. Its sections match the tables below:

| `prune_check.md` section | README table |
|---|---|
| 1. Removal rate by clade | clade / in tree / removed / rate / vs overall |
| 2. Terminal branches >= 250 vs. removals | median terminal / tips >= 250 / pruned |
| 3. Prunings at internal nodes | tip terminals / stem / internal below |
| 4. Surviving parents of pruned subtrees | the parents in "Where the pruned samples sat" |
| 5. Unlabelled samples | the unlabelled-parent counts and run table |

It lists every clade, while the README tables leave out the small ones.
Clades are collapsed to their first dot-separated part (`--depth 0` keeps full
names). `--threshold` only sets the `tips >=` column, so give it the same
value you passed to the prune script. `--input_tree` needs `bte` and adds the
exact stems (section 3) and the parent table (section 4); without it, the stem
and the internal branches below it are reported together as one remainder.

An absolute mutations:tips threshold is a divergence filter, not a quality
filter, and on this tree the difference is measurable. Of the 60 samples
pruned at 250 (smaller clades with no removals omitted):

| clade | in tree | removed | rate | vs overall |
|---|---|---|---|---|
| lineage8 | 3 | 3 | 100% | — |
| None | 867 | 14 | 1.6% | 35x |
| La2 | 184 | 1 | 0.54% | 12x |
| lineage1 | 11,500 | 37 | 0.32% | 6.9x |
| lineage6 | 967 | 1 | 0.1% | 2.2x |
| La1 | 8,278 | 1 | 0.012% | 0.26x |
| lineage4 | 56,543 | 3 | 0.005% | 0.1x |
| lineage2 | 38,960 | 0 | 0 | — |
| lineage3 | 11,187 | 0 | 0 | — |

The overall rate was 0.046%. (`None` excludes `ReconstructedAncestor`, a
placeholder sample this tree carries on its root: it was given that name with
`matUtils mask --rename-samples` when the tree was rerooted. The command
above excludes it with `--exclude`; another tree will have its own such
names, or none.) Read the ratio before explaining the gradient,
because it is easy to get wrong: `total_mutations / total_tips` sums each
branch in the subtree *including the node's own stem*, and no root-to-tip path
is computed anywhere. It is a local measure. It is not distance from the
reference, so the choice of reference and the rerooting of the tree do not
enter -- except as below, which is the case that matters.

Compare each clade's pruned count against its count of terminal branches at or
above the threshold, and split the removals by where the threshold was crossed
-- at the sample's own terminal branch, or at an internal node whose whole
subtree went:

| clade | in tree | median terminal | tips >= 250 | pruned | via terminal | via internal node |
|---|---|---|---|---|---|---|
| lineage4 | 56,543 | 2 | 3 | 3 | 3 | 0 |
| lineage2 | 38,960 | 3 | 0 | 0 | 0 | 0 |
| lineage1 | 11,500 | 31 | 30 | 37 | 21 | 16 |
| lineage3 | 11,187 | 12 | 0 | 0 | 0 | 0 |
| None | 867 | 6 | 10 | 14 | 7 | 7 |
| lineage8 | 3 | 1 | **0** | **3** | 0 | 3 |

Where the pruned count exceeds the long tips, samples went with a subtree.
The internal-node prunings show which of two things drove the crossing:

| child | tips | clade | tip terminals | stem | internal below |
|---|---|---|---|---|---|
| node_3 | 3 | lineage8 | 0, 1, 105 | **863** | 87 |
| node_5009 | 2 | None | 190, 198 | **410** | 0 |
| node_5010 | 5 | None | 6, 15, 528, 584, 909 | 27 | 561 |
| node_8504 | 6 | lineage1.1.2 | 218, 223, 235, 247, 261, 297 | 24 | 53 |
| node_13065 | 2 | lineage1.2.1 | 250, 261 | 95 | 0 |
| node_13084 | 2 | lineage1.2.2.2 | 248, 254 | 12 | 0 |
| node_13100 | 2 | lineage1.2.2.2 | 231, 278 | 20 | 0 |
| node_13666 | 2 | lineage1.2.2.2 | 282, 298 | 9 | 0 |
| node_14547 | 2 | lineage1.2.2.2 | 236, 319 | 13 | 0 |

Lineage 1's internal-node removals are small groups of tips that are each long
-- most of them near or over 250 on their own -- joined by short stems; the
tips drove the crossing, and the result is close to what a per-tip filter
would have done. `node_5010` is the same, carried by its 528-, 584- and
909-mutation tips. Two prunings are different. Lineage 8's three genomes
have terminal branches of 0, 1 and 105, none of them over the threshold, and
it was removed as a *clade*, at an ancestral node whose subtree ratio crossed
250 on the strength of its **stem**: 863 mutations, against 193 in the
subtree below it. `node_5009`, two unlabelled genomes at 190 and 198 inside
the lineage 6 neighbourhood, went the same way on a 410-mutation stem.

So one threshold produces two behaviours. Firing on a terminal branch, it
removes an individually divergent genome, which is roughly the intent. Firing
on an internal node, it removes a small clade on a long stem -- and for a clade
near the root, the stem length is set by where the root was placed. This tree
was rerooted to a putative ancestral node, with a new reference built by
applying the rerooting substitutions, which redistributes mutations among the
branches around the root. **Whether lineage 8 survives `--threshold 250` is
therefore a consequence of the rerooting, not a property of those genomes.**
Deeply nested clades are unaffected by the rerooting, which is why lineage 1's
37 removals need no such caveat; `node_5009` is nested too, so its stem is a
property of those two genomes, but it was still removed as a pair rather than
for either genome's own branch.

The quantitative damage was slight either way: 37 of 11,500 lineage 1 tips is
0.3%, far too few to shift a spectrum built on 164,489 mutations, and the bias
runs against the eventual finding rather than for it, since dropping the most
divergent tips of a clade lowers its `CT|GA`.

And the lineage 8 loss costs nothing in practice: 3 samples could never clear
`--min_mutations` at any threshold, and removing 3 tips changes neither the
topology elsewhere nor the global spectrum. The reason to record it is that it
is the clearest instance of stem-driven pruning, and it was visible only
because lineage 8's own branches are short. Without splitting removals by
terminal versus internal node, `node_5009` would have gone unnoticed among the
unlabelled removals.

This is what to fix if you replace the criterion, and it sets what "relative"
should mean: compare a node's **stem against the subtree below it**, not
against the global distribution. A contaminated genome or a junk small clade
has a long stem *and* short internal branches; a genuinely old lineage has a
long stem with proportionate internal diversity. Lineage 8, with a stem of 863
over 193 mutations across 3 tips, is the second pattern, and a
stem-versus-subtree test keeps it.

So: run the check (the commands are at the top of this section), and read
its first two tables this way. A spread across unrelated clades at roughly the overall rate means the filter
removed bad genomes. Concentration in the basal clades, or an entire small
clade going at once, means it removed biology, and the threshold wants raising
-- or replacing with a clade-relative criterion, flagging a tip only when it
is an outlier against its own clade's divergence rather than against the whole
tree. That would have kept lineage 8 while still catching the 909-mutation
genome sitting among close relatives.

Unlabelled samples deserve a look too. The 14 `None` here were removed at 35x
the overall rate and cannot be placed from the clade file alone, so whether
they were poor genomes or simply unclassified divergent ones is unresolved.

#### Where the pruned samples sat, and what it turned up instead

Follow up the per-clade rates by asking *where* in the tree the removals were,
which the prune report gives for free: its first column is the parent of each
pruned subtree, and that parent survives. The tips under each parent, joined
against the `--sample-clades` table, say what the removed samples' neighbours
were; `prune_clade_check.py --input_tree` lists them (by hand, it is one
`matUtils extract -I <parent> -u` per distinct parent).

On the MTBC tree the 14 unlabelled removals sat under 9 parents, and 8 of the 9
were large ordinary clades -- 2,601 tips of `lineage1.2.1.2.1`, 11,385 of
`lineage3`, 36,716 of `lineage4`, and so on -- with the pruned samples as
isolated oddballs inside normal diversity. That is the pattern the filter is
*for*: it triggers on a spike against near-identical neighbours, and these were
spikes. So the pruning question closes. It also settles the lineage 8 worry in
the other direction: 3 samples could never clear `--min_mutations` at any
threshold, and removing 3 tips changes neither the topology elsewhere nor the
global spectrum. The filter's logic is wrong in principle -- step 2 shows it
fires on a clade's stem, which near the root is set by the rerooting rather
than by the genomes -- and on this tree it cost nothing.

The ninth parent was a 91-tip clade whose tips are **all unlabelled**, which
raised the question of the other 853 unlabelled samples that survive pruning.
The sample table answers it without touching the tree, since `matUtils summary
-s` reports each sample's parent node, and usher numbers nodes in traversal
order so a clade occupies a contiguous range:

| | |
| --- | --- |
| parents holding at least one unlabelled sample | 492 |
| ...that also hold a labelled sample | 0 |
| contiguous unlabelled runs of 5+ samples | 11, covering 857 of 867 |

(Counted on the unpruned tree, as is everything here; the pruned tree gives
480 parents.)

**Do not read the zero as a discovery.** It is close to tautological: lineage
barcodes define monophyletic clades, so a sample is unlabelled exactly when its
whole clade lacks a barcode, which makes unlabelled samples monophyletic by
construction. The zero describes the labelling scheme, not the data. (Said
plainly because the first reading here was that it implied a hidden taxon, and
it does not.)

What does answer the question is the flanking lineage of each run (the six
largest; `prune_clade_check.py` lists all 11):

| size | node range | below | above |
| --- | --- | --- | --- |
| 290 | 65852..66056 | `lineage4.3.4.2` | `lineage4.3.3` |
| 169 | 39797..39944 | `lineage2.2.2` | `lineage3` |
| 142 | 62105..62175 | `lineage4.2.1.1` | `lineage4.2` |
| 92 | 59244..59329 | `lineage4.1.3` | `lineage4.2.2` |
| 56 | 70026..70066 | `lineage4.3.3` | `lineage4.3` |
| 28 | 4172..4195 | `La3` | `lineage6.2.3` |

Every run is embedded between named lineages, most of them lineage 4
sublineages. Nothing sits basal, nothing sits outside the named diversity. So
these are ordinary MTBC genomes in sublineages the barcode set does not name --
an annotation gap, not an outgroup and not contamination. The 633 unlabelled
samples inside splits are contributing normally, and the splits need no
revision on their account.

Worth doing anyway, because it is three commands and it is the check that would
have caught the bad case: a run flanked by nothing, or sitting below the
basal-most named lineage, would be an outgroup inside clades whose spectra you
are comparing.

#### Prune known laboratory strains here too

If you already know which samples are laboratory material — vaccine strains,
lab-adapted reference stocks, experimental-evolution lines — prune them now,
in the same step. Their mutations accumulated in culture, not in hosts, so a
split containing them averages two different processes. Step 11 would find
the ones a split contains, but by then every later step has run on the mixed
tree, and pruning renumbers every node, so steps 3 to 10 would have to be
redone. Pruning here costs nothing extra.

Build the list from the clade labels and keep the pattern explicit, since
which labels mark lab strains depends on the organism and the annotation. On
the MTBC tree it is BCG, the vaccine strain, labelled `La1.2.BCG`:

```
awk -F'\t' '$2 ~ /BCG/ {print $1}' clades_raw.txt > lab_samples.txt
matUtils extract -i tree.pruned.pb.gz -s lab_samples.txt --prune \
    -o tree.pruned.nolab.pb.gz
mv tree.pruned.nolab.pb.gz tree.pruned.pb.gz
```

`--prune` takes no argument: the sample list still goes to `-s/--samples`,
and `--prune` inverts that selection so the listed samples are the ones
removed. Without it the same two arguments would keep *only* the lab
strains, which is the opposite tree and will not announce itself, so check
the output's sample count with `matUtils summary`: it should have fallen by
the number of lab samples the divergence prune had not already removed.
(`-I` selects internal node IDs and `-c` annotated clades; all of these are
easy to confuse.) Keep `lab_samples.txt` and the pattern with the run, so
the prune can be repeated on a newer tree.

A label only finds what was labelled. A lab stock with no distinguishing
label, such as an H37Rv derivative filed under the lineage 4.9 it belongs
to, stays in. Step 11 checks the splits for any lab label the pattern
missed.

*This walkthrough did not prune BCG here*, which is why step 11 finds it.

### 3. Build and apply a region mask

```
annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions resistance_genes.txt --promoter_flank 200 \
    --mask_bed selection_loci.bed --mask_positions selection_pos.txt

annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions alignment_loci.txt --mask_pattern '^(PE|PPE|esx)' \
    --mask_product_pattern '(?i)(transposase|integrase|insertion element|^maturase$)' \
    --mask_pad 150 \
    --mask_bed alignment_loci.bed --mask_positions alignment_pos.txt

sort -un selection_pos.txt alignment_pos.txt \
    | awk '{ print "N" $1 "N" }' > matUtils_mask.txt
```

Two passes, because the two lists want different widening: upstream-only for
promoter mutations at selected loci, symmetric for artifacts that spill past a
feature's edges. The `sort -un` merges and deduplicates the two position lists,
and the `awk` puts them in the `N<position>N` form `matUtils mask -m` reads.
With one pass, `--mask_matutils` writes that form directly and the `awk` is
unnecessary.

```
matUtils mask -i tree.pruned.pb.gz -m matUtils_mask.txt \
    -o tree.pruned.masked.pb.gz
```

Expect this to take a while: `matUtils mask` traverses the tree once per masked
position.

#### Masking comes in rounds

The commands above build the finished mask, but nobody could have written it
in one go. Masking is iterative: build a mask, apply it, re-profile (step 4),
see which sites still recur most, and add what they point to. Each
build-apply-re-profile cycle is a *round*. The MTBC tree took three:

| round | adds | mask | prompted by |
| --- | --- | --- | --- |
| 1 | `tb_resistance_genes.txt` with a 200 bp upstream flank, and the `^(PE\|PPE\|esx)` name pattern | 348 kb, 7.9% | Prior knowledge: loci under selection, and the best-known paralogous families |
| 2 | The first section of `tb_alignment_loci.txt` (mobile elements, prophages, repeats), the transposase / integrase / insertion element product pattern, and `--mask_pad 100` | 469 kb | Round 1's surviving top sites, many of them just outside a feature |
| 3 | The round 3 section of `tb_alignment_loci.txt`, `^maturase$`, and `--mask_pad 150` | 504 kb, 11.4% | Round 2's surviving top sites: synonymous-only hotspots, clusters such as `Rv2019`'s just past the 100 bp pad, the `Rv0071` maturase site, and `Rv2082` under a clade-defining branch |

The commands above are round 3, and following this walkthrough you can run
them directly; you don't need to repeat the history. The earlier rounds are
here because the worked examples compare against them, and because on a
different tree you will need rounds of your own: you can't write an artifact
list until you have seen which sites recur. Step 4's "Three rounds, and when
to stop" has what each round bought and how to tell when you are done. After
round 3 the BCG samples were also pruned (step 12); that is a prune, not a
mask round.

*Worked example.* In round 1, the resistance list resolved to 41 intervals
over 81 kb and the PE/PPE pattern to 151 intervals over 266 kb, 348 kb of the
4,411,532 bp reference together. Round 3's commands, as written above, cover
504 kb.

Bases removed is not the number that matters, though; the next step has it.

### 4. Re-profile, and check what the mask bought

```
profile_tree.py --input_tree tree.pruned.masked.pb.gz \
    --reference NC_012345.1.fna --annotation NC_012345.1.gff3.gz > profile_masked.txt

annotate_sites.py --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --from_profile profile_masked.txt --output surviving_sites.tsv \
    --near_mask selection_loci.bed,alignment_loci.bed
```

`--near_mask` checks each surviving site against the BED files the mask was
built from; "The surviving sites are mostly flank" below says what it is for.

*Worked example.* Judge the mask by recurrence, not by total count. "Usable
SNVs" is every A/C/G/T substitution on every branch, so a site that mutated
on 200 branches counts 200 times; it is the tree's parsimony score. The
round 3 column is what step 3's commands produce:

| | raw | round 1 | round 3 |
| --- | --- | --- | --- |
| usable SNVs | 3,249,554 | 2,972,123 (−8.5%) | 2,846,170 (−12.4%) |
| distinct mutated sites | 1,988,164 | 1,865,580 | 1,807,975 |
| mutations per site (mean) | 1.634 | 1.593 | 1.574 |
| worst single site | 6,102 | 1,989 | 225 |
| top 20 sites' share | 1.24% | 0.50% | 0.12% |

Round 1 masked 7.9% of the genome and removed 8.5% of the mutations, so by
total count it looks like it bought nothing. It bought a great deal: 13 of the raw top 20
sites are gone, including all seven largest, and the worst hotspot fell by
two thirds. Those are the sites that break the independence assumption. By
round 3, 11.4% of the genome masked has removed 12.4% of the mutations and
96% of the worst hotspot.

The reason the total barely moved is worth knowing, because it is easy to
misread. Only 122,584 of the ~348,000 masked positions had ever mutated — 35%,
against 45% genome-wide — because the upstream variant caller had already
excluded most of the PE/PPE families, which were 76% of the mask. So the
masked regions were below average in how many positions varied at all, and far
above average in how often those positions recurred.

#### Three rounds, and when to stop

Step 3 lists what each round added. On the MTBC tree the returns fell off
sharply:

| | raw | round 1 | round 2 | round 3 |
| --- | --- | --- | --- | --- |
| usable SNVs | 3,249,554 | 2,972,123 | 2,886,050 | 2,845,156 |
| worst single site | 6,102 | 1,989 | 318 | 225 |
| top 20 sites' share | 1.24% | 0.50% | 0.16% | 0.12% |
| `CT\|GA` | 0.5011 | 0.5058 | 0.5067 | 0.5080 |
| `TC\|AG` | 0.1605 | 0.1581 | 0.1578 | 0.1575 |

These are the figures from the original runs. The round 3 column counts
1,014 fewer mutations than step 3's commands give, because that run used an
unanchored `maturase`, which also masked `Rv0693` (see below); 2,846,170
above is the anchored pattern's figure. Later steps' worked examples come
from the original round 3 tree.

Round 2 was the one that mattered: the worst hotspot fell 84%, for 2.9% of the
mutations. Round 3 cost another 1.4% to take it 29% further. Note the last two
rows, which are the real test of whether a mask is removing noise or reshaping
the data — the folded spectrum moves by ~0.001 per round and the direction
never reverses. A mask that shifted a class by a percent would be deleting
signal, not artifact.

The decision to run round 3 at all was not made on recurrence. The
permutation null on the round 2 tree had already shown recurrence immaterial:
the sites null peaked at 62.7 against an observed 18,283. It was made on
*interpretation* — `Rv2082` held 11 of the 48 mutations defining `node_6`, a
[tier A](#tiers) clade (step 10), inside 917 bp, and no amount of null
calibration addresses a clade whose defining branch sits on top of a suspect
locus. That is a separate
reason to mask, it arrives only once you have splits to look at, and it is
worth a round on its own. Once a round was happening anyway, the cheap fixes
round 2's survivors pointed to went into it too: the wider pad and the
maturase site, below.

#### The surviving sites are mostly flank, which is a pad problem

Annotating what is left is how you decide there is nothing cheap remaining,
and the commonest leftover is a site just outside a locus that is already
masked. The `distance` column cannot show that, because it is measured to the
nearest gene, masked or not. `--near_mask` measures it to the mask itself: it
adds a `mask_gap` column (bp outside the nearest masked interval: positive
past its right end, negative before its left end, in genome coordinates) and
a `mask_interval` column naming the masked genes, for every site within
`--near_mask_within` bp (default 200). On stderr it summarizes, per masked
interval, how many flank sites it has, how many mutations they carry, and how
far it would have to widen to cover them.

**Why the pad went from 100 to 150.** Round 2's survivors showed `Rv2019`'s
recurrence cluster sitting just past the 100 bp pad. To see the same thing on
your own tree, run step 3 with round 3's list but `--mask_pad 100` and
without `^maturase$`. Six of the top 20 surviving sites are then flanks:

| masked interval | gene | sites | mutations | `mask_gap` |
| --- | --- | --- | --- | --- |
| 2,265,889–2,266,505 | `Rv2019` | 5 | 906 | +3, +5, +6, +45, +48 |
| 887,037–888,736 | `Rv0794c` | 1 | 190 | +38 |

Both genes were on the list and both were masked; the recurrence clusters
simply extend past the pad. `Rv2019` is the clear case: three sites within 3
bp of each other, 103 to 106 bp past the gene's end, holding 585 mutations
between them, all three among the ten largest sites left. Every
gap was under 50 bp, so `--mask_pad 150` covers all six. The cost was 17 kb more
masked sequence, 0.4% of the genome, since the pad widens every feature on
the alignment list, not just these two. Whether a global pad was the right
tool is taken up below, under "Hindsight".

Not every survivor of that run was a flank. The second largest, `79558` with
224 mutations, is codon 25 of `Rv0071`, whose only annotation is the product
"maturase". Maturases are typically encoded by mobile group II introns, the
same class of element as the transposases the product pattern already
caught, but no name or product pattern reached it, and the two maturase locus
tags `tb_alignment_loci.txt` had guessed are not in NC_000962.3. So round 3
added it to `--mask_product_pattern`, at a cost of 1.0 kb with the pad.
The pattern is anchored, `^maturase$`, because the product pattern matches
anywhere in a product: a bare `maturase` also catches `Rv0693`,
"mycofactocin radical SAM maturase MftC", an enzyme rather than a mobile
element. The original round 3 run used the bare form and masked it: 1.5 kb
with the pad, holding 1,014 mutations at 663 sites.

On the round 3 tree, at 150, the same check comes back with the flanks
further out:

- Four sites at 2,266,583–2,266,624, 28 to 69 bp past the padded `Rv2019`
  interval (2,265,839–2,266,555; the gene itself is 2,265,989–2,266,405). The
  recurrence cluster continues past each new boundary.
- `889003`, 217 bp past the padded `Rv0794c` interval -- too far out to be a
  flank, and outside the default `--near_mask_within`.
- Five sites at 4,385,916–4,386,142 in `Rv3900c`, annotated in the H37Rv GFF
  as a possible MTBC-specific genomic island. That is the class of the
  prophages and mobile elements on `tb_alignment_loci.txt`, and it was simply
  not on the list.

**Hindsight: a global pad for a local problem.** The 50 bp that round 3
added to every alignment locus came to 17,556 bp, holding 15,491 mutations
on the pruned tree:

- 1,669 were the `Rv2019` and `Rv0794c` flanks the change was aimed at.
- About 430 were lesser hotspots at `Rv1313c`, `PE29`, `PE6` and `whiB6`.
- The other ~13,400 ran at 0.77/bp, about the tree's background rate
  (0.73/bp).

It also didn't finish the job. `Rv2019`'s recurrence sits 82 to 219 bp past
the gene's end: about 2,200 mutations in one block, back to background by
245 bp. A pad of 150 caught only half of it.

The two widenings answered different evidence. Round 1's survivors were
flanks 1 to 56 bp outside many unrelated features. That is a genome-wide
pattern, and a global pad is the right fix. Round 2's were two specific loci,
and the right fix for a specific locus is a targeted interval. About 245 bp
past `Rv2019` and a few tens of bp past `Rv0794c` would have cost under
0.3 kb instead of 17 kb. Neither choice moves the results much, since the
pad's cost is 0.5% of mutations, all ordinary. So the walkthrough keeps 150
rather than spend another round of masking on it.

The rule, then:

- If the flank survivors sit at many unrelated loci, widen `--mask_pad`.
- If they cluster at one or two, extend those loci by hand.

`tb_alignment_loci.txt` takes names, not coordinates, so a targeted extension
means adding positions to the mask in your reference's coordinates. For
`Rv2019` on NC_000962.3, run this before the `sort -un` in step 3:

```
seq 2266406 2266650 >> alignment_pos.txt
printf 'NC_000962.3\t2266405\t2266650\n' >> alignment_loci.bed
```

The BED line keeps `--near_mask` in step 4 consistent with the mask. Record
the extension next to the commands, because it won't carry over to another
reference.

Either way, don't chase a cluster round after round. The four `Rv2019`
sites left at 150 hold ~600 mutations out of 2.8M, and the whole top 20
holds 0.12%. Each round costs hours, over six on a slower server, for
`matUtils mask` alone. Record the finding, put the targeted extension into
the next rebuild from scratch, and stop.

The signal to stop is this: when the surviving top sites are flanks of loci you
already masked, and the folded spectrum has stopped moving, the mask is done.

If the new top sites in `surviving_sites.tsv` still annotate to resistance
loci, widen `resistance_genes.txt` instead of masking those sites one by one.
Watch in particular for a site where every possible substitution is
synonymous, or where a third codon position recurs hundreds of times:
selection cannot produce that, so it is a mapping artifact, and paralogous
gene families are the usual source.

Node IDs changed when the tree was masked, so compare clades between the two
profiles by sample membership, not by ID.

### 5. Decide whether to fold

Read this off `profile_masked.txt` from step 4. Its GLOBAL SPECTRUM section
has the 12-type and folded spectra; the STRAND ASYMMETRY section, which
`--reference` turns on, has what actually decides the question.

The tempting check doesn't work. Comparing each type with its complement in
the global spectrum (C>T against G>A, and so on) looks like a test for a
strand-asymmetric process, but every count is taken on the reference's forward
strand. That is the leading strand on one replication arm and the lagging
strand on the other, and the template strand for about half the genes. A
process tied to the leading or the transcribed strand pushes the forward-strand
ratio one way on half the genome and the other way on the other half, so it
cancels. Near-equal global pairs are what you'd see with or without one.

*Worked example.* On the round 3 tree the global pairs agree closely: C>T
0.2542 against G>A 0.2537, C>A 0.0601 against G>T 0.0604, and every other
pair within 2%. The STRAND ASYMMETRY section shows what that hides. Its
numbers are log2(first/second), where 0 means symmetric:

| | C>T / G>A | C>G / G>C | T>G / A>C | T>C / A>G |
| --- | --- | --- | --- | --- |
| forward strand leading | +0.559 | +0.349 | +0.283 | −0.095 |
| forward strand lagging | −0.485 | −0.313 | −0.242 | +0.080 |
| + strand genes | +0.132 | +0.032 | −0.194 | −0.162 |
| − strand genes | −0.131 | −0.033 | +0.203 | +0.155 |
| intergenic | +0.023 | +0.020 | +0.020 | +0.015 |

The terminus was estimated from GC skew at 2.04 Mb. On the leading strand,
C>T outnumbers G>A by about 47% (2^0.56). In genes, the coding strand carries a
smaller excess of C>T and deficits of T>C and T>G. Both effects flip sign
cleanly and cancel genome-wide. They also vary by clade: in the leading-strand
orientation, C>T/G>A is +0.34 to +0.39 in the lineage 1 sublineages and
mostly +0.5 to +0.67 in lineages 2, 3 and 4. Its heterogeneity chi2/df is 55,
far above the 1 that sampling noise would give.

So fold. But the reason is not that the genome is free of strand-asymmetric
processes. spectrumSplits.py counts forward-strand types, so a clade with a
stronger leading-strand C>T excess gains forward C>T on one arm and forward G>A
on the other, and the two cancel inside the clade just as they do globally.
That difference is invisible with or without `--fold`. The six degrees of
freedom that folding removes carry mostly noise and artifacts from clades whose
mutations happen to sit more on one arm. Folding doubles every cell count and
takes the degrees of freedom from 11 to 5.

**When not to fold.** The 12-type test sees something the folded one can't
only if there is an asymmetry tied to the reference's forward strand itself,
and only if that asymmetry *differs between clades*. A double-stranded
genome's biology doesn't produce one. A single-stranded genome would, and so
could an artifact of the reference or the alignment in part of the tree.

The FOLD TEST at the end of the STRAND ASYMMETRY section checks for it. It
splits each pair's forward-strand log ratio into three parts:

- F, the part tied to the forward strand.
- R, the replication part.
- T, the transcription part.

F, R and T are the mean, the arm contrast and the gene-strand contrast of the
four arm × gene-strand cells, so each cancels the others. A constant F does
nothing to splits, because both sides of every split share it. So the
verdict rests on how much F varies between clades. That is measured two
ways: as chi2/df around the clades' weighted mean, and as the z of the most
extreme clade against the rest. Each clade's F is the mean of its two arms,
so a clade whose mutations sit more on one arm is not penalised.

The clades are the tree's annotated lineages, each holding the mutations
below it that no smaller annotated lineage claims, and only those big enough
to measure enter. "Big enough" is set by precision, not by a count:
`--strand_max_se` (default 0.15) is the largest standard error allowed on a
clade's log2 ratio for the rarest pair, and the report prints the clade size
that implies for this tree's spectrum. Any lower and the smallest rows are
mostly noise. On the round 3 tree that is 24,006 mutations, with TA|AT the
rarest pair, which admits 26 clades. A larger SE admits more, smaller
clades: the chi2/df already allows for their noise, but the single-clade z
gets more chances to cross its limit by chance. `--strand_min_clade` sets a
fixed size instead.

*Worked example.* On the round 3 tree:

| pair | F | z | R | T | F chi2/df across 26 clades | most extreme clade, z |
| --- | --- | --- | --- | --- | --- | --- |
| C>A / G>T | −0.015 | −2.8 | −0.081 | +0.017 | 1.57 | −2.9 |
| C>G / G>C | +0.018 | +3.4 | +0.345 | −0.031 | 0.78 | −1.9 |
| C>T / G>A | +0.035 | +13.5 | +0.525 | +0.037 | 1.55 | −2.8 |
| T>A / A>T | +0.027 | +1.8 | +0.021 | −0.085 | 1.66 | −2.8 |
| T>C / A>G | −0.007 | −1.6 | −0.076 | −0.146 | 1.19 | −3.2 |
| T>G / A>C | +0.023 | +3.5 | +0.300 | −0.247 | 0.93 | +2.1 |

C>T/G>A has a highly significant F (z 13.5). That is because the leading-arm
skew (+0.56) is a little larger than the lagging-arm one (−0.49), probably an
arm difference in the replication effect or a slightly misplaced terminus.
Neither is a reason not to fold. Its F is the same in every clade, with
chi2/df 1.55, against 55 for the replication part in the per-clade table
above. The verdict is fold.

The defaults flag a pair at chi2/df ≥ 3 (`--fold_max_dispersion`) or a
single clade at |z| ≥ 5 (`--fold_max_z`). Recurrent sites inflate both, and
on this tree pairs with no strand effect at all come in at chi2/df 1.2 to 1.7
and |z| up to about 3. A synthetic tree with an extra 10% of one clade's
G>A turned into forward C>T gave chi2/df 26 and z 13 for that clade. The test
reports the clade, so if it fires, look there for a reference or alignment
problem before you unfold.

Settle this now: it changes both the floor in step 6 and the null in step 7.

### 6. Pick `--min_mutations` and `--max_branch_length`

*`--min_mutations`.* The masked profile reports a rarest 12-type class at
p = 0.00796 (floor 629) and a rarest folded class at p = 0.01608 (floor 311).
Folding, so the binding floor is 311. But the floor is not the binding
constraint here — `smallest_eligible_subtree_tips` is:

| `--min_mutations` | eligible nodes | smallest eligible subtree |
| --- | --- | --- |
| 2,500 | 2,510 | 17 tips |
| 5,000 | 1,411 | 62 tips |
| 10,000 | 840 | 130 tips |

A 17-tip clade reaching 2,500 mutations is running 147 mutations per tip
against a tree-wide 23: that is a divergence outlier, not a clade, and
`--min_mutations 2500` would let it be scored as a "split". Take 10,000. It
still leaves 840 candidates, far more than any real analysis needs, clears the
folded floor 32 times over, and puts a 130-tip minimum under every candidate.
Nothing worth finding is near the boundary — the strongest split on this tree
scored 20,568 on an 11,468-tip clade.

*`--max_branch_length`.* Branch counts run median 4, p90 44, p99 106,
p99.9 177, max 886, and the cap table reads:

| cap | branches excluded | mutations excluded |
| --- | --- | --- |
| 50 | 16,811 | 43.72% |
| 100 | 2,577 | 11.24% |
| 250 | 17 | 0.22% |
| 500 | 3 | 0.08% |

There is no knee in that curve: the distribution is broad rather than
bulk-plus-outliers, so any cap that removes a meaningful number of branches
removes a meaningful fraction of the data. Cap at 250, which deletes 17 absurd
branches and nothing else. A cap of 50 would discard 44% of the tree's
mutations, which on a diverged, sparsely sampled species is real divergence.

Note that 63% of all mutations sit on terminal branches here, with a terminal
p90 of 49. That much private variation is more likely per-genome mapping error
than biology, and error has its own spectrum. It will not manufacture a
clade-shaped split, but it dilutes every real one, so read effect sizes with it
in mind.

### 7. Calibrate `--min_chi`

```
permutation_null.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_mutations 10000 --max_branch_length 250 --mode types \
    --replicates 20 --include_observed --seed 1 --nthreads 8 \
    --output null_types.tsv

permutation_null.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_mutations 10000 --max_branch_length 250 --mode sites \
    --max_splits 25 --min_chi 1 --replicates 20 --include_observed \
    --seed 1 --nthreads 8 --output null_sites_cascade.tsv
```

The `types` run scores one split per replicate, over the whole tree. The
`sites` run follows the greedy cascade 25 splits deep, so the same run also
measures the null for the deeper splits ("The threshold should scale with
subtree size" below). Its first iteration is the whole-tree pass. Each
replicate's permutation depends only on `--seed` and the replicate number,
so those rows are exactly what a `--max_splits 1` run would give, and a
separate one adds nothing. The printed summary reports each replicate's
maximum. That is in practice its first split: on the MTBC tree no replicate
had a deeper split that scored higher. So the summary is the whole-tree
calibration. To see the whole-tree rows alone, run
`awk -F'\t' '$2 == 1' null_sites_cascade.tsv`. In this run `--min_chi 1` only
stops the search: a higher value would cut the null off at exactly the depth
being measured.

Use `--output` for the per-split TSV; the summary goes to stdout, so
redirecting stdout saves the summary and silently discards the table. Pass the
same `--fold`, `--min_mutations` and `--max_branch_length` you settled above —
a threshold calibrated under different settings does not transfer — and the
same `--seed` in both, so the replicates are paired and the only difference
between the runs is the null model. `--include_observed` matters: it scores the
real tree through the identical code path, so the comparison is not against a
number from a different pass.

*Worked example*, the round 3 tree from step 4, folded, at
`--min_mutations 10000` and `--max_branch_length 250`:

| | min | median | max |
|---|---|---|---|
| `types` | 15.5 | 19.8 | 31.5 |
| `sites` | 15.7 | 21.8 | 32.0 |

Observed: **17,826.8**, at `node_5011`, splitting 586,079 mutations below from
2,254,631 above.

Three things to read off that.

**Which null to believe.** Always calibrate against `sites`. It keeps the
clumping that recurrent sites cause, so it can only come out at or above
`types`, and a threshold taken from `types` is exceeded by chance whenever
the two differ. Here they barely differ: `sites` tops out 1.6% above `types`,
and only 1 of its 20 replicates exceeds the whole `types` range. That was not
true on the round 2 tree, where `sites` reached 62.7 against 23.9 for
`types`, 2.6x higher, and 11 of 20 `sites` replicates exceeded the entire
`types` range. Round 2's top recurrent sites were still setting the null.
The round 3 mask removed them.

**Whether to keep masking.** The gap between the two nulls' maxima is what
residual recurrent-site clumping costs. Weigh it against another round of
masking. Here it is 0.5, so another round could buy nothing. On the round 2
tree it was 38.8. Even that was immaterial against an observed maximum
hundreds of times larger, which is why step 4 says round 3 was not justified
by recurrence. Stop masking when the gap between the nulls is small compared
with the observed value, not when the recurrence table looks tidy.

**The threshold.** A 2x margin on the `sites` maximum gives `--min_chi 65`,
the "with 2x margin" line of the printed summary. The observed maximum is
274x that. When the headroom is this large the exact choice does not matter;
when it is not, prefer the larger margin.

One caveat on transferring this number: it is the threshold the *first*
split must beat. Deeper splits in the greedy cascade are scored on smaller
subtrees, where both the observed chi-squared and the null are smaller.
Check the cascade before trusting one threshold to govern every split.

#### The threshold should scale with subtree size

That caveat turns out to matter more than a caveat usually does. The printed
summary can't answer it, because it reports only each replicate's maximum,
which is its first split. The answer is in the cascade TSV, and
`null_by_size.py` reads it:

```
null_by_size.py --null null_sites_cascade.tsv --min_chi 65
```

It groups the null splits by the size of the table each was scored on
(`mutations_below + mutations_above`), giving the whole-tree pass a row of its
own and putting the rest in factor-of-2 bins. It then compares each observed
split with the highest null from subtrees within 1.5x of its size. That is the
same ceiling `split_confidence.py` uses in step 10, so the two agree. With
`--min_chi`, it flags the splits that threshold rejects even though they clear
their size-matched null by `--margin` (default 2x), and any it accepts that
don't clear their null at all. It ends with the hotspot check described at the
end of this section.

*Worked example*, same tree and parameters as above, 500 null splits:

| size of table scored | null splits | null max chi-squared |
|---|---|---|
| 2.84M (whole tree) | 20 | 32.0 |
| 1.46M-2.83M | 243 | 23.4 |
| 839k-1.04M | 4 | 13.1 |
| 354k-586k | 12 | 16.0 |
| 181k-345k | 23 | 14.8 |
| 90k-170k | 26 | 14.2 |
| 45k-88k | 72 | 12.8 |
| 22k-44k | 98 | 14.9 |

The second row is the root's residual: the whole tree minus splits already
carved out of it, still scored as one table. Below about 1M the ceiling
settles at 13 to 16, half the whole-tree figure, because chi-squared scales
with the counts in the table. A single `--min_chi` calibrated on the first
split is therefore too strict everywhere below it: 65 against a real ceiling
of about 14 is more than 4x.

That costs real splits. Because `--min_chi 1` scores everything, the cascade
run also shows what a production run at 65 would reject. Four splits between
25.8 and 64.5 clear their size-matched null by 2.0x to 4.5x. The largest is
`node_12073` at 64.5 against a ceiling of 14.2, which misses the threshold by
0.5. Only one scored split genuinely fails, `node_13090` at 9.0 against 14.9.
Two more, at 18.9 and 19.3, clear their null by less than 1.5x.

Two ways to use this. The simple one is to set `--min_chi` from the
*smallest* subtree you intend to split rather than from the root. On this
tree, 2x a ceiling of about 15 means something near 30 rather than 65. The
better one is for the threshold to scale with the candidate subtree's
mutation count. That is a change to `spectrumSplits.py` rather than a tuning
choice; until it exists, a fixed threshold is a known limitation and worth
stating alongside results.

Two caveats on the size-matching itself. The 1.5x windows overlap, so the
ceilings for nearby sizes are not independent estimates. And the null
subtrees are the permuted tree's own, which do not correspond one-to-one
with the observed subtrees being compared against them — the matching is on
size, not identity.

If one node is the maximum in most `sites` replicates, a single hotspot is
setting your threshold: find and mask it rather than paying for it forever.
`null_by_size.py`'s HOTSPOT CHECK counts how often each node is the
whole-tree maximum, and flags any node that is the maximum in more than half
the replicates (`--hotspot_fraction`). Node IDs are comparable across
replicates because the permutation changes only mutation types, never the
tree. Here no node is the maximum in more than 2 of the 20 replicates, which
is what an absence of such a hotspot looks like.

### 8. Run the search with bootstrap, and check where it saturates

```
spectrumSplits.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_chi 30 --min_mutations 10000 --max_branch_length 250 \
    --bootstrap_splits 100 --nthreads 8 --output_spectrum spectra.tsv \
    2> splits.log
```

`--min_chi 30` is 2x the size-matched null ceiling for deep subtrees (about
15 in step 7), not 2x the whole-tree ceiling (65). A threshold set at the
root would reject splits that clear their own null several times over. The
low threshold costs nothing in rigour, because step 10 scores every accepted
split against the null for its own size and ranks it by that margin.

**One run does both jobs, in order.** It first does the ordinary search: a
`New split found at node_N with x2 X` line in `splits.log` for each accepted
split, then `spectra.tsv` and a `Spectra written to spectra.tsv` line. Only
then does it print `Bootstrapping splits with 100 replicates using 8
threads.` and start the replicates. Each replicate is a full search on
resampled sites and writes its own `bootstrap_N_splits_output.tsv` in the
current directory. The replicates are where the time goes: with 8 threads,
about 12 times the ordinary search. So check saturation on the lines above
that marker as soon as it appears, while the replicates run. If the check
says the threshold should change, kill the run and restart it, and you lose
only the replicates already started. The ordinary search is deterministic, so
a separate run without `--bootstrap_splits` would find the same splits and
buy nothing. Step 10 likewise reads only the lines above the marker.

*Worked example (hypothetical — this run has not been done).* Suppose the
lines above the marker record each accepted split's chi-squared as 17827,
4100, 1900, 640, 310, 120, 60, 33. The run is still finding splits just above
the threshold, so the curve has not bent yet — either raise `--min_chi` and
keep the splits that clear the higher bar, or keep the rest and label them
marginal. A run that instead stopped at 17827, 4100, 1900 with nothing under
1000 has found its own ceiling, and the threshold is not what is limiting it.

One `--min_chi` governs every split, but the chance ceiling depends on
subtree size: 32 for the whole tree and 13 to 16 below about 1M (step 7). So
any single value is a compromise. Set from the root (65), it is too strict
for every split but the first few. Set from deep subtrees (30), it is lenient
at the root, but the splits found there clear 30 by orders of magnitude, and
step 10 checks each against its own size anyway. Whatever the value, a run
whose last split sits a fraction above it has not saturated. It has hit the
bar, and splits below it are going unexamined. `null_by_size.py --min_chi`
in step 7 lists which ones.

#### Then interpret

When the run finishes, `split_drivers.py` reads `spectra.tsv` the way you
would by eye, and lists the splits worth explaining:

```
split_drivers.py --spectra spectra.tsv --splits_log splits.log \
    --bootstrap bootstrap_*_splits_output.tsv \
    --input_tree tree.pruned.masked.pb.gz --max_branch_length 250 \
    --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --highlight repair_genes.txt \
    --output split_drivers.tsv > split_drivers.txt
```

A split driven almost entirely by one mutation class is usually an artifact
or a single over-mutated site; one spread across several classes that
survives the bootstrap is the thing worth explaining. The script tests each
of those per split, and labels it `interpret` only when it passes all three:

- **one class**: the split's row is set against its parent's row (the
  subtree it was found in, from `splits.log`) as a chi-squared table, and the
  statistic is divided among the classes. Flagged when the top class carries
  80% or more (`--one_class_share`).
- **one site**: each of the split's five most recurrent sites is removed from
  both rows in turn. Flagged when one removal takes away half the
  chi-squared or more (`--one_site_loss`).
- **weak support**: fewer than half the bootstrap replicates split the node,
  or one within 10 node numbers of it (`--min_support`, `--neighbourhood`,
  the same window as step 10).

Its `branch_mutations` column lists the mutations on each split's defining
branch, comma-separated, which is the form `annotate_sites.py --positions`
takes. Annotating them with `--highlight repair_genes.txt` is the first place
to look for a mechanism. Do it for the `interpret` splits only. A one-class
or one-site split is more likely an artifact, and its defining branch is the
wrong place to look: the `site` column names the recurrent site to annotate
instead, as step 4 did for the sites that survived the mask. Skip, too, any
branch longer than `--max_branch_length`: the search left it out of the
spectrum, so it cannot have caused the split, and its mutations only label
the lineage. The report ends with a loop that does all of this, filled in
with your file names:

```
awk -F'\t' 'NR>1 && $14=="interpret" && split($15, m, ",") <= 250 {print $1}' split_drivers.tsv |
while read n; do
  echo "== $n" >&2
  annotate_sites.py --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
      --highlight repair_genes.txt --output ${n}_branch.tsv \
      --positions "$(awk -F'\t' -v n=$n '$1==n {print $NF}' split_drivers.tsv)"
done 2> branch_annotation.log
awk '/^== /{n=$2} /HIGHLIGHTED/{print n":"$0}' branch_annotation.log
```

It writes one `node_N_branch.tsv` per split, and the last line lists the
splits whose branch hits a highlighted gene.

**Expect hits by chance.** A defining branch with hundreds of mutations will
hit a gene list covering 1-2% of the genome whatever the mechanism, so a hit
means something only against its expectation. `branch_hits.py` computes it:

```
branch_hits.py --branches node_*_branch.tsv \
    --annotation NC_012345.1.gff3.gz --highlight repair_genes.txt \
    --mask selection_loci.bed,alignment_loci.bed --reference NC_012345.1.fna \
    --output branch_hits.tsv > branch_hits.txt
```

Only hits that change the protein count, by default: a synonymous change
cannot disable a repair gene, and counting it lets a mutational hotspot pass
for a recurrent hit. The chance rate per mutation is the genes' length
outside the mask over the genome's length outside the mask (masked sites
cannot carry a mutation), times the share of the branches' coding mutations
that are nonsynonymous. `--count all` counts every hit instead.

Every P comes with a "by chance" count: how many of the tests made would
reach that P by chance, from each test's own expectation. Read it first.
Near or above 1, a result is what the best of that many draws looks like,
and its rank among the others means nothing. That matters most on short
branches, where any single hit gives P of about the branch's mutations times
the rate, so the per-branch ranking there mostly sorts branches by length.
Three kinds of result can carry weight: a branch with more hits than its
length predicts, one gene hit several times on one branch, and one gene hit
on several unrelated branches. The last is the strongest, because a single
branch cannot separate a mechanism from everything else that defines its
lineage.

*Worked example.* The round 3 search accepted 68 splits below the root. 63
were labelled `interpret`. Three were one class, all driven by C>A|G>T:
`node_7984` (83% of its chi-squared), `node_76478` (85%) and `node_72955`
(80%). Two had weak support: `node_69997` (48%) and `node_56835` (41%). No
split hung on one site. The largest loss from removing one was 18%
(`node_21060`, position 2,544,938), so the mask had done its job. The script
rebuilds each split's mutation total from the tree, and every total matched
`spectra.tsv`. That is the check that the tree and `--max_branch_length` are
the ones the search used. So the screen removes the obvious cases, and step
10 ranks the rest.

Three `interpret` splits had branches longer than 250 and were skipped:
`node_5011` (413 mutations), `node_3594` (833) and `node_4609` (591). The
other 60 were annotated, and `branch_hits.py` reported:

```
  chance rate per mutation: 70,579 / 3,907,786 x 973/1,519 nonsynonymous coding mutations = 1.157%

ALL BRANCHES
  60 branches, 1,718 mutations, 21 nonsynonymous hits, 19.9 expected, P(>=21) = 0.43
  nonsynonymous share: 21 of 40 hits (52%), against 973 of 1,519 (64%) coding mutations on these branches
```

Nothing here beats chance. Repair-gene hits change an amino acid *less* often
than the branches' other coding mutations, and the protein-changing ones
number what the genes' length predicts. The best single results:

| result | observed | expected | P | by chance |
|---|---|---|---|---|
| `polA` on `node_14051` | T186P, R188G | 0.047 | 0.0011 | 0.67 of 3,060 branch-gene pairs |
| `node_13078` (`polA` G795D) and `node_14377` (`mutT4` R48G) | 1 hit in 4 mutations | 0.05 | 0.045 | 1.19 of 60 branches |
| `recA` (`node_3821`, `node_12343`) | 2 branches | 0.66 | 0.14 | 2.88 of 51 genes |

The `polA` pair is the only result below 1, and only just. Its two changes
are 6 bp apart, which is also what a single multi-base event or a local
misalignment looks like, so check it in the reads before anything else.
`node_14377`'s `mutT4` R48G is one hit on a four-mutation branch, which
chance produces about once in a run of this size. What makes it notable is
outside knowledge the script cannot use, and its lineage explains it (see
below). Counting synonymous hits too (`--count all`) gives a different and
misleading picture: `ligC` then looks recurrent on three branches (0.25 by
chance), but all three of its hits are synonymous.

**Check each split's lineage before reading a mechanism into it.** If the tree
carries clade annotations, this is a few commands and it changes the
interpretation completely:

```
matUtils summary -i tree.pruned.masked.pb.gz --sample-clades clades.txt
matUtils extract -i tree.pruned.masked.pb.gz -I node_N -u node_N_samples.txt
grep -Fwf node_N_samples.txt clades.txt | cut -f2 | sort | uniq -c | sort -nr
```

`matUtils extract -I` lists every sample below the node, including those in
splits nested inside it. `spectra.tsv`'s `Number_Tips` counts only what is
left after nested splits are carved out, so the two differ whenever a split
has splits below it.

*Worked example.* Every split carrying a repair-gene result in
`branch_hits.py`, and the two that step 9 points to, turned out to be a
recognised lineage. Three pairs of them are nested (checked with `comm -12`
on the sorted sample files):

| split | clades | samples | `Number_Tips` | hit on its branch |
|---|---|---|---|---|
| `node_14051` | all of lineage 2: 2.1, 2.2, 2.2.1, 2.2.2 | 38,960 | 2,586 | `polA` T186P, R188G |
| `node_14377`, inside `node_14051` | lineage 2.2.1 (Beijing) | 35,044 | 1,966 | `mutT4` R48G |
| `node_12343` | lineage 1.2.2 | 2,132 | 444 | `recA` R540Q |
| `node_13078`, inside `node_12343` | lineage 1.2.2.2 | 565 | 565 | `polA` G795D |
| `node_3821` | lineage 6, with 27 lineage 9 samples | 993 | 209 | `uvrA`, `recA`, `nudC`, `mutY`: one each in 192 mutations |
| `node_6` | the animal-adapted lineages: La1 (77%), La1.2.BCG, La2, La3 | 8,729 | 3,178 | `recF` L158V |
| `node_3594`, inside `node_6` | La3 | 269 | 269 | none counted: its 833-mutation branch is longer than `--max_branch_length` |

Every one of those repair-gene hits is lineage-defining, so in every case
the lineage is the explanation to beat.

The repair-gene reading collapses for `node_14377`. One of only four
mutations on its defining branch is `mutT4` R48G, and it is lineage 2.2.1,
so that is the documented Beijing `mutT4` polymorphism. The variant and the
spectrum are both properties of lineage 2.2.1, and a single branch cannot
separate them. Any deep branch carries hundreds of lineage-defining
mutations, so co-occurrence is guaranteed. What would be evidence is
independent recurrence: the same gene disrupted on an unrelated branch with
the same shift. `branch_hits.py` found none that beats chance.

The `polA` pair, the one result `branch_hits.py` could not put down to
chance, sits on the branch that defines all of lineage 2, so about 39,000
genomes share it. That cuts both ways. A misalignment would not normally
reach every genome in a lineage, but a systematic one would: if lineage 2
carries an indel or rearrangement near `polA` relative to the reference,
every lineage 2 sample mismaps the same way, and the result is exactly this,
two changes 6 bp apart on one deep branch. So the first check is the reads
of a few lineage 2 samples across 1,831,200-1,831,300, or a known lineage 2
indel nearby. And if the pair is real, it is fixed across lineage 2 along
with everything else that defines the lineage, so this branch cannot
separate it from lineage either.

`node_3821` is *M. africanum* lineage 6, with lineage 9 inside it. A
192-mutation defining branch is what a deeply diverged lineage looks like,
and its four protein-changing hits against 2.2 expected are well within
chance. Its spectrum is the lineage's spectrum.

The quality reading collapses too, more thoroughly than step 9 could show.
The clade that broke the divergence correlation — high mutations per tip, but
a spectrum matching its low-divergence parent — is `node_3594`, and all 269
of its samples are La3, a distinct animal-adapted species. Its long branches are sparse sampling of a diverged
species, and its spectrum matches its parent because it belongs to the
animal-adapted group. The anomaly has a phylogenetic explanation, which is
the strongest single sign that these spectra track lineage.

One caveat on `node_6`: 336 of its 8,729 samples are `La1.2.BCG`, a
laboratory-passaged vaccine strain. That is 3.9%, so it will not dominate,
but serial passage is not natural evolution and it is worth excluding before
making a quantitative claim about the animal-adapted spectrum.

Extract the split's samples before doing anything else to the tree. Node IDs
are assigned at read time, so the first mask, prune or collapse renumbers
every node and the ID in `splits.log` stops meaning anything (see "Node IDs
are not stable across trees" above). The sample set is the only durable
handle:

```
matUtils extract -i tree.pruned.masked.pb.gz -I node_14051 -u split_samples.txt
```

`-I` selects by internal node ID, which is what `splits.log` reports. `-c` is
for annotated clades — the lineage labels stored as node annotations in a MAT,
such as Pango names — and will not find a `node_N` ID.

### 9. Check the splits against sequencing quality

A split can be real biology or it can be a quality gradient, and the two look
identical in `spectra.tsv`. On a sparsely sampled tree most mutations sit on
terminal branches — 64% on the tree described here — and a terminal branch
carries everything wrong with one genome: mapping error in repeats,
deamination damage, contamination, assembly artifacts. Those have their own
spectrum, usually enriched for C>T and G>A. So a clade whose genomes are
individually worse shows a shifted spectrum for reasons that have nothing to
do with its mutational processes.

Internal branches cannot carry per-genome error. A mutation is placed on an
internal branch only because the samples below it share it, so one bad genome
cannot reach there. Shared artifacts still can — a repeat that misaligns the
same way in every genome of a lineage — which is why this comes after masking
rather than instead of it.

```
branch_class_spectra.py --input_tree tree.pruned.masked.pb.gz \
    --spectra spectra.tsv --fold --max_branch_length 250 \
    --output branch_class.tsv > branch_class_summary.txt
```

Pass the same `--fold` and `--max_branch_length` as the search, and the
`spectra.tsv` from the run you are checking. Split roots are read from its
`Node_ID` column, and each node's mutations are attributed to the nearest
split root at or above it, which is how `spectrumSplits.py` reports spectra:
once a child split is carved out, the parent's row covers only what is left.
Because node IDs are assigned at read time, the `spectra.tsv` and the tree
must be from the same file; the script warns rather than guessing if the IDs
do not match.

**What to read.** The between-split correlations are the test. If the spectrum
tracks lineage, the ordering across splits survives on internal branches
alone. If it tracks quality, it is strong on terminal branches and fades on
internal ones. The within-split internal-versus-terminal chi-squared is *not*
the test: terminal branches genuinely differ on any real tree, so a large
value there is expected. It matters when it varies a lot between splits, which
says the gradient is uneven across clades.

A single clade can also settle it informally, if one exists with high
divergence per tip but a spectrum resembling its low-divergence parent rather
than the other high-divergence clades. That is one clade rather than a
pattern, so treat it as a reason to run this step, not a substitute for it.

*Worked example*, the 69 rows of the round 3 `spectra.tsv` (68 splits and
the root), from `branch_class_summary.txt`. The correlations with
divergence per tip are the same on both branch classes:

| class | internal r | terminal r | internal rho | terminal rho |
|---|---|---|---|---|
| `CA\|GT` | 0.493 | 0.533 | 0.503 | 0.496 |
| `CG\|GC` | -0.304 | -0.323 | -0.233 | -0.235 |
| `CT\|GA` | **0.595** | **0.578** | 0.587 | 0.579 |
| `TA\|AT` | -0.090 | -0.246 | -0.221 | -0.358 |
| `TC\|AG` | **-0.528** | **-0.547** | -0.616 | -0.660 |
| `TG\|AC` | **-0.624** | **-0.635** | -0.600 | -0.602 |

Over the whole tree, the terminal spectrum is within 3% of the internal one
in every class. The summary's internal-versus-terminal chi-squared of 212.4
is large only because the counts are: 1.0M internal and 1.8M terminal
mutations. Within each clade, internal and terminal spectra differ by a
median of 0.12 to 0.55 percentage points per class, with a single maximum
of 2.31 (`CT|GA` in `node_74367`). So the spectrum is a property of the
clade and not of the branch class: terminal branches are not carrying a
distinct artifact spectrum, even though they hold 64% of the mutations. The
quality explanation is ruled out for this tree.

The between-clade separation therefore stands on internal branches alone,
where per-genome error cannot reach: `TC|AG` runs 0.109 to 0.276 across
clades, a 2.5x spread on counts of 3,000 to 99,000 internal mutations each.
And one clade breaks the divergence correlation and holds up. `node_3594`
has a divergence ratio of 45.8, among the highest, but its internal-only
`TC|AG` of 0.229 matches its low-divergence parent `node_6` (0.250, ratio
16.4), not the 0.109 to 0.157 of 22 of the other 23 splits above a ratio of
40 (the exception, `node_3979`, is at 0.191). Step 8's lineage check says
why: all 269 of its samples are La3, inside the animal-adapted `node_6`.

Note what this step does *not* establish. It rules out per-genome artifacts,
not shared ones: a repeat that misaligns the same way in every genome of a
lineage places its mutations on internal branches too. Masking is what
addresses those, which is why it comes first. A clean result here means the
splits are not a sequencing-quality gradient; it does not by itself make them
biology.

### 10. Rank the splits before interpreting them

`split_confidence.py` reads the files steps 7, 8 and 9 already produced and
answers two questions that fail independently. It runs nothing on the tree.

Two of its inputs need assembling first. The bootstrap replicates, one file
each from step 8, become one `replicate<TAB>node` table; and each split's
samples, from the tree, become one `node<TAB>sample` table:

```
for f in bootstrap_*_splits_output.tsv; do n=${f#bootstrap_}; n=${n%%_*};
  cut -f1 "$f" | grep ^node_ | sed "s/^/$n\t/"; done > bootstrap_by_replicate.tsv

for n in $(cut -f1 spectra.tsv | grep ^node_); do
  matUtils extract -i tree.pruned.masked.pb.gz -I $n -u split_samples.txt >& /dev/null
  sed "s/^/$n\t/" split_samples.txt
done > split_clade_samples.tsv
```

The second loop runs `matUtils` once per split, so it takes a while on a big
tree. Run both before anything renumbers the tree's nodes (see step 8).

```
split_confidence.py --splits_log splits.log --spectra spectra.tsv \
    --null null_sites_cascade.tsv \
    --bootstrap_replicates bootstrap_by_replicate.tsv \
    --clade_samples split_clade_samples.tsv \
    --output split_confidence.tsv
```

Pass `--clade_samples` whenever you have it. Both of the other inputs are
reported per node ID, and node IDs alone do not say which clade contains
which; sample membership does, and getting the nesting wrong distorts every
margin below it. See "Why the sizes need membership" below.

**Margin** is the split's chi-squared over the highest chi-squared the
permutation null reached on subtrees of comparable size. This is what makes it
safe to run the search at a low `--min_chi` and filter afterwards, which is
the right way round: a single threshold cannot serve both ends of the size
range, because chi-squared scales with the counts in the table. On the round
3 tree, `--min_chi 30` accepted 68 splits and *every one* cleared its
size-matched ceiling, by 2.44x at the weakest (`node_56835`) and 556x at the
strongest (`node_5011`).

The script needs the log, not just the TSV, because the TSV does not record
which subtree each split was tested in. The log gives it from the `(sub)tree
node_P` lines that precede each `New split found at`. Note that the log also
contains every bootstrap replicate's identical-looking output, since
replicates print to the same stderr; the script stops reading at the bootstrap
marker, and a parser that does not will invent splits the real tree never had.

#### Why the sizes need membership

The search root in the log is not the size of the table. The search carves
accepted splits out as it goes, so a split was tested against its search root
*minus everything already accepted inside it*, and that shrinks from one
iteration to the next. Two things follow, and both bite.

First, the subtraction has to respect nesting: subtract only the maximal prior
splits, since a prior split sitting inside another prior split left with its
ancestor and must not be removed twice.

Second, the nesting cannot be read off the log. A split accepted late can turn
out to be an *ancestor* of one accepted earlier, because both were found while
searching the same root. On the round 3 tree, 13 of the 68 splits sit inside
a clade other than the subtree they were found in. `node_14051`, all of
lineage 2, was accepted after `node_14377` and turned out to contain it, and
`node_12343` contains `node_13078` and `node_13635`, both accepted before it.
Sample membership is the only durable statement of containment, which is what
`--clade_samples` is for. On this tree the log's own nesting happened to give
the same sizes for all 68 splits, so `--clade_samples` changed nothing. That
is the luck of the search order, not something to rely on.

Getting this wrong is not a rounding error and not conservative in any useful
sense: using the search root's full size inflates the size, and so the null
ceiling, for every *late-found* split — which penalises precisely the deep
splits a low `--min_chi` was meant to reach. On the round 3 tree it
understated `node_54711` by 2.4x (15.9x instead of 38.9x), and `node_49931`
and `node_44527` by 2.0x, which drops both from tier A to B. It moved five
splits in all, three of them down and two up: the null ceiling is not
monotonic in size, so a wrong size can err in either direction. `searched_size()` models the carving and is validated
against the exact table totals the cascade null records in its
`mutations_below` and `mutations_above` columns: 21 of 21 agree exactly. If
you change that function, that check is how you find out.

**Support** is the fraction of `--bootstrap_splits` replicates that recover
the split. Those replicates resample alignment positions, so this asks
whether the split depends on a handful of recurrent sites — which is exactly
the artifact that survives masking. The two measures disagree usefully. A
large margin with low support says the signal is real but rests on few
positions. A small margin with high support says a consistent weak
difference.

#### Tiers

The `tier` column combines the two measures into one label per split:

| tier | margin | support | how to treat it |
|---|---|---|---|
| A | >= 5x | >= 70% | Interpret it. Chance does not come close, and the split does not hinge on a few positions |
| B | >= 5x | < 70% | Interpret with a caveat: the signal is real but rests on few positions |
| B | < 5x | >= 70% | Interpret with a caveat: a weak but consistent difference |
| C | < 5x | < 70% | Do not build on it. Report it, if at all, as marginal |

The thresholds are `--strong_margin` (default 5) and `--strong_support`
(default 70). Support is the neighbourhood figure when
`--bootstrap_replicates` is given and the exact one otherwise, and an input
that is missing never passes: without bootstrap replicates no split can be
tier A. The tiers are a reading aid, not a test. A split at 4.9x and 69% is
tier C, one at 5.0x and 70% is tier A, and both deserve the same caution.
What makes them worth having is that they predict reproducibility: when the
search was repeated on a changed tree in step 12, three of the four splits
that failed to come back were tier C, and the fourth was lost to the
`--min_mutations` floor, which the tiers do not measure.

On the round 3 tree the 68 splits came out **35 tier A, 27 tier B and 6 tier
C**, with correct sizing and wobble-tolerant support. The support correction
matters most: exact node IDs alone give 13 tier A and 22 tier C. Search-root
sizing gives 34 tier A, which looks close, but only because five splits move
in opposite directions. Neither correction is a refinement of the other, so
run both.

#### The placement wobble, and a trap

Replicates often place a real split one node up or down, so counting exact
node IDs understates support, sometimes drastically. Given
`--bootstrap_replicates` the script asks, per replicate, whether *any* node
within `--neighbourhood` was split, which is a fraction of replicates and
cannot exceed 100%.

How much this matters: `node_14377`, the Beijing clade, has a margin of 315x
yet appears under its own ID in only 52% of replicates. Allowing +-10 in node
number takes it to 100%. Half the replicates were splitting the same clade
one node away, and on exact IDs alone it would have been filed as a
large-margin, low-support split -- the artifact pattern -- when it is
nothing of the kind.

Do not get there by summing the per-node counts from the aggregate file. That
double-counts any replicate which split two nearby nodes, and on the round 3
tree those sums ran to 192% for `node_14377` and 152% for `node_36360` —
which does demonstrate the wobble is real, but is an upper bound and not a
support value. The aggregate `uniq -c` file is enough for exact support only;
the per-replicate `bootstrap_by_replicate.tsv` made at the top of this step is
what the neighbourhood calculation needs.

Two caveats on the neighbourhood number. Node-ID adjacency is a proxy for
topological adjacency, not the same thing — usher assigns IDs in traversal
order at read time, so neighbours are usually close but need not be. It can
err either way. It can credit a split with a neighbour that is really a
different clade, so read it as an upper estimate, and treat exact support as
the floor. And it can miss a wobble that lands farther away. `node_4609` is
missing from 40 of the 100 replicates, and in 19 of those the split went to
`node_4691` or `node_4693` instead: 82 to 84 node numbers away, outside the
window, but both inside `node_4609`'s clade (all 300 and all 240 of their
samples). Its 60% is an underestimate.

So a split with low exact support and *no* neighbour picking up the signal
is not necessarily a site-driven artifact. `node_3594`, La3, has 69% on both
measures, and in none of the 31 replicates that lack it was any node within
300 numbers of it split instead, except `node_3821` (lineage 6, a different
clade). It is not the `--min_mutations` floor either: its weighted total
across the replicates that found it stayed between 12,079 and 12,687,
against a floor of 10,000. So La3, a distinct species, is genuinely missing
from about a third of the position resamples. Low support without wobble
says the split rests on fewer positions than its margin suggests. It is
worth a look, not a verdict.

The 6 splits in tier C have neighbourhood support of 41% to 69% on margins of
2.4x to 4.2x, so they are weak on both axes rather than rescued or condemned
by either. Two of them, `node_69997` and `node_56835`, are the two that step
8's `split_drivers.py` flagged for weak support. Unlike round 2's tier C, they
do not cluster in node number.

Note what this step does *not* establish. Both measures ask whether a split
is statistically robust, not whether it is biology. A shared alignment
artifact spanning a whole lineage is present in most position-resamples and
clears any size-matched null; steps 3 and 9 are what address it.

### 11. Check what each split *contains* before claiming a mechanism

Steps 9 and 10 ask whether a split is robust. Neither asks what is inside it,
and on a tree that mixes natural populations with laboratory material that is
the question that decides whether a spectrum means anything.

Run it on every split you would interpret, which means tiers A and B: the
tiers measure robustness, not composition, and a tier A split can still be a
third laboratory strain. It needs no tree access once you have
`split_clade_samples.tsv` from step 10 and a `--sample-clades` table
(`matUtils summary -i tree.pruned.masked.pb.gz --sample-clades clades.txt`):

```
split_composition.py --clade_samples split_clade_samples.tsv \
    --clades clades.txt --confidence split_confidence.tsv --lab_pattern BCG \
    --output split_composition.tsv > split_composition.txt
```

For each split it counts the clades of the samples the split's *own*
spectrum covers: its samples minus those of every split nested inside it,
since `spectrumSplits.py` carves nested splits out of their parent's row.
`matUtils extract -I` lists every sample below a node, so counting a parent
that way borrows its children's clades: `node_14051`, all of lineage 2,
would report 36,260 lineage 2.2.1 (Beijing) samples, but 35,044 of them
belong to `node_14377`, nested inside it. Two flags pick out what needs reading:

- **lab**: a clade matching `--lab_pattern` holds 1% or more of the split
  (`--lab_share`). Part of its spectrum accumulated in culture rather than
  in hosts. This is the one that calls for action. There is no default
  pattern, because which labels mark laboratory strains depends on the
  organism and on how the tree was annotated: for MTBC it is `BCG`, the
  vaccine strain. Without it the flag is not tested.
- **mixed**: 20% or more of the split lies outside the top clade's own
  lineage (`--mixed_share`), judged by the dotted names, so a sublineage or
  an ancestor label counts as the same lineage. A deep split often is
  mixed: what remains after its sub-splits are carved out is several sister
  sublineages. That is not a defect, but it says the spectrum is an average
  over them, and a claim about one sublineage cannot rest on it.

*Worked example.* Of the 62 tier A and B splits on the round 3 tree, 11 are
mixed and one is lab. The mixed ones are mostly deep: `node_57318` is
lineage 4.3.3 and 4.3.4 sublineages with 67% outside the top one, and
`node_6`, the animal-adapted group, is La1.8.1, La1.7.1, La2 and others.
The lab one is `node_10`:

```
  node     tier  samples  residual  clades  other  top 3 clades                            flags
  node_10  A     913      913       4       13%    La1.2 49%, La1.2.BCG 37%, La1.3 13%     lab
```

To look at one split in detail, list all its clades from the same files
(this counts every sample below the node, nested splits included):

```
awk -F'\t' '$1=="node_N"{print $2}' split_clade_samples.tsv | sort > n.txt
grep -Fwf n.txt clades.txt | cut -f2 | sort | uniq -c | sort -nr
```

On the round 3 tree `node_10` is 913 tips, all *M. bovis* (La1) — and 336 of
them (36.8%) are `La1.2.BCG`, the vaccine strain serially passaged in vitro since
1908. Every BCG genome in the whole tree falls inside this one split and
nowhere else. Its mutations accumulated under lab propagation and relaxed
selection, not in a host, so `node_10`'s distinctive spectrum (the highest
`TC|AG` and lowest `CA|GT` in the set) is measured on a clade that is over a
third artificial.

**The split's own absence at a finer level is not reassurance.** BCG's residual
is roughly 3,760 mutations (913 tips carry 10,190, so 11.2 per tip), well under
`--min_mutations 10000`. The search could not have separated BCG however
different it is, so only membership could have found this. Any `--min_mutations`
floor creates a band of clades that are invisible to the search rather than
tested and rejected by it, and lab strains, outbreak clusters and recently
sequenced collections all tend to sit in it.

The fix is to prune, not to caveat afterwards — a spectrum computed over
mixed material cannot be decomposed later. The place for that is step 2,
when the lab strains are known in advance. Then this step is the check that
the prune caught them: with the same `--lab_pattern`, the `lab` flag should
not fire. If it fires anyway, for a strain you did not know about or one
labelled differently, prune it as in step 2 and rerun from step 3. Do that
in the same pass as a re-mask if one is pending, since both renumber every
node (see "Node IDs are not stable across trees"), and the rerun is then
also step 12's replicate.

#### Size the confound by mutations, not by tips

Before pruning anything, measure what it is worth, and measure it the way the
test does. `spectrumSplits.py` builds its tables from mutations; tip counts
appear nowhere in the chi-squared. The two can differ by an order of magnitude.

On the MTBC tree the BCG prune removed 336 of 913 tips from `node_10` -- 36.8%
-- and **535 mutations** out of 2,845,156. That is 1.6 mutations per pruned tip
against a tree-wide terminal median of 2, because repeated sequencing of a
near-clonal vaccine stock produces tips that are nearly identical to each other
and to their relatives. Estimating instead from the clade's average rate of
11.2 mutations/tip gives ~3,760, seven times too high. So BCG was about 5% of
`node_10`'s mutations, not 37% of its signal, and the clade's spectrum was
never meaningfully driven by laboratory material.

The honest way to size it is the profile delta: prune, re-profile, subtract.
That is one cheap run and it is exact. A rate multiplication is not, and it
fails hardest exactly where replicate sequencing has inflated a tip count --
which is the case that made you suspicious in the first place.

The prune is still worth making, because it makes every tip-weighted statement
about the clade true and removes a question that would otherwise need a caveat
in every later claim. Just do not expect it to move a spectrum.

More generally: treat the top of the per-clade composition table as part of
reading any split. A clade that is one named sublineage plus one named lab
derivative is not a population.

### 12. Replicate the search on a changed tree, and match by membership

A mask round or a prune gives you a second, independent search over nearly the
same data. That is the strongest check available on this pipeline, and it is
free once you have the sample sets, because node IDs are gone and membership is
all that carries across.

Any change that leaves the data nearly the same will do: a further mask
round, a prune that step 11 called for, or, when neither is pending, a
deliberate one, such as dropping a random tenth of the samples (`matUtils
summary -s`, then `shuf` the names, then `matUtils extract -s ... --prune`).
The point is reproducibility under a change, not the particular change. The
example below uses a mask round plus a BCG prune, because that is the change
this tree went through; the `_nobcg` file names are that run's.

```
split_confidence.py --splits_log splits_nobcg.log --spectra spectra_nobcg.tsv \
    --null null_sites_cascade_nobcg.tsv \
    --bootstrap_replicates bootstrap_by_replicate_nobcg.tsv \
    --clade_samples split_clade_samples_nobcg.tsv \
    --output split_confidence_nobcg.tsv
```

Then match each old clade to its best new clade by Jaccard on the sample sets,
after removing any samples the prune deleted so both are measured on the same
universe. On the MTBC tree, the round-3 mask plus the BCG prune gave 59 splits
against the previous 62, and **58 of the 62 matched at Jaccard 1.000** -- 56 of
them exactly, two off by a handful of tips. Tier counts went 36/20/6 to
33/22/4. That is the result that says the search is finding structure rather
than fitting the artifact landscape, and no amount of null calibration or
bootstrapping could have established it.

#### What failed to replicate is what the tiering predicted

Three of the four lost splits -- `node_74277`, `node_73308`, `node_73473` --
were tier C, and were the three *lowest* margins in the entire old set (2.37x,
2.41x, 2.83x). The three tier C splits that did survive came back tier C. So
the tier was doing its job: it is a prediction about reproducibility, and here
it was tested against an actual replication and held.

#### The fourth loss was the `--min_mutations` floor, not the biology

`node_10` was tier A at 8.9x margin and 100% support, and it is gone. The
reason is arithmetic. Its residual held 10,190 mutations; pruning BCG removed
535 of them, all necessarily inside it since BCG sits nowhere else; 9,655 is
below `--min_mutations 10000`. The clade did not stop being different -- it
stopped being *testable*.

Watch what the search did instead, because it is the reason this is worth
documenting rather than filing as a loss. It carved out the **complement**: a
new 7,356-tip split (`node_286`) inside the same parent, leaving a 768-tip
residual that contains all 577 surviving `node_10` samples and has `node_10`'s
spectrum (`CA|GT` 0.0792 and `TC|AG` 0.2635, against the old 0.0776 and
0.2682). `node_286` shares not one sample with `node_10`.

So the signal is fully intact and has moved from a split to a residual. Two
things follow. A clade dropping off the split list is not evidence it stopped
being distinctive, and you cannot tell the difference by reading node IDs or
split lists -- only by asking where its samples went. And a `--min_mutations`
floor does not merely hide small clades, it can *relocate* a finding from the
split you would report into a residual you might not think to look at.

## Script reference

In walkthrough order.

### `profile_tree.py`

Used in [step 1](#1-profile-the-raw-tree) and again in
[step 4](#4-re-profile-and-check-what-the-mask-bought).

```
profile_tree.py --input_tree tree.pb.gz
profile_tree.py --input_tree tree.pb.gz --node_table nodes.tsv
profile_tree.py --input_tree tree.pb.gz --reference ref.fna --annotation genes.gff3.gz
```

Reports:

| Section | What to do with it |
| --- | --- |
| Tree size, total usable SNVs | Sets the ceiling on how much splitting is possible |
| Global 12-type and strand-folded 6-type spectrum | Gives the rarest-class floor on `--min_mutations` (5 / p_rarest, so expected cell counts stay at or above 5). Don't use it to decide whether to fold: strand asymmetries cancel in it (see [step 5](#5-decide-whether-to-fold)) |
| Per-branch mutation counts, and what each `--max_branch_length` cap would exclude | Pick a cap only if there is a clear outlier tail; on a diverged, sparsely sampled tree long branches are usually real |
| Subtree mutation counts, and eligible candidate nodes at each `--min_mutations` | If only a handful of nodes clear your threshold, the analysis will be coarse whatever `--min_chi` is |
| Mutations:tips ratio, with the `prune_mutation_sample_ratio.py` auto-threshold | Shows what that QC script would prune before you run it |
| Site recurrence (mutations per site, top sites' share) | The independence violation that forces `--min_chi` above any nominal cutoff; `qc/mask_site_splits.py` targets it |
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

`--node_table` writes one row per node (branch and subtree counts, tips,
mutations:tips, per-branch 12-type spectrum) if you want to plot distributions.

### `prune_clade_check.py`

Used in [step 2](#2-prune-samples-with-implausible-divergence). Tabulates
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

Used in [step 3](#3-build-and-apply-a-region-mask) to build the mask, and in
[step 4](#4-re-profile-and-check-what-the-mask-bought) to identify what survived it.

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

### `permutation_null.py`

Used in [step 7](#7-calibrate---min_chi).

Destroys any association between mutation type and position in the tree, then
runs the same search and the same statistic, and reports the chi-squared values
it still finds. Set `--min_chi` comfortably above the maximum.

```
# whole-tree calibration and the cascade null, in one run
permutation_null.py --input_tree tree.pb.gz --min_mutations 500 \
    --replicates 20 --mode sites --include_observed \
    --max_splits 25 --min_chi 1
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
  smaller and not comparable to the first pass. A raised run includes the
  `--max_splits 1` run: each replicate's permutation depends only on `--seed`
  and the replicate number, so its iteration-1 rows are identical, and so is
  in practice the printed summary, which reports each replicate's maximum. Keep
  `--min_chi 1` with a raised `--max_splits`, so the threshold doesn't cut the
  cascade short.
- `--fold` scores the strand-folded 6-type spectrum (6 categories, 5 degrees of
  freedom) instead of 12. This must match how you intend to run
  `spectrumSplits.py`, which takes the same `--fold` flag: the statistic changes
  scale and the null changes shape, so a threshold calibrated one way is
  meaningless the other. See `spectrumSplits/README.md` for when folding is
  appropriate.
- `--include_observed` scores the unpermuted tree through the same code path, so
  the real maximum and the null maximum are directly comparable. If the ratio is
  under about 2x, the strongest real split is close to what chance produces.
- `--nthreads` runs replicates concurrently where `fork` is available, and falls
  back to serial otherwise; running several instances with different `--seed`
  also works.
- `--selftest` checks the chi-squared implementation against a naive reference
  and against `scipy.stats.chi2_contingency`, and exits. Run it once in the
  environment where `spectrumSplits.py` runs, not on a laptop: where scipy is
  not importable the scipy comparison is SKIPPED with a note, so a pass there
  confirms only the naive reference. Confirmed on the MTBC server, where both
  agree on 200 random tables to 1.8e-12.

The statistic is reimplemented in pure Python rather than called through scipy,
so all candidate nodes can be scored in one vectorless pass without scipy's
per-call overhead. It matches `chi2_contingency` exactly for tables larger than
2 x 2 (no continuity correction); columns whose total is zero are skipped, where
scipy would raise.

### `null_by_size.py`

Used in [step 7](#7-calibrate---min_chi).

```
null_by_size.py --null null_sites_cascade.tsv
null_by_size.py --null null_sites_cascade.tsv --min_chi 65 --margin 2
```

Reads a cascade null from `permutation_null.py --max_splits 25 --min_chi 1`
and reports:

| Section | What to do with it |
| --- | --- |
| Null ceiling by subtree size | How far the chance maximum falls below the whole-tree figure as subtrees shrink. A `--min_chi` set from the root is that much too strict for deep splits |
| Observed splits against a size-matched null | Each observed split's margin over the highest null from subtrees within `--null_window` (default 1.5x) of its size. With `--min_chi`, it flags splits that threshold rejects despite a margin of at least `--margin`, and accepted splits with a margin under 1 |
| Hotspot check | Whether one node is the whole-tree null maximum in more than `--hotspot_fraction` (default half) of the replicates. If so, mask the recurrent sites inside it |

The observed rows come from the cascade run's own `--min_chi 1` search, which
accepts every split and so carves subtrees differently from a production run
at a higher threshold. Treat its deep rows as approximate. For the production
run's own splits, `split_confidence.py` models the carving exactly. `ceiling`
and the window come from `split_confidence.py`, so the two scripts give the
same margin for the same size. `--min_nulls` (default 10) marks a ceiling drawn
from too few null splits to trust. `--selftest` checks the binning, the
ceilings, the flags and the hotspot rule on a constructed cascade.

### `split_drivers.py`

Used in [step 8](#8-run-the-search-with-bootstrap-and-check-where-it-saturates).

```
split_drivers.py --spectra spectra.tsv --splits_log splits.log \
    --bootstrap bootstrap_*_splits_output.tsv \
    --input_tree tree.pruned.masked.pb.gz --max_branch_length 250 \
    --output split_drivers.tsv > split_drivers.txt
```

Prints a table of every split, with the splits to interpret listed below it.
`--output` writes the same per split as a TSV:

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
`--max_branch_length`, and is reported. `--annotation`, `--reference` and
`--highlight` are only copied into the loop the report ends with, which
annotates the defining branch of every `interpret` split short enough to be
in its spectrum.
`--selftest` checks the chi-squared arithmetic, the residual walk and the
one-site test on a constructed tree.

### `branch_hits.py`

Used in [step 8](#8-run-the-search-with-bootstrap-and-check-where-it-saturates).

```
branch_hits.py --branches node_*_branch.tsv \
    --annotation NC_012345.1.gff3.gz --highlight repair_genes.txt \
    --mask selection_loci.bed,alignment_loci.bed --reference NC_012345.1.fna \
    --output branch_hits.tsv > branch_hits.txt
```

Reads the `annotate_sites.py --highlight` output for each branch, which
`split_drivers.py`'s loop writes as `node_N_branch.tsv`. `--annotation` and
`--highlight` must be the ones that run used: the gene lengths come from
resolving the same names the same way. A hit is a site inside a named gene;
sites `--highlight_flank` labels as `name@+N` are not counted. Only
nonsynonymous hits, stop gains included, count unless `--count all` is
given, and every result's "by chance" column, not its rank, says whether it
stands out.
`--genome_size` can replace `--reference`.

| Section | What to do with it |
| --- | --- |
| Chance rate | Highlighted genes' length outside the mask over the genome's, times the branches' nonsynonymous share of coding mutations. Check the names reported as not found |
| All branches | Counted hits against expectation summed over every branch, and the nonsynonymous share of all hits beside the share for every coding mutation on the branches |
| Per branch | Each branch with counted hits, by P. "by chance" is how many of the branches would reach that P by chance, from each one's own length, so short branches that cannot reach it do not inflate it |
| One gene, one branch | Genes with several counted hits on one branch, with the mutations. Hits a few bases apart may be one event; check them in the reads before believing them |
| One gene, several branches | Genes hit on two or more branches, against the branches expected to hit them from the gene's length and each branch's; "by chance" corrects across the gene list |

`--output` writes the per-branch numbers as a TSV. `--selftest` checks the
Poisson tails, the chance counts, the hit counting and the mask arithmetic on
constructed files.

### `branch_class_spectra.py`

Used in [step 9](#9-check-the-splits-against-sequencing-quality). Splits
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
  [Node IDs are not stable across trees](#node-ids-are-not-stable-across-trees)).
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

Used in [step 10](#10-rank-the-splits-before-interpreting-them), and again in
step 12. Ranks the accepted splits by two independent measures, margin over a
size-matched null and bootstrap support. It reads only files earlier steps
produce, runs nothing on the tree, and needs nothing beyond the standard
library.

```
split_confidence.py --splits_log splits.log --spectra spectra.tsv \
    --null null_sites_cascade.tsv \
    --bootstrap_replicates bootstrap_by_replicate.tsv \
    --clade_samples split_clade_samples.tsv \
    --output split_confidence.tsv
```

Inputs:

- `--splits_log`: the search's stderr. Bootstrap replicates print to the same
  stream; only the lines before the bootstrap marker are read.
- `--spectra`: the `--output_spectrum` TSV from the same run.
- `--null`: the cascade null from `permutation_null.py --max_splits`, as in
  [step 7](#7-calibrate---min_chi).
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
- `--neighbourhood` (default 10): the node-number window within which a
  bootstrap replicate's split counts as the same split (0 disables).
- `--strong_margin` (default 5) and `--strong_support` (default 70%): tier A
  needs both, tier B one, tier C neither. An input that is missing counts as a
  failure.

Output columns: `node`, `search_root`, `clade_parent`, `chi2`,
`tested_in_mutations`, `null_ceiling`, `null_replicates`, `margin`,
`support_exact`, `support_neighbourhood`, `tier`. `--selftest` runs the
self-tests and exits.

### `split_composition.py`

Used in [step 11](#11-check-what-each-split-contains-before-claiming-a-mechanism).

```
split_composition.py --clade_samples split_clade_samples.tsv \
    --clades clades.txt --confidence split_confidence.tsv --lab_pattern BCG \
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
dotted names), the `--top` (default 3) clades with shares, and the flags
`mixed` (`other` >= `--mixed_share`, default 0.2) and `lab` (clades matching
`--lab_pattern` >= `--lab_share`, default 0.01). `--lab_pattern` has no
default: give the regex for this tree's laboratory-strain labels, such as
`BCG` for MTBC, or the lab flag is not tested. `--tiers` (default `A,B`)
picks which tiers to list; `--tiers all` also lists tier C and the root, and
works without `--confidence`. `--output` writes the same per split as a TSV, with the top
clades as `clade:count`. `--selftest` checks the residual, the lineage
naming and the flags on a constructed tree.

## Gene lists

Three lists of *M. tuberculosis* loci ship with these scripts, for
`annotate_sites.py --regions` and `--highlight`. Two are masked and one is
deliberately not.

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

This list carries gene names and **no locus tags at all** — the Rv numbers are
in comments as hints to check. A name missing from your annotation is reported
as NOT FOUND, which is harmless, whereas a misremembered locus tag would resolve
to the wrong gene and mislead you silently. Verify any you rely on against
Mycobrowser or the NCBI H37Rv annotation.

### `tb_alignment_loci.txt`

`tb_alignment_loci.txt` is also masked, but for the opposite reason to
`tb_resistance_genes.txt`. Resistance loci recur because the same mutation is
selected repeatedly — the variants are real. These loci recur because short
reads cannot be placed unambiguously: paralogous families (PE/PPE, `esx`, the
`Rv1318c`-`Rv1320c` adenylate cyclases, `lldD1`/`lldD2`), repeats (IS6110, the
13E12 family) and prophages (phiRv1, phiRv2). There the variants are largely
not real at all.

Keeping the two lists apart keeps each auditable, and they want different
widening, so [step 3](#3-build-and-apply-a-region-mask) resolves them in two
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
