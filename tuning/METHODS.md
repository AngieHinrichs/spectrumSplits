# Methods

How the tuning workflow works and why: what each step measures, the
statistics behind it, its caveats, and the worked examples from the
*M. tuberculosis* (MTBC) tree it was developed on. [HOWTO.md](HOWTO.md) is
the step-by-step guide for running it on a new tree, and each step here
links to its counterpart there. [SCRIPTS.md](SCRIPTS.md) documents every
script.

## Contents

- [Summary](#summary)
- [Background](#background)
- [About the worked examples](#about-the-worked-examples)
- [Step 1. Profile the raw tree](#step-1-profile-the-raw-tree)
- [Step 2. Prune samples with implausible divergence](#step-2-prune-samples-with-implausible-divergence)
- [Step 3. Build and apply a region mask](#step-3-build-and-apply-a-region-mask)
- [Step 4. Re-profile, and check what the mask bought](#step-4-re-profile-and-check-what-the-mask-bought)
- [Step 5. Decide whether to fold](#step-5-decide-whether-to-fold)
- [Step 6. Pick `--min_mutations` and `--max_branch_length`](#step-6-pick---min_mutations-and---max_branch_length)
- [Step 7. Calibrate `--min_chi`](#step-7-calibrate---min_chi)
- [Step 8. Run the search with bootstrap, and check where it saturates](#step-8-run-the-search-with-bootstrap-and-check-where-it-saturates)
- [Step 9. Check the splits against sequencing quality](#step-9-check-the-splits-against-sequencing-quality)
- [Step 10. Rank the splits before interpreting them](#step-10-rank-the-splits-before-interpreting-them)
- [Step 11. Check what each split contains before claiming a mechanism](#step-11-check-what-each-split-contains-before-claiming-a-mechanism)
- [Step 12. Test a candidate mechanism genotype-first](#step-12-test-a-candidate-mechanism-genotype-first)
- [Step 13. Ask whether the splits' spectra show in their resistance mutations](#step-13-ask-whether-the-splits-spectra-show-in-their-resistance-mutations)
- [Record: the four rounds of tree modification](#record-the-four-rounds-of-tree-modification)
- [Testing the method](#testing-the-method)

## Summary

`spectrumSplits.py` finds clades whose mutation spectrum differs from the
rest of the tree, scoring every candidate with a chi-squared statistic and
keeping those above `--min_chi`. Its parameters can't be copied from one
tree to another, and its output can't be read as mechanism without checks.
The workflow has thirteen steps in four groups:

- **Clean the tree (steps 1 to 4).** Profile the tree, prune genomes with
  implausible divergence and any laboratory strains, and mask regions under
  selection or prone to misalignment. Recurrent sites break the
  independence the chi-squared assumes, so the mask is judged by how much
  recurrence it removes, and built in rounds until the surviving hotspots
  are only flanks of loci already masked.
- **Choose the parameters (steps 5 to 7).** Decide whether to fold
  strand-complementary classes, from a test that separates forward-strand,
  replication and transcription effects. Set `--min_mutations` and
  `--max_branch_length` from the tree's own distributions. Calibrate
  `--min_chi` against permutation nulls that keep recurrent-site clumping.
  Chance chi-squared shrinks with subtree size, so one threshold can't serve
  every split, and the threshold only decides which candidates are admitted.
- **Run and screen (steps 8 to 10).** Run the search with bootstrap, and
  screen each split for one-class or one-site drivers. Rule out a
  sequencing-quality gradient by comparing internal and terminal branches.
  Rank each split by its margin over a size-matched null and its bootstrap
  support.
- **Interpret (steps 11 to 13).** Check what each split contains, and
  whether it differs in trinucleotide context more than random clades of
  its size. Then test candidate genes genotype-first, comparing every hit
  child with its siblings across all independent hits, so that lineage
  can't stand in for mechanism. Finally, ask whether the splits' spectra
  change how often, or with which variants, they acquire resistance.

On the MTBC tree (129,314 samples) this pruned 396 genomes and masked 11.5%
of the genome in four rounds. The search at `--min_chi 30` accepted 66
splits, all clearing their size-matched null (33 tier A, 28 B and 5 C).
Their spectra track lineage, not sequencing quality. Eight differ in
trinucleotide context beyond what random clades show. Among 56 repair genes,
only `dnaE2` is followed by a consistent spectrum shift (more C>T and G>A,
less C>G and G>C), with every hit and with internal hits only. The splits'
spectra change their supply of resistance mutations by at most 11%, and
tilt which resistance variants arise only slightly.

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
  clades near the root. [Step 2](#step-2-prune-samples-with-implausible-divergence) has the worked case.
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

## About the worked examples

The worked examples are measurements from a real *M. tuberculosis* tree
(129,314 samples, 4.4 Mb reference), not invented figures, so the reasoning
can be checked against the reports it came from. They are not
recommendations: another species, or another sampling scheme, will give
different numbers and sometimes different decisions.

The worked numbers in steps 3 to 12 come from the tree after step 2's
pruning, with BCG removed, and after all four masking rounds, unless a
passage says otherwise. That final tree is called the *walkthrough tree*
below: the tree you get by following [HOWTO.md](HOWTO.md). Several passages compare against earlier versions of
the tree, with BCG still in it or before the last masking rounds; the
[record of the four rounds](#record-the-four-rounds-of-tree-modification)
says what each version was.

## Step 1. Profile the raw tree

How-to: [Step 1. Profile the raw tree](HOWTO.md#step-1-profile-the-raw-tree).

Step 2 prunes any subtree whose mutations per tip reach a threshold. Choose
it from the MUTATIONS:TIPS RATIO section of this report rather than
accepting the prune script's automatic one, which is mean + 2 SD. The report
gives both the distribution and how many nodes that automatic threshold
would catch.

*Worked example.* The report shows mean 18.6, median 7, p90 54, p99 114,
p99.9 192, max 909, and an automatic threshold of 71.7 with 10,930 nodes at
or above it. Taking it would prune 10,930 nodes — 5% of the tree — on a
species where a high ratio is often real divergence. 250, well above p99.9,
aims only at the extreme tail.

Read the distribution before choosing. A median of 7 with a p90 of 54 is not a
tight distribution with outliers, it is a broad one, and mean + 2 SD on a broad
distribution cuts into the bulk.

## Step 2. Prune samples with implausible divergence

How-to: [Step 2. Prune samples with implausible divergence](HOWTO.md#step-2-prune-samples-with-implausible-divergence).

`--threshold` is the value chosen in step 1. Keep stdout: the script prints
one line per pruned node with the tips it took, which is the only record of
what was removed. Note that it prunes the whole subtree at the first node
that crosses the threshold and does not recurse into it, so one crossing can
remove several samples.

*Worked example.* At 250 it removed only 60 samples and still cut the worst
terminal branch from 909 mutations to 234.

### Laboratory strains

If you know which samples are laboratory material — vaccine strains,
lab-adapted reference stocks, experimental-evolution lines — prune them now,
in the same step. Their mutations accumulated in culture, not in hosts, so a
split containing them averages two different processes, and a spectrum
computed over mixed material cannot be decomposed afterwards. Pruning later
renumbers every node and means redoing everything from the mask on; pruning
here costs nothing extra.

The search will not do this for you. On the MTBC tree all 336 BCG genomes
sit inside one 913-tip split, 37% of its tips, and BCG's own mutations fall
well under `--min_mutations`, so the search could never have separated them
however different they are. Any `--min_mutations` floor leaves a band of
small clades the search cannot see rather than tests and rejects, and lab
strains sit in it.

*Worked example.* The BCG pattern selects 336 samples, all `La1.2.BCG`, and
none of them was among the divergence prune's 60. `matUtils summary` on the
result reports 128,918 samples, down from 129,254: exactly 336 fewer.

A label only finds what was labelled. A lab stock with no distinguishing
label, such as an H37Rv derivative filed under the lineage 4.9 it belongs
to, stays in, so look through the clade labels for anything else that names
a lab or vaccine strain before settling the pattern.

Do not expect the prune to move spectra much. `spectrumSplits.py` builds its
tables from mutations, and tips appear nowhere in the chi-squared. Pruning
BCG from the MTBC tree removed 336 tips but only **535 mutations** out of
2,845,156: 1.6 per tip, because repeated sequencing of a near-clonal vaccine
stock produces tips nearly identical to each other and to their relatives.
So BCG was about 5% of its split's mutations, not 37%. Estimating from the
clade's average of 11.2 mutations per tip would have said ~3,760, seven
times too many. To size a prune exactly, re-profile and subtract. The prune
is still worth making: it makes every statement about the clade true and
removes a caveat from every later claim.

### What the divergence prune removed, by clade

`prune_clade_check.py` writes, as markdown, every table in this step and the
next one. Its sections match the tables below:

| `prune_check.md` section | Table here |
|---|---|
| 1. Removal rate by clade | clade / in tree / removed / rate / vs overall |
| 2. Terminal branches >= 250 vs. removals | median terminal / tips >= 250 / pruned |
| 3. Prunings at internal nodes | tip terminals / stem / internal below |
| 4. Surviving parents of pruned subtrees | the parents in "Where the pruned samples sat" |
| 5. Unlabelled samples | the unlabelled-parent counts and run table |

It lists every clade, while the tables here leave out the small ones.
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

So: run the check ([HOWTO step 2](HOWTO.md#step-2-prune-samples-with-implausible-divergence) has the commands), and read
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

### Where the pruned samples sat, and what it turned up instead

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

## Step 3. Build and apply a region mask

How-to: [Step 3. Build and apply a region mask](HOWTO.md#step-3-build-and-apply-a-region-mask).

The mask is built from two gene lists, loci under selection and loci that
align badly, each resolved against the annotation by `annotate_sites.py`,
then applied with `matUtils mask`. Two passes, because the two lists want different widening: upstream-only for
promoter mutations at selected loci, symmetric for artifacts that spill past a
feature's edges. The `sort -un` merges and deduplicates the two position lists,
and the `awk` puts them in the `N<position>N` form `matUtils mask -m` reads.
With one pass, `--mask_matutils` writes that form directly and the `awk` is
unnecessary.

Expect this to take a while: `matUtils mask` traverses the tree once per masked
position.

Some artifacts turn up only later, a site or two at a time, from step 4's
flank check or step 7's hotspot check. They are added to the build by
position, with their intervals in the BED file so that step 4's
`--near_mask` agrees with the mask. A mask can't be written in one go: it is
built in rounds, each prompted by what the previous one left. The
[record of the four rounds](#record-the-four-rounds-of-tree-modification)
has what each MTBC round added, what prompted it, and what it bought.

## Step 4. Re-profile, and check what the mask bought

How-to: [Step 4. Re-profile, and check what the mask bought](HOWTO.md#step-4-re-profile-and-check-what-the-mask-bought).

*Worked example.* Judge the mask by recurrence, not by total count. "Usable
SNVs" is every A/C/G/T substitution on every branch, so a site that mutated
on 200 branches counts 200 times; it is the tree's parsimony score. The raw
column is the tree as downloaded, before step 2's prunes; the round 4 column
is what steps 2 and 3's commands produce:

| | raw | round 1 | round 4 |
| --- | --- | --- | --- |
| usable SNVs | 3,249,554 | 2,972,123 (−8.5%) | 2,843,878 (−12.5%) |
| distinct mutated sites | 1,988,164 | 1,865,580 | 1,806,793 |
| mutations per site (mean) | 1.634 | 1.593 | 1.574 |
| worst single site | 6,102 | 1,989 | 225 |
| top 20 sites' share | 1.24% | 0.50% | 0.12% |

Round 1 masked 7.9% of the genome and removed 8.5% of the mutations, so by
total count it looks like it bought nothing. It bought a great deal: 13 of the raw top 20
sites are gone, including all seven largest, and the worst hotspot fell by
two thirds. Those are the sites that break the independence assumption. By
round 4, 11.5% of the genome masked has removed 12.5% of the mutations and
96% of the worst hotspot. (Step 2's prunes are in that 12.5% too, but they
are small: BCG was 535 mutations.)

The reason the total barely moved is worth knowing, because it is easy to
misread. Only 122,584 of the ~348,000 masked positions had ever mutated — 35%,
against 45% genome-wide — because the upstream variant caller had already
excluded most of the PE/PPE families, which were 76% of the mask. So the
masked regions were below average in how many positions varied at all, and far
above average in how often those positions recurred.

### Reading the surviving sites

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

The MTBC record of how the pad was chosen, and what a global pad cost
against a targeted extension, is in
[the record of the four rounds](#the-pad-from-100-to-150-and-the-maturase).
What it led to is the rule:

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
already masked, and the folded spectrum has stopped moving, the mask is done,
as far as the profile can tell. Round 4 shows the limit of that. Its two
sites carry 117 and 109 mutations, below the 145 at which the top 20 ends,
so the recurrence table never listed them. What made them matter was not how
often they mutated but where: all but two of those mutations sit in one clade
of 12,828, enough to set the `sites` null's maximum whenever they drew a rare
class. Only the null sees that, so step 7's hotspot check can reopen the
mask after this step has closed it.

If the new top sites in `surviving_sites.tsv` still annotate to resistance
loci, widen `resistance_genes.txt` instead of masking those sites one by one.
Watch in particular for a site where every possible substitution is
synonymous, or where a third codon position recurs hundreds of times:
selection cannot produce that, so it is a mapping artifact, and paralogous
gene families are the usual source.

Node IDs changed when the tree was masked, so compare clades between the two
profiles by sample membership, not by ID.

## Step 5. Decide whether to fold

How-to: [Step 5. Decide whether to fold](HOWTO.md#step-5-decide-whether-to-fold).

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

*Worked example.* On the walkthrough tree the global pairs agree closely: C>T
0.2542 against G>A 0.2537, C>A 0.0601 against G>T 0.0604, and every other
pair within 2%. The STRAND ASYMMETRY section shows what that hides. Its
numbers are log2(first/second), where 0 means symmetric:

| | C>T / G>A | C>G / G>C | T>G / A>C | T>C / A>G |
| --- | --- | --- | --- | --- |
| forward strand leading | +0.559 | +0.349 | +0.283 | −0.095 |
| forward strand lagging | −0.485 | −0.314 | −0.244 | +0.081 |
| + strand genes | +0.132 | +0.031 | −0.196 | −0.161 |
| − strand genes | −0.131 | −0.033 | +0.203 | +0.155 |
| intergenic | +0.023 | +0.020 | +0.021 | +0.015 |

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
mostly noise. On the walkthrough tree that is 24,006 mutations, with TA|AT the
rarest pair, which admits 26 clades. A larger SE admits more, smaller
clades: the chi2/df already allows for their noise, but the single-clade z
gets more chances to cross its limit by chance. `--strand_min_clade` sets a
fixed size instead.

*Worked example.* On the walkthrough tree:

| pair | F | z | R | T | F chi2/df across 26 clades | most extreme clade, z |
| --- | --- | --- | --- | --- | --- | --- |
| C>A / G>T | −0.015 | −2.7 | −0.081 | +0.018 | 1.61 | −3.0 |
| C>G / G>C | +0.018 | +3.3 | +0.345 | −0.031 | 0.79 | −1.9 |
| C>T / G>A | +0.035 | +13.4 | +0.525 | +0.037 | 1.55 | −2.8 |
| T>A / A>T | +0.027 | +1.8 | +0.021 | −0.085 | 1.67 | −2.8 |
| T>C / A>G | −0.007 | −1.5 | −0.077 | −0.145 | 1.15 | −3.3 |
| T>G / A>C | +0.022 | +3.3 | +0.301 | −0.249 | 0.92 | +2.1 |

C>T/G>A has a highly significant F (z 13.4). That is because the leading-arm
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

Make the decision here, on the masked tree, not earlier. On the raw MTBC tree,
and on the pruned tree before masking, the same test says "consider NOT
folding": TC|AG's F has chi2/df 3.87, with lineage 2.2.1 at z −7.5.
Pruning doesn't change that; masking does, and after masking the verdict is
fold. So what fired was in the masked regions, most likely alignment
artifacts, which is what the test's warning tells you to look for.

The decision then holds for every script that takes `--fold`:
`spectrumSplits.py`, `permutation_null.py`, `branch_class_spectra.py`,
`gene_hit_contrasts.py`, `split_contexts.py`, and `profile_tree.py`'s
context spectrum below. Unfolding carries through to the contexts too. The 96 context classes
fold the strands together, so they can't see the forward-strand signal that
made you unfold. Without `--fold`, the context scripts use the 192
forward-strand classes (12 types x 16 flank pairs) instead. On a
single-stranded genome this matters most: C>U and G>A, say, are different
processes with different contexts.

Settle this now: it changes both the floor in step 6 and the null in step 7.

### Profiling the contexts

The flanking bases come from each mutation's parent genome, reconstructed
down the tree, so **the FASTA must be the root's genome**. On a rerooted tree
that is the rerooted reference, not the published assembly (see "Know how
your tree was rooted"). The report counts the SNVs whose own reference base
disagrees with the reconstructed parent base, and warns above 1%.

Check four things. Disagreements should be near zero. The unknown and
doublet counts say what the context analysis cannot see. No context should
be dominated by one site: the report lists any where a single site carries
20% or more, which would be a hotspot standing in for a context. And the
floors say how many mutations a context test needs.

*Worked example.* On the masked walkthrough tree, using the rerooted
reference, 2,829,239 of 2,843,878 SNVs (99.5%) get a context, and **none**
disagree with the parent base. 310 have an unknown flank, and 14,329 (0.5%)
sit next to another change on the same branch and are left out as doublets.
1,442 contexts (0.05%) would have been miscalled using the reference's
flanks; reconstructing the parent genome costs nothing extra, but on a tree
this shallow it barely matters. No context has a site carrying even 5% of its
count. The context floors explain why the search itself stays at 6 classes:
all 96 at once need about 18,000 mutations on the smaller side, and the
rarest within-class test (`TG|AC`'s 16 contexts) needs about 850 `TG|AC`
mutations, about 11,000 in total. Unfolded, all 192 classes at once would
need about 39,000.

## Step 6. Pick `--min_mutations` and `--max_branch_length`

How-to: [Step 6. Pick `--min_mutations` and `--max_branch_length`](HOWTO.md#step-6-pick---min_mutations-and---max_branch_length).

*`--min_mutations`.* The masked profile reports a rarest 12-type class at
p = 0.00765 (floor 654) and a rarest folded class at p = 0.01541 (floor 325).
Folding, so the binding floor is 325. But the floor is not the binding
constraint here — `smallest_eligible_subtree_tips` is:

| `--min_mutations` | eligible nodes | smallest eligible subtree |
| --- | --- | --- |
| 2,500 | 2,301 | 17 tips |
| 5,000 | 1,284 | 62 tips |
| 10,000 | 777 | 130 tips |

A 17-tip clade reaching 2,500 mutations is running 147 mutations per tip
against a tree-wide 22: that is a divergence outlier, not a clade, and
`--min_mutations 2500` would let it be scored as a "split". Take 10,000. It
still leaves 777 candidates, far more than any real analysis needs, clears the
folded floor 30 times over, and puts a 130-tip minimum under every candidate.
The strongest split (step 8) scored 17,809.5 at `node_4965`, with 585,749
mutations below it.

*`--max_branch_length`.* Branch counts run median 4, p90 43, p99 102,
p99.9 171, max 832, and the cap table reads:

| cap | branches excluded | mutations excluded |
| --- | --- | --- |
| 50 | 15,559 | 41.77% |
| 100 | 2,207 | 10.02% |
| 250 | 14 | 0.19% |
| 500 | 3 | 0.08% |

There is no knee in that curve: the distribution is broad rather than
bulk-plus-outliers, so any cap that removes a meaningful number of branches
removes a meaningful fraction of the data. Cap at 250, which deletes 14 absurd
branches and nothing else. A cap of 50 would discard 42% of the tree's
mutations, which on a diverged, sparsely sampled species is real divergence.

Note that 64% of all mutations sit on terminal branches here, with a terminal
p90 of 47. That much private variation is more likely per-genome mapping error
than biology, and error has its own spectrum. It will not manufacture a
clade-shaped split, but it dilutes every real one, so read effect sizes with it
in mind.

## Step 7. Calibrate `--min_chi`

How-to: [Step 7. Calibrate `--min_chi`](HOWTO.md#step-7-calibrate---min_chi).

One run does all of step 7. It runs two nulls with the same `--seed`, so
their replicates are paired and the only difference between them is the null
model:

- **`types`** scores one split per replicate, over the whole tree.
- **`sites`** follows the greedy cascade 25 splits deep (`--max_splits`), so
  it also measures the null for the deeper splits ("The threshold should
  scale with subtree size" below). Its first iteration is the whole-tree
  pass. Each replicate's permutation depends only on `--seed` and the
  replicate number, so those rows are exactly what a whole-tree-only run
  would give.

The unpermuted tree is scored through the same code path, so the observed
maximum and the null are directly comparable. The cascade stops only at
`--min_chi 1`, the default: a higher value would cut the null off at exactly
the depth being measured.

The report, `null_report.txt` here, has four sections, and the rest of this
step reads them in turn:

1. **PERMUTATION NULL SUMMARY:** each null's replicate maxima side by side,
   and the recommended `--min_chi`.
2. **NULL CEILING BY SUBTREE SIZE:** how the chance maximum falls as
   subtrees shrink.
3. **OBSERVED SPLITS AGAINST A SIZE-MATCHED NULL:** which real splits the
   recommended `--min_chi` would wrongly reject.
4. **HOTSPOT CHECK:** whether one clade is setting the threshold.

Each replicate's maximum is in practice its first split: on the MTBC tree no
replicate had a deeper split that scored higher. So the summary is the
whole-tree calibration. To see the whole-tree rows alone, run
`awk -F'\t' '$2 == 1' null_sites_cascade.tsv`.

The summary recommends `--min_chi` from the 95th percentile of the
replicates' maxima (`--quantile`), not from the largest. The largest of n
replicates grows with n and is set by a single replicate, so a threshold
built on it rises as replicates are added and jumps between runs. A
quantile stays put and only gets more precise. The largest of 20 replicates
sits near the 95th percentile, so 0.95 gives about what the largest of 20
used to; it takes about 100 replicates to resolve, and the summary warns
with fewer. On the MTBC tree, 100 replicates took under two minutes for each
run with 8 threads.

The per-split records go to `null_types.tsv` and `null_sites_cascade.tsv`
(`--prefix`, default `null`); step 10 reads the cascade. The report goes to
stdout. Pass the same `--fold`, `--min_mutations` and `--max_branch_length`
you settled above: a threshold calibrated under different settings does not
transfer.

### The hotspot check

The summary's thresholds are multiples of the single largest null maximum,
so whatever sets that maximum sets the threshold. Under the `sites` null each
site keeps its branches, so a site that mutated many times inside one clade
puts all of those mutations into one class in every replicate. The class
changes from replicate to replicate, so such a clade need not win often; it
wins big when it draws a rare class. The HOTSPOT CHECK flags either
pattern:

- a clade that is the maximum in more than half the replicates
  (`--hotspot_fraction`);
- a clade that holds the two largest maxima (`--top_maxima`).

It counts a node and a child holding nearly all of it as one clade
(`--same_clade`), which would otherwise split the clade's count between two
IDs. For each flagged clade it lists the most recurrent sites and prints an
`annotate_sites.py` command for them. (The rest of the report is read in
"The threshold should scale with subtree size" below. If the check flags a
hotspot and you mask it, rerun `permutation_null.py` and take `--min_chi`
from the new summary.)

On the MTBC tree this check found the two `lppN` sites that became round 4
of the mask; [the record](#round-4-the-lppn-sites) has the details.

### The threshold

*Worked example*, the round 4 tree from step 4, folded, at
`--min_mutations 10000` and `--max_branch_length 250`, 100 replicates:

| | min | median | p95 | max |
|---|---|---|---|---|
| `types` | 14.0 | 18.7 | 25.1 | 29.1 |
| `sites` | 14.1 | 20.7 | 26.1 | 29.3 |

Observed: **17,809.5**, at `node_4965`, splitting 585,749 mutations below from
2,252,671 above.

Three things to read off that.

**Which null to believe.** Always calibrate against `sites`. It keeps the
clumping that recurrent sites cause, so it can only come out at or above
`types`, and a threshold taken from `types` is exceeded by chance whenever
the two differ. Here they barely differ: the `sites` p95 is 4% above the
`types` p95, and only 1 of its 100 replicates exceeds the whole `types`
range. Before round 4, with 20 replicates, `sites` reached 45.4 against
21.4 for `types`, 2.1x higher, set by the `lppN` sites. On the round 2 tree it reached 62.7 against 23.9, 2.6x
higher, with 11 of 20 `sites` replicates above the entire `types` range;
round 2's top recurrent sites were still setting the null, and round 3
removed them.

**Whether to keep masking.** The gap between the two nulls' maxima is what
residual recurrent-site clumping costs. Weigh it against another round of
masking. Here the p95s differ by 1.0, so another round could buy little.
Before round 4 the maxima differed by 24.0, and on the round 2 tree by 38.8. Neither mattered against the
observed maximum, hundreds of times larger, which is why step 4 says round 3
was not justified by recurrence. But the observed maximum is not the only
split the threshold governs. Round 4 took the recommended `--min_chi` from
91 to 53, and the cascade below has five splits in between, at 54.3 to 85.2,
each clearing its size-matched null by 3.2x to 5.4x. A threshold of 91 would
have rejected all five. So weigh the gap against the smaller splits as well
as the largest, not against the recurrence table: stop masking when it is
small next to the splits you would lose.

**The threshold.** A 2x margin on the `sites` p95 gives `--min_chi 53`,
the summary's "Recommended" line (`--margin`, default 2). The observed maximum is
336x that. When the headroom is this large the exact choice does not matter;
when it is not, prefer the larger margin. The reason for a quantile shows
in how much a 20-replicate maximum moves. The `types` null, which the `lppN`
sites cannot touch, gave maxima of 31.5, 28.7, 21.4 and 23.3 over 20
replicates on four versions of this tree: the earlier one, the earlier one
with BCG pruned, and the walkthrough tree before and after round 4. Their
medians were 19.8, 18.9, 18.5 and 18.6. The maxima differ by single high
replicates. On the final tree, the `sites` replicates 1 to 20 are the same in
the 100-replicate run, since each depends only on `--seed` and its number:
their maximum was 25.7, the 100 replicates' p95 is 26.1, and their maximum
is 29.3.

One caveat on transferring this number: it is the threshold the *first*
split must beat. Deeper splits in the greedy cascade are scored on smaller
subtrees, where both the observed chi-squared and the null are smaller.
Check the cascade before trusting one threshold to govern every split.

### The threshold should scale with subtree size

That caveat turns out to matter more than a caveat usually does. The summary
can't answer it, because it reports only each replicate's maximum, which is
its first split. The answer is in the cascade, and the report's middle two
sections read it.

NULL CEILING BY SUBTREE SIZE groups the null splits by the size of the table
each was scored on (`mutations_below + mutations_above`), giving the whole-tree pass a
row of its own and putting the rest in factor-of-2 bins. The second compares
each observed split with a ceiling: the 99th percentile (`--null_quantile`) of the null
splits from subtrees within 1.5x of its size. That is the same ceiling
`split_confidence.py` uses in step 10, so the two agree. It is 0.99 rather
than the summary's 0.95 because each replicate puts several null splits in
a window, so a window held dozens even with 20 replicates, and their maximum
sat near the 99th percentile. On this tree the 20-replicate maxima and the
100-replicate p99s agree within about 20%, so step 10's margins and tiers
keep their meaning. The next section, OBSERVED SPLITS AGAINST A SIZE-MATCHED
NULL, flags the splits the recommended `--min_chi` rejects even though they
clear their size-matched null by `--margin` (default 2x), and any it accepts
that don't clear their null at all. `--flag_min_chi` flags against another
value instead.

*Worked example*, same tree and parameters as above, 2,500 null splits:

| size of table scored | null splits | median | p99 | max | margin at 30 |
|---|---|---|---|---|---|
| 2.84M (whole tree) | 100 | 20.7 | 28.9 | 29.3 | 1.0x |
| 1.43M-2.83M | 1,089 | 13.3 | 21.4 | 28.1 | 1.4x |
| 712k-1.38M | 79 | 11.1 | 19.2 | 19.7 | 1.6x |
| 354k-705k | 100 | 10.2 | 15.8 | 18.2 | 1.9x |
| 177k-351k | 128 | 7.9 | 15.9 | 18.5 | 1.9x |
| 89k-176k | 182 | 7.3 | 15.4 | 19.3 | 1.9x |
| 44k-88k | 300 | 5.9 | 15.3 | 20.0 | 2.0x |
| 22k-44k | 487 | 4.5 | 12.1 | 14.2 | 2.5x |
| 21k-22k | 35 | 3.1 | 8.6 | 8.8 | 3.5x |

The last column is explained below.

The second row is the root's residual: the whole tree minus splits already
carved out of it, still scored as one table. Below about 700k the ceiling
settles at 12 to 16, about half the whole-tree figure, because chi-squared
scales with the counts in the table. A single `--min_chi` calibrated on the
first split is therefore too strict everywhere below it: 53 against a real
ceiling of about 15 is 3.5x. The maxima column shows why the ceiling is a
quantile: in the bins with hundreds of null splits, one replicate's outlier
sets the maximum up to 30% above the p99.

Whether that costs real splits depends on how strict. Because `--min_chi 1`
scores everything, the cascade run also shows what a production run would
reject. At 53 it costs nothing here: the three scored splits below 53 that
clear their size-matched null do so by less than 2x, at 25.8 (1.7x), 19.0
(1.1x) and 18.7 (1.2x), and one genuinely fails, `node_13044` at 9.0
against 12.3. The cost shows at 91, the threshold before round 4. It would
reject five splits between 54.3 and 85.2 that clear their null by 3.2x to
5.4x, among them `node_12027` at 64.1 against a ceiling of 14.9. The
earlier tree's 65 rejected a split of the same size there, `node_12073` at
64.5.

Two ways to use this. The simple one is to set `--min_chi` as low as the
whole-tree null allows, and read off what that buys deeper down. A split of
the whole tree has to beat the whole-tree null outright, so `--min_chi` must
exceed the largest value in the table's whole-tree row: 29.3 here, so 30.
The report prints that number below the table ("Lowest --min_chi that
clears the whole-tree null: 30"), and its last column is that threshold's
margin over each bin's ceiling. Here 30 gives about 2x for every bin below
700k (1.9x to 2.0x, and 2.5x below 44k), 1.6x at 712k to 1.38M, and none to
spare at the whole tree. That is acceptable, because the real
splits there score in the hundreds and thousands, and step 10 scores every
accepted split against the null for its own size. The fixed threshold only
decides which candidates are admitted; the margin decides which are
believed. The better way is for the threshold to scale with the candidate
subtree's mutation count. That is a change to `spectrumSplits.py` rather than a tuning
choice; until it exists, a fixed threshold is a known limitation and worth
stating alongside results.

Two caveats on the size-matching itself. The 1.5x windows overlap, so the
ceilings for nearby sizes are not independent estimates. And the null
subtrees are the permuted tree's own, which do not correspond one-to-one
with the observed subtrees being compared against them — the matching is on
size, not identity.

## Step 8. Run the search with bootstrap, and check where it saturates

How-to: [Step 8. Run the search with bootstrap, and check where it saturates](HOWTO.md#step-8-run-the-search-with-bootstrap-and-check-where-it-saturates).

2x the size-matched ceiling for subtrees under 700k (about 15, step 7), not
2x the whole-tree p95 (53). A threshold set at the
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

A run that stops with nothing near the threshold, say 17,809, 4,100 and
1,900 and nothing under 1,000, has found its own ceiling, and the threshold
is not what is limiting it. A run still finding splits just above the
threshold has not bent, and the threshold is deciding what it finds.

*Worked example.* The lines above the marker record 66 accepted splits, from
17,809.5 at `node_4965` down to 31.0 at `node_56725`. Seven scored over
1,000, but 36 scored under 100, 12 under 50 and 5 under 40, so this run has
not saturated: it is still finding splits at the bar. That is expected here,
and it does not call for a change. 30 is already the lowest value that clears
the whole-tree null, so it cannot go lower, and raising it would discard
small-subtree splits that clear their own, lower ceiling. Step 10 rates those
near-threshold splits against the null for their own size.

One `--min_chi` governs every split, but the chance ceiling depends on
subtree size: 29.3 at most for the whole tree, and about 15 for subtrees
under 700k (step 7). So any single value is a compromise. Set from the root
with a 2x margin (53), it is too strict for most splits. Set at 30, it is
lenient at the root, but the splits found there clear 30 by orders of
magnitude, and step 10 checks each against its own size anyway. Whatever the value, a run
whose last split sits a fraction above it has not saturated. It has hit the
bar, and splits below it are going unexamined. Step 7's report lists which
ones; `permutation_null.py --flag_min_chi` checks a value other than the one
it recommends.

### Interpreting the splits

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

**The highlight list.** The GENE NAME CHECK section checks `--highlight`
the way `check_gene_names.py` checks step 3's mask lists. A gene the
annotation names differently drops out of the highlight without anything
looking wrong, so fix any entry under "To fix" and rerun.

**The defining branches.** The mutations on a split's defining branch are
the first place to look for a mechanism, and the DEFINING-BRANCH ANNOTATION
section annotates them for the `interpret` splits only. A one-class or
one-site split is more likely an artifact, and its defining branch is the
wrong place to look: the `site` column names the recurrent site to annotate
instead, as step 4 did for the sites that survived the mask. Branches longer
than `--max_branch_length` are skipped too, and named: the search left them
out of the spectrum, so they cannot have caused the split, and their
mutations only label the lineage. The section lists each branch that hits a
highlighted gene, and `node_N_branch.tsv` has every site's annotation
(`--branch_dir` sets where; `''` skips them). The `branch_mutations` column
of `split_drivers.tsv` has the same mutations, comma-separated for
`annotate_sites.py --positions`.

**Expect hits by chance.** A defining branch with hundreds of mutations will
hit a gene list covering 1-2% of the genome whatever the mechanism, so a hit
means something only against its expectation. The BRANCH HITS section
computes it. `--mask` takes step 3's mask BEDs, for the chance rate.

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

*Worked example.* The search accepted 66 splits below the root. 61 were
labelled `interpret`. Four were one class: three driven by C>A|G>T,
`node_7938` (83% of its chi-squared), `node_76364` (85%) and `node_72843`
(81%), and one by T>G|A>C, `node_47761` (82%). One had weak support:
`node_69899` (48%). No split hung on one site. The largest loss from removing
one was 9% (`node_61730`, position 2,366,891), so the mask had done its job.
Before step 7's two `lppN` sites were masked, the largest was 18%, at one of
them (2,544,938). The script
rebuilds each split's mutation total from the tree, and every total matched
`spectra.tsv`. That is the check that the tree and `--max_branch_length` are
the ones the search used. So the screen removes the obvious cases, and step
10 ranks the rest.

Three `interpret` splits had branches longer than 250 and were skipped:
`node_4965` (413 mutations), `node_3548` (832) and `node_4563` (591). The
other 58 were annotated, and the BRANCH HITS section reported:

```
  chance rate per mutation: 77,356 / 3,904,767 x 929/1,451 nonsynonymous coding mutations = 1.268%

ALL BRANCHES
  58 branches, 1,643 mutations, 20 nonsynonymous hits, 20.8 expected, P(>=20) = 0.602
  nonsynonymous share: 20 of 39 hits (51%), against 929 of 1,451 (64%) coding mutations on these branches
```

Nothing here beats chance. Repair-gene hits change an amino acid *less* often
than the branches' other coding mutations, and the protein-changing ones
number what the genes' length predicts. The best single results:

| result | observed | expected | P | by chance |
|---|---|---|---|---|
| `polA` on `node_14005` | T186P, R188G | 0.047 | 0.0011 | 0.74 of 3,248 branch-gene pairs |
| `node_13032` (`polA` G795D) and `node_14331` (`mutT4` R48G) | 1 hit in 4 mutations | 0.05 | 0.049 | 1.32 of 58 branches |
| `recA` (`node_3775`, `node_12297`) | 2 branches | 0.63 | 0.13 | 2.87 of 56 genes |

The `polA` pair is the only result below 1, and only just. Its two changes
are 6 bp apart, which is also what a single multi-base event or a local
misalignment looks like, so check it in the reads before anything else.
`node_14331`'s `mutT4` R48G is one hit on a four-mutation branch, which
chance produces about once in a run of this size. What makes it notable is
outside knowledge the script cannot use, and its lineage explains it (see
below). Counting synonymous hits too (`--count all`) gives a different and
misleading picture: `ligC` then looks recurrent on three branches (0.23 by
chance), but all three of its hits are synonymous.

**Check each split's lineage before reading a mechanism into it.** If the tree
carries clade annotations, this is a few commands and it changes the
interpretation completely:

```
matUtils summary -i tree.pruned.masked.pb.gz --sample-clades clades.txt
matUtils extract -i tree.pruned.masked.pb.gz -I node_N -u node_N_samples.txt
grep -Fwf node_N_samples.txt clades.txt | cut -f2 | sort | uniq -c | sort -nr
```

`-I` selects by internal node ID, which is what `splits.log` reports; `-c`
is for annotated clades, such as Pango names, and will not find a `node_N`
ID. `matUtils extract -I` lists every sample below the node, including those
in splits nested inside it. `spectra.tsv`'s `Number_Tips` counts only what is
left after nested splits are carved out, so the two differ whenever a split
has splits below it.

*Worked example.* Every split carrying a repair-gene result in
BRANCH HITS, and the two that step 9 points to, turned out to be a
recognised lineage. Three pairs of them are nested (checked with `comm -12`
on the sorted sample files):

| split | clades | samples | `Number_Tips` | hit on its branch |
|---|---|---|---|---|
| `node_14005` | all of lineage 2: 2.1, 2.2, 2.2.1, 2.2.2 | 38,960 | 2,586 | `polA` T186P, R188G |
| `node_14331`, inside `node_14005` | lineage 2.2.1 (Beijing) | 35,044 | 1,966 | `mutT4` R48G |
| `node_12297` | lineage 1.2.2 | 2,132 | 444 | `recA` R540Q |
| `node_13032`, inside `node_12297` | lineage 1.2.2.2 | 565 | 565 | `polA` G795D |
| `node_3775` | lineage 6, with 27 lineage 9 samples | 993 | 209 | `uvrA`, `recA`, `nudC`, `mutY`: one each in 192 mutations |
| `node_6` | the animal-adapted lineages: La1 (95%), La2, La3 | 8,393 | 768 | `recF` L158V |
| `node_285`, inside `node_6` | La1 only (La1.8.1 84%) | 7,356 | 2,987 | none |
| `node_3548`, inside `node_6` | La3 | 269 | 269 | none counted: its 832-mutation branch is longer than `--max_branch_length` |

These are the same clades as on the earlier run with BCG still in the tree,
with the same sample counts, under new node IDs. The exception is `node_6`,
which lost its 336 BCG samples to step 2's prune. Its row has only 768 tips
because `node_285`, nested inside it, carves most of La1 out of it.

Every one of those repair-gene hits is lineage-defining, so in every case
the lineage is the explanation to beat.

The repair-gene reading collapses for `node_14331`. One of only four
mutations on its defining branch is `mutT4` R48G, and it is lineage 2.2.1,
so that is the documented Beijing `mutT4` polymorphism. The variant and the
spectrum are both properties of lineage 2.2.1, and a single branch cannot
separate them. Any deep branch carries hundreds of lineage-defining
mutations, so co-occurrence is guaranteed. What would be evidence is
independent recurrence: the same gene disrupted on an unrelated branch with
the same shift. BRANCH HITS found none that beats chance.

The `polA` pair, the one result BRANCH HITS could not put down to
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

`node_3775` is *M. africanum* lineage 6, with lineage 9 inside it. A
192-mutation defining branch is what a deeply diverged lineage looks like,
and its four protein-changing hits against 2.4 expected are well within
chance. Its spectrum is the lineage's spectrum.

The quality reading collapses too, more thoroughly than step 9 could show.
The clade that breaks the divergence correlation in step 9 — high mutations
per tip, but a spectrum close to its low-divergence parent's — is
`node_3548`, and all 269 of its samples are La3, a distinct animal-adapted
species. Its long branches are sparse sampling of a diverged species, and its
spectrum is close to its parent's because it belongs to the animal-adapted
group. The anomaly has a phylogenetic explanation, which is
the strongest single sign that these spectra track lineage.

## Step 9. Check the splits against sequencing quality

How-to: [Step 9. Check the splits against sequencing quality](HOWTO.md#step-9-check-the-splits-against-sequencing-quality).

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

*Worked example*, the 67 rows of step 8's `spectra.tsv` (66 splits and the
root), from `branch_class_summary.txt`. The correlations with
divergence per tip are the same on both branch classes:

| class | internal r | terminal r | internal rho | terminal rho |
|---|---|---|---|---|
| `CA\|GT` | 0.492 | 0.529 | 0.501 | 0.490 |
| `CG\|GC` | -0.329 | -0.348 | -0.277 | -0.282 |
| `CT\|GA` | **0.580** | **0.565** | 0.560 | 0.562 |
| `TA\|AT` | -0.078 | -0.240 | -0.199 | -0.342 |
| `TC\|AG` | **-0.506** | **-0.525** | -0.588 | -0.635 |
| `TG\|AC` | **-0.616** | **-0.629** | -0.594 | -0.602 |

Over the whole tree, the terminal spectrum is within 3% of the internal one
in every class. The summary's internal-versus-terminal chi-squared of 210.8
is large only because the counts are: 1.0M internal and 1.8M terminal
mutations. Within each clade, internal and terminal spectra differ by a
median of 0.12 to 0.51 percentage points per class, with a single maximum
of 2.32 (`CT|GA` in `node_74254`). So the spectrum is a property of the
clade and not of the branch class: terminal branches are not carrying a
distinct artifact spectrum, even though they hold 64% of the mutations. The
quality explanation is ruled out for this tree.

The between-clade separation therefore stands on internal branches alone,
where per-genome error cannot reach: `TC|AG` runs 0.109 to 0.276 across
clades, a 2.5x spread on counts of 3,200 to 99,000 internal mutations each.
And one clade breaks the divergence correlation and holds up. `node_3548`
has a divergence ratio of 45.8, among the highest, but its internal-only
`TC|AG` of 0.229 is near its low-divergence parent `node_6` (0.268, ratio
20.3) and La1's `node_285` (0.249, ratio 15.4), not the 0.109 to 0.157 of 23
of the other 24 rows above a ratio of 40 (the exception, `node_3933`, is at
0.191). Step 8's lineage check says
why: all 269 of its samples are La3, inside the animal-adapted `node_6`.

Note what this step does *not* establish. It rules out per-genome artifacts,
not shared ones: a repeat that misaligns the same way in every genome of a
lineage places its mutations on internal branches too. Masking is what
addresses those, which is why it comes first. A clean result here means the
splits are not a sequencing-quality gradient; it does not by itself make them
biology.

## Step 10. Rank the splits before interpreting them

How-to: [Step 10. Rank the splits before interpreting them](HOWTO.md#step-10-rank-the-splits-before-interpreting-them).

`split_confidence.py` reads the files steps 7, 8 and 9 already produced and
answers two questions that fail independently: does the split clear the
chance ceiling for subtrees of its size (**margin**), and do bootstrap
replicates recover it (**support**)? It runs nothing on the tree.

Pass `--clade_samples` whenever you have it. Both of the other inputs are
reported per node ID, and node IDs alone do not say which clade contains
which; sample membership does, and getting the nesting wrong distorts every
margin below it. See "Why the sizes need membership" below.

**Margin** is the split's chi-squared over the highest chi-squared the
permutation null reached on subtrees of comparable size. This is what makes it
safe to run the search at a low `--min_chi` and filter afterwards, which is
the right way round: a single threshold cannot serve both ends of the size
range, because chi-squared scales with the counts in the table. On this
tree, `--min_chi 30` accepted 66 splits and *every one* cleared its
size-matched ceiling, by 2.10x at the weakest (`node_56725`) and 721x at the
strongest (`node_4965`).

The script needs the log, not just the TSV, because the TSV does not record
which subtree each split was tested in. The log gives it from the `(sub)tree
node_P` lines that precede each `New split found at`. Note that the log also
contains every bootstrap replicate's identical-looking output, since
replicates print to the same stderr; the script stops reading at the bootstrap
marker, and a parser that does not will invent splits the real tree never had.

### Why the sizes need membership

The search root in the log is not the size of the table. The search carves
accepted splits out as it goes, so a split was tested against its search root
*minus everything already accepted inside it*, and that shrinks from one
iteration to the next. Two things follow, and both bite.

First, the subtraction has to respect nesting: subtract only the maximal prior
splits, since a prior split sitting inside another prior split left with its
ancestor and must not be removed twice.

Second, the nesting cannot be read off the log. A split accepted late can turn
out to be an *ancestor* of one accepted earlier, because both were found while
searching the same root. On this tree, 12 of the 66 splits sit inside a
clade other than the subtree they were found in. `node_14005`, all of
lineage 2, was accepted after `node_14331` and turned out to contain it, and
`node_12297` contains `node_13032` and `node_13589`, both accepted before it.
Sample membership is the only durable statement of containment, which is what
`--clade_samples` is for. On this tree the log's own nesting happened to give
the same sizes for all 66 splits, so `--clade_samples` changed nothing. That
is the luck of the search order, not something to rely on.

Getting this wrong is not a rounding error and not conservative in any useful
sense: using the search root's full size inflates the size, and so the null
ceiling, for every *late-found* split — which penalises precisely the deep
splits a low `--min_chi` was meant to reach. On this tree it understated
`node_44419` by 1.4x (3.7x instead of 5.3x), and `node_49823`, `node_4563`
and `node_52706` by nearly as much. It moved three splits down a tier:
`node_44419` and `node_4339` from A to B, and `node_7221` from B to C. It can
also move a split up, since the null ceiling is not monotonic in size; on the
earlier run with BCG it moved two up and three down. `searched_size()` models
the carving and is validated
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

### Tiers

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
search was repeated after a re-mask ([Testing the
method](#testing-the-method)), three of the four splits that failed to come
back were tier C, and the fourth was lost to the `--min_mutations` floor,
which the tiers do not measure.

On this tree the 66 splits came out **33 tier A, 28 tier B and 5 tier C**,
with correct sizing and wobble-tolerant support. The support correction
matters most: exact node IDs alone give 16 tier A and 22 tier C. Search-root
sizing gives 31 tier A and 6 tier C. Neither correction is a refinement of
the other, so run both.

### The placement wobble, and a trap

Replicates often place a real split one node up or down, so counting exact
node IDs understates support, sometimes drastically. Given
`--bootstrap_replicates` the script asks, per replicate, whether *any* node
within `--neighbourhood` was split, which is a fraction of replicates and
cannot exceed 100%.

How much this matters: `node_14331`, the Beijing clade, has a margin of 409x
yet appears under its own ID in only 49% of replicates. Allowing +-10 in node
number takes it to 100%. Half the replicates were splitting the same clade
one node away, and on exact IDs alone it would have been filed as a
large-margin, low-support split -- the artifact pattern -- when it is
nothing of the kind.

Do not get there by summing the per-node counts from the aggregate file. That
double-counts any replicate which split two nearby nodes, and on this tree
those sums ran to 188% for `node_14331` and 138% for `node_36254` —
which does demonstrate the wobble is real, but is an upper bound and not a
support value. The aggregate `uniq -c` file is enough for exact support only;
the per-replicate `bootstrap_by_replicate.tsv` that [HOWTO step 10](HOWTO.md#step-10-rank-the-splits-before-interpreting-them)'s first loop makes is
what the neighbourhood calculation needs.

Two caveats on the neighbourhood number. Node-ID adjacency is a proxy for
topological adjacency, not the same thing — usher assigns IDs in traversal
order at read time, so neighbours are usually close but need not be. It can
err either way. It can credit a split with a neighbour that is really a
different clade, so read it as an upper estimate, and treat exact support as
the floor. And it can miss a wobble that lands farther away. `node_4563`,
lineage 5, is missing from 41 of the 100 replicates, and in 20 of those the
split went to `node_4644`, `node_4645` or `node_4647` instead: 81 to 84 node
numbers away, outside the window, but all inside `node_4563`'s clade (all
376, 300 and 240 of their samples). Its 59% is an underestimate.

So a split with low exact support and *no* neighbour picking up the signal
is not necessarily a site-driven artifact. `node_3548`, La3, has 69% on both
measures, and in none of the 31 replicates that lack it was any node within
300 numbers of it split instead, except `node_3775` and `node_3777` inside it
(lineage 6, a different clade). It is not the `--min_mutations` floor either:
its weighted total across the replicates that found it stayed between 12,105
and 12,561, against a floor of 10,000. So La3, a distinct species, is genuinely missing
from about a third of the position resamples. Low support without wobble
says the split rests on fewer positions than its margin suggests. It is
worth a look, not a verdict.

The 5 splits in tier C have neighbourhood support of 48% to 68% on margins of
2.1x to 4.4x, so they are weak on both axes rather than rescued or condemned
by either. One of them, `node_69899`, is the split step 8's
`split_drivers.py` flagged for weak support. All five lie inside
`node_54601`, a lineage 4 clade of about 41,000 samples: three are its
direct children, and `node_69899` and `node_56725` are one level further
down. So the weak end of the search is one region of lineage 4, subdivided
more finely than the data support, rather than marginal splits scattered
across the tree.

Note what this step does *not* establish. Both measures ask whether a split
is statistically robust, not whether it is biology. A shared alignment
artifact spanning a whole lineage is present in most position-resamples and
clears any size-matched null; steps 3 and 9 are what address it.

## Step 11. Check what each split *contains* before claiming a mechanism

How-to: [Step 11. Check what each split contains before claiming a mechanism](HOWTO.md#step-11-check-what-each-split-contains-before-claiming-a-mechanism).

Steps 9 and 10 ask whether a split is robust. Neither asks what is inside it:
one population, or several averaged together. That decides what a claim
about its spectrum can be a claim about.

For each split it counts the clades of the samples the split's *own*
spectrum covers: its samples minus those of every split nested inside it,
since `spectrumSplits.py` carves nested splits out of their parent's row.
`matUtils extract -I` lists every sample below a node, so counting a parent
that way borrows its children's clades: `node_14005`, all of lineage 2,
would report 36,260 lineage 2.2.1 (Beijing) samples, but 35,044 of them
belong to `node_14331`, nested inside it.

A split is flagged **mixed** when 20% or more of it lies outside the top
clade's own lineage (`--mixed_share`), judged by the dotted names, so a
sublineage or an ancestor label counts as the same lineage. A deep split
often is mixed: what remains after its sub-splits are carved out is several
sister sublineages. That is not a defect, but it says the spectrum is an
average over them, and a claim about one sublineage cannot rest on it.

*Worked example.* Of the 61 tier A and B splits, 11 are mixed, mostly deep
ones:

```
  node        tier  samples  residual  clades  other  top 3 clades                                                flags
  node_4965   A     11,468   548       5       38%    lineage1.1.3.3 31%, lineage1.1.3 31%, lineage1.2.1.1 20%    mixed
  node_6      A     8,393    768       5       41%    La1.2 58%, La2 24%, La1.3 16%                               mixed
  node_57208  A     12,388   8,005     5       67%    lineage4.3.3 28%, lineage4.3.4.2 26%, lineage4.3.4.2.1 22%  mixed
```

`node_57208` is what is left of lineage 4.3 once its sublineage splits are
carved out. `node_4965`, the strongest split on the tree, is lineage 1, and
its own row is the 548 samples its many sub-splits leave behind, from both
1.1 and 1.2. And `node_6`, the animal-adapted group, is La1.2, La2 and La1.3
once `node_285` (La1.8.1, La1.7.1 and other La1) and `node_3548` (La3) are
carved out. Their spectra describe those groups, not any one member.

Ten of the 11 are the same clades as on the earlier run with BCG still in
the tree. The eleventh, `node_285`, is new to this run. One clade dropped
out: `node_44419` (lineage 4.1) was mixed then and is 14% other now.

To look at one split in detail, list all its clades from the same files
(this counts every sample below the node, nested splits included):

```
awk -F'\t' '$1=="node_N"{print $2}' split_clade_samples.tsv | sort > n.txt
grep -Fwf n.txt clades.txt | cut -f2 | sort | uniq -c | sort -nr
```

That is also where something step 2 missed would show: a split whose top
clades include a strain you would not call a natural population.

### Within the classes: trinucleotide context

The search scores six classes. Two clades can share those six proportions and
still differ inside a class, for example more C>T at CpG and less elsewhere.
A split found on a class shift may or may not carry a context shift as well.
`split_contexts.py` asks that of every split, after the fact, without
changing the search. Use the tree, `--fold`, `--max_branch_length` and masks
the search used, and the root's genome as `--reference` (on the MTBC tree,
the rerooted `modified_NC_000962.3.fa`). With `--fold` it counts the 96
folded context classes. Without it, for an unfolded search, it counts the 192
forward-strand ones, and the classes below are the 12 types instead of 6.

It tests every split twice: on all branches, and on internal branches only,
which keeps out per-genome error as in
[step 9](#step-9-check-the-splits-against-sequencing-quality). The report has a
section for each run, then a summary of the splits that pass on both. The
TSVs are `split_contexts.tsv` and `split_contexts_internal.tsv` (one row per
split), `split_contexts_96.tsv` and `split_contexts_internal_96.tsv` (one row
per split and context), `split_contexts_table.tsv` (the summary) and
`split_contexts_controls.tsv` (the random clades, for the figures below).

Each split's residual row is set against its clade parent's as a 2 x 96 table
(2 x 192 unfolded), scored with the likelihood-ratio statistic G:

  G = 2 Σ O ln(O / E)

O is each cell's observed count and E its expected count if the two rows
did not differ. This is the G-test of Sokal and Rohlf's *Biometry*, written
G² or called the likelihood-ratio chi-squared or deviance in Agresti's
*Categorical Data Analysis*. It is the likelihood-ratio counterpart of
Pearson's chi-squared, which `spectrumSplits.py` uses, and on large tables
the two nearly agree. G is used here because, unlike Pearson's statistic, it
adds up exactly over nested tables, so the 2 x 96 table's G splits into two
parts that sum to it:

- **G(between)**, the 2 x 6 table of classes: the shift the search found.
- **G(within)**, the sum over the six classes of a 2 x 16 table of that
  class's contexts. It is zero whenever the context shares inside every
  class agree, whatever the class proportions, so it is the test. In
  statistical terms it tests conditional independence: is a mutation's
  context independent of which row it is in, given its class? In a
  log-linear model, that is the context x group | class term.

The p values do not come from the chi-squared distribution G would follow
if mutations were independent, because they are not: recurrent sites and
shared ancestry tie them together. They come from the two nulls below.

Two nulls are needed, and they answer different questions:

- **Branch shuffle (`p`, `q`).** Whole branches are reassigned between the
  split and its comparator, each keeping all its mutations. This asks whether
  the split differs from its comparator. `excess` is G(within) over the
  shuffles' mean, about 1 under the null, so splits of different sizes can be
  compared.
- **Random clades (`ctrl_p`, `ctrl_q`).** Non-split clades within 2-fold of
  each split's size are tested the same way, each against the rest of its
  residual group. `ctrl_p` is the share of them whose excess reaches the
  split's. This asks whether the split differs more than any clade would.
  `ctrl_q` is `ctrl_p` corrected for the number of splits tested
  (Benjamini–Hochberg), as q corrects p, and is what a split must pass.

The second matters. If context differs between sublineages all over the
tree, a split can pass the shuffle for the same reason any clade would.

*Worked example*, 1,000 shuffles per split and 30 controls per split (about
900 per run), in 75 seconds and 1 GB. The all-branch section reads:

```
  66 splits tested; 44 with q < 0.05.
  911 random non-split clades, size-matched to the splits: 239 (26%) with p <= 0.05;
  median excess 1.09, against 1.50 for the splits.
  8 splits have q < 0.05 and ctrl_q <= 0.05 (18 with uncorrected ctrl_p <= 0.05).
  11 splits with q < 0.05 have too few controls in their size window to reach ctrl_q <= 0.05
```

So 26% of clades chosen at random pass the shuffle. In the internal-branch
section the share falls to 11%, so terminal branches add much of it. On this
tree the q values overstate the case on their own; read `ctrl_q`. Without
the correction, 18 splits have `ctrl_p` at or below 0.05, where chance alone
would give about 3.

**Splits that can't be judged.** `ctrl_p` can't fall below 1/(controls + 1),
so a split with few controls in its size window can't reach `ctrl_q` ≤ 0.05
however strong its shift. Each run names those splits and says which fix
applies. Here all 11 are the largest splits, with 31 to 122 controls each,
and every clade of their size on the tree has already been drawn: raising
`--controls` from 30 to 100 adds none. Only a wider `--control_window` would,
at the cost of looser size matching. They include `node_57208`, whose
all-branch excess of 5.0 is the third largest of any split. So for the
largest splits this test is silent, not negative.

The summary takes the splits that pass both nulls on both runs. There are
eight, sorted by excess. This table is `split_contexts_table.tsv`, with the
clades shortened:

| Split | Clade | Excess, all / internal | Contexts, split vs clade parent, share of class |
| --- | --- | --- | --- |
| `node_7246` | lineage 1.1.3.1 | 7.7 / 3.6 | A[T>G]G 0.210 vs 0.024; A[C>G]G 0.120 vs 0.049 |
| `node_6` | La1.2, La2, La1.3 | 5.9 / 3.3 | G[T>C]C 0.237 vs 0.135; C[T>G]G 0.079 vs 0.166 |
| `node_14331` | lineage 2.2.1 | 4.3 / 2.8 | C[T>C]C 0.087 vs 0.043; A[T>C]C 0.112 vs 0.076 |
| `node_76364` | lineage 4.8 | 3.3 / 1.6 | T[C>A]G 0.133 vs 0.070; T[C>A]A 0.065 vs 0.028 |
| `node_3775` | lineage 6 and 9 | 3.0 / 1.9 | G[C>T]C 0.069 vs 0.112; C[C>T]A 0.089 vs 0.062 |
| `node_21970` | lineage 2.2.1 | 2.7 / 1.8 | G[T>G]G 0.189 vs 0.129; G[C>T]A 0.074 vs 0.059 |
| `node_71721` | lineage 4.6.1 | 2.7 / 1.5 | G[T>G]G 0.166 vs 0.103; G[T>G]C 0.185 vs 0.138 |
| `node_3548` | La3 | 2.3 / 1.6 | G[T>C]G 0.137 vs 0.176; G[C>A]C 0.125 vs 0.078 |

The summary also names the splits that pass on only one run: here five
(`node_13032`, `node_4965`, `node_52706`, `node_54608` and `node_57207`)
pass on internal branches only. `--control_alpha 0.1` adds one split to the
table, `node_54608`. `--summary_top N` trims the table to the N with the
lowest all-branch `ctrl_p`, if a tree gives more than you want to read.

None of these contexts comes from a single site. The summary checks each
listed context for one site carrying 20% or more of the split's count
(`--one_site`), and here the largest is 3.9%, `node_21970`'s G[T>G]G; the
`top_site_shares` column has them all. Three are leads in their
own right:

- **`node_7246`.** A[T>G]G is a fifth of the clade's T>G, against 2% in the
  rest of lineage 1.1.3, and it holds on internal branches. A ninefold shift
  in one context across many sites suggests one process with a sharp
  sequence preference. It could also be a shared artifact in that clade's
  sequencing, so check the reads before reading biology into it.
- **`node_76364` and `node_76366`.** C>A at TpC, in all three contexts,
  T[C>A]G, A and C. It rises in `node_76364`'s own row and falls back in
  `node_76366`, the lineage 4.8 split nested inside it. So it belongs to the
  branches between the two nodes, not to lineage 4.8 as a whole.
- **`node_14331`.** Beijing doubles C[T>C]C.

**Figures.** `plot_split_contexts.py` draws three from the same files. It
needs only matplotlib, so it can run outside the `bte` environment:

```
plot_split_contexts.py --prefix split_contexts
```

- `split_contexts_overview.pdf`: excess against split size, one panel per
  run, over the grey cloud of random clades. This is the `ctrl_p` comparison
  drawn out. The summary splits stand clear of the cloud. The crosses on the
  right are the large splits that can't be judged, `node_57208` at excess 5
  among them, where the cloud runs out.
- `split_contexts_profiles.pdf`: one row per summary split, the 96 contexts
  in signature order. Bars are the change in each context's share of its
  class on all branches, dots the same on internal branches only. A class
  shifted in proportion draws flat. `node_7246`'s A[T>G]G stands out in its
  row, and `node_76364`'s rise is all C>A after a T. Triangles mark contexts
  where one site carries 20% or more of the split's count. None of the named
  contexts has one, but several others do, mostly T>A: read those bars as one
  site each.
- `split_contexts_heatmap.pdf`: every split by every context, coloured by
  signed √(G contribution), with cells below nominal significance drawn
  neutral. It shows what the profiles can't: contexts that recur across
  splits, such as G[T>G]G in `node_21970`, `node_71721`, `node_54603` and
  `node_52706`. That is either a process shared across lineages or a shared
  artifact; either is worth checking.

What this cannot tell you is whether a search scoring all the contexts would
find splits this one did not. That needs the search itself to score
contexts.

## Step 12. Test a candidate mechanism genotype-first

How-to: [Step 12. Test a candidate mechanism genotype-first](HOWTO.md#step-12-test-a-candidate-mechanism-genotype-first).

Step 8 showed why a repair-gene hit on a split's defining branch proves
nothing: that branch carries every lineage-defining mutation, so co-occurrence
is guaranteed. The evidence that would count is independent recurrence -- the
same gene hit on unrelated branches, each followed by the same shift. This
step looks for exactly that, starting from the genotype instead of the
splits. One run of `gene_hit_contrasts.py` tests one gene list three ways:
with every hit counted, with internal hits only, and against a control list
drawn by `control_genes.py`.

`--permutations 10000` because the smallest p a permutation test can give is
1 over the number of permutations plus one, and with dozens of genes the bar
for any one of them is far below 0.001: 0.05 over 56 genes is 0.0009. The
default of 1,000 cannot reach it.

It finds every branch with a coding SNV in a listed gene, classifies it
against the codon as it stood in the *parent's* genome (not the reference,
which is wrong wherever an ancestral branch already changed the codon), and
compares what accumulated in the hit child's subtree with what accumulated in
its siblings'.

**There is no sister clade on an UShER tree.** Parents are polytomies, often
with hundreds of children, so the hit child is compared with *all* of its
parent's other children pooled. Each sibling descends from the same ancestral
genome, so each is a valid contrast and the lineage is held fixed within
every comparison. The null is that the hit child is exchangeable with its
siblings. A permutation redraws, for every event, which child carried the hit
and recomputes the pooled statistic. For a binary parent that is a paired
sign-flip; for a polytomy it is the same test with more choices, and the
polytomy's size does not affect its validity.

Two details matter for the null being honest:

- **The draw is size-biased.** A mutation lands on a long branch more often
  than a short one, so real hit branches are longer than random siblings.
  Uniform draws would make every hit child look mutation-rich. Each child is
  drawn in proportion to its top-branch SNV count plus one -- the chance the
  hit would have landed there.
- **The hit mutation itself is removed** before counting. Other hits of the
  same gene are cut out of each other's comparisons, and a hit nested inside
  another is tested as its own event (`nested_under` names the outer one).

`--top_branch` sets whether each child's own branch is counted, which
decides whether terminal hits contribute at all; see the
[script reference](SCRIPTS.md#gene_hit_contrastspy) for the trade-off, and for what
`--selftest` checks.

There is no rate test. A mutator should also raise the mutation count below
the hit, and an earlier version compared tip depth below hit children with
their siblings'. On the MTBC tree that comparison failed the synonymous
control in 31 of 56 genes: real hit children sit shallower than the
size-biased draw expects, whatever the gene. The draw gets each child's
chance of carrying the hit right, but not the extra branch time a real hit
implies, and depth depends on that time. The spectrum statistics are
proportions, which barely depend on depth, and they pass the same control.

**Two more readings to read it against.** The first repeats the test on
internal hits only. Per-genome error cannot reach an internal branch, so a
shift that survives there is not a sequencing artifact, at the cost of the
terminal hits, which are most of them. Each reading starts its random draws
from `--seed`, so it gives exactly what a run of it alone would.

The second asks whether the excess is the list's own. The synonymous
control says the comparison is unbiased. It cannot say whether a
nonsynonymous excess belongs to the genes tested or to nonsynonymous changes
anywhere. Ordinary genes
answer that: `control_genes.py` draws length-matched protein-coding genes,
none of them masked, in the lists given, or with a product matching
`--exclude_product`, and the same test runs on them. Give `--exclude` every
other candidate list as well, such as `replication_genes.txt` below: the
product pattern catches most of their genes but not all (on NC_000962.3 it
misses `dcd`, `nrdH`, `nrdI` and `nrdR`), and a candidate gene does not belong
in a control. Look at `control_genes.txt` before trusting it; that is why
drawing it stays a separate step.

**What to read.** The report reads the three summaries together and prints
the punchline. In order: whether the synonymous control looks like chance; whether the
nonsynonymous hits exceed chance and the control genes; each gene that clears
the bar, with every class at p < 0.05 shown for all hits, internal hits and
its own synonymous hits; and the leads below the bar. It also says when the
smallest p is the permutation floor, and when a control rests on few
effective events. The columns underneath, for reading the file itself: one
row per gene and effect set, plus a pooled `ALL`. `omnibus_p` says whether to look further; `diff_<class>` is the
weighted mean shift in that class's proportion (hit child minus siblings),
with `p_<class>` its permutation p. `effective_events` says how many equally
weighted events a test is worth. The weights grow with subtree size, so one
deep branch can carry most of a test: read a result with a handful of
effective events as a statement about those few branches, whatever its event
count.

**Read the synonymous control first.** It runs by default (`--effects
nonsynonymous,synonymous`). Synonymous hits in the same genes should show no
shift. If they do, the comparison is biased and the nonsynonymous rows cannot
be read either. `--effects lof` restricts to nonsense and start-loss, the
changes most likely to remove function.

`ALL` pools genes with different mechanisms; a gene that raises `CT|GA` and
one that raises `CA|GT` cancel there, so treat it as a screen. With dozens of
genes and six or twelve classes, a handful of per-class p < 0.05 is expected
by chance -- trust consistent results across classes and effect sets over a
single small p.

*Worked example.* `tb_repair_genes.txt` resolves to 56 genes and 42,433
hit events, 25,233 nonsynonymous and 17,200 synonymous, of which 18,185 and
13,857 were informative; about 60% of those are terminal. The control drew
three genes per repair gene, 168 from 3,189 eligible. The run took 29
minutes and 1.3 GB: 8 for every hit, 2 for internal hits, 20 for the control
and under a minute for `dnaE2`'s context test.

The synonymous control is clean with all hits and with internal hits only:
genes with omnibus p < 0.05 number 3 and 2, against 2.8 expected by chance,
and per-class p < 0.05 come to 17 and 12 against 16.8. The nonsynonymous hits
show about twice the chance count, 9 genes and 37 class tests with every hit
counted, and 5 and 26 with internal hits only. Most of that excess points one
way, more `CT|GA` below the hit:

| gene | events | `CT\|GA` shift, all hits | p | internal hits | shift | p |
|---|---|---|---|---|---|---|
| `dnaE2` | 1,109 | +0.93 points | 0.0006 | 218 | +1.00 | 0.006 |
| `udgB` | 376 | +2.11 | 0.0001 | 75 | +1.83 | 0.10 |
| `ogt` | 284 | +1.82 | 0.008 | 57 | +1.47 | 0.12 |
| `lexA` | 254 | +2.80 | 0.04 | 55 | +3.91 | 0.013 |

Only `dnaE2` clears the bar for 56 genes on the omnibus test, at p 0.0001,
the smallest 10,000 permutations can give. Its shift is in more than one
class, and `CT|GA` is not the strongest:

```
  dnaE2: omnibus p 1.0e-04 (the permutation floor; the true p may be smaller)
    class    all hits               internal hits         synonymous
    CG|GC    -0.75 pts, p 1.0e-04   -0.87 pts, p 0.0043   -0.58 pts, p 0.0953
    CT|GA    +0.93 pts, p 6.0e-04   +1.00 pts, p 0.006    +0.25 pts, p 0.444
    TG|AC    -0.29 pts, p 0.041     -0.56 pts, p 0.0101   -0.14 pts, p 0.136
    events   1,109 (54 effective)   218 (14 effective)    771 (2 effective)
```

Below a `dnaE2` hit, C>G and G>C fall as C>T and G>A rise. It is also the
gene that keeps its full shift with terminal hits dropped, on hits where
per-genome error cannot reach. Its 1,109 events are worth 54 effective ones, and no few of them carry
it: dropping the heaviest 1, 5, 10 or 20 leaves the `CT|GA` shift at +1.0
to +1.1 points. Its own synonymous hits are flat for `CT|GA` (+0.25 points,
p 0.44) but lean the same way for `CG|GC` (-0.58, p 0.10), and either way
that control is weak: one deep branch carries 77% of their weight, so they
are worth 2 effective events. The control that carries weight is the one
across genes. `udgB`'s `CT|GA` p of 0.0001 clears the bar too, but as one class of
six, its omnibus p is 0.014, and its internal hits alone do not reach
significance. That is the shape of evidence step 8 could not provide: over a
thousand independent hits, each against its own siblings, so no single
lineage carries it. `udgB` and `ogt`, whose loss would be expected to raise
C>T and G>A, shift the same way by more but on fewer hits, and keep the
direction but lose significance when the terminal hits are dropped; they are
worth a closer look, not a claim. The shift is small throughout: `CT|GA` is half of all
mutations, so +0.93 points is a 2% relative change.

The control genes' nonsynonymous hits behave as chance predicts: 9 of 167
genes (5.4%) at omnibus p < 0.05, against 9 of 56 (16%) for the repair genes,
a difference a one-sided Fisher test puts at p 0.016 (2×10⁻⁴ over the class
tests). No control gene clears 0.05 / 167. So the repair genes' excess is
theirs.

The control's pooled rows show something its per-gene rows are too small
to. A synonymous hit in an ordinary gene is followed by slightly more
`CT|GA`, +0.32 points (p 0.007; omnibus p 0.015, about 3,800 effective
events); a nonsynonymous one by +0.03, which is nothing. A possible reason is
that a hit marks its branch as one whose process produced that kind of
change, and most coding changes are C>T or G>A, while the size-biased draw
matches branch length, not branch spectrum; why it shows for synonymous hits
and not nonsynonymous ones is not clear. An earlier draw of the control,
before the replication genes were excluded, gave +0.32 and +0.14, so read
these as a few tenths of a point, not as exact. `dnaE2`'s +0.93 is well
above the nonsynonymous baseline, and the per-gene false-positive rate in
the control stays at the nominal 5%, but the comparison for a single gene is
with that baseline, not with zero.

**Replication genes.** `tb_replication_genes.txt` lists the replisome and
dNTP-supply genes not already in the repair list, as a separate family: its
17 genes set their own bar rather than raising the repair list's.

```
control_genes.py ... --match tb_replication_genes.txt \
    --exclude resistance_genes.txt,alignment_loci.txt,tb_repair_genes.txt \
    --output control_replication.txt
gene_hit_contrasts.py ... --genes tb_replication_genes.txt \
    --control control_replication.txt --prefix hit_replication \
    > hit_replication.txt
```

A second list is a second run with its own `--prefix` and its own control,
read on its own.

On this tree, in 3 minutes, nothing in it reaches its bar: one gene in each
effect set has omnibus p < 0.05 (`nrdE` nonsynonymous, `dcd` synonymous), as
chance predicts.

**Pathways: gene sets and joint tests.** A single gene often has too few
hits to show anything, and `ALL` mixes genes whose loss should do opposite
things, so they cancel. A pathway sits between the two, and fewer pathways
than genes means a lower multiple-testing bar. `--gene_sets` groups the
genes, and every gene and every pathway is then permuted *together*: each
hit branch gets one random number per permutation, shared by every test
that contains it, and each test maps it through its own size-biased draw.
Each test's own null is unchanged (a one-gene set is permuted with exactly
its gene's draws), but the tests are now coupled, which buys four things:

- **The count of significant genes**, judged against the count expected
  with every gene permuted together rather than a binomial. Genes are not
  independent: one branch can hit two of them, and nested hits share
  subtrees.
- **A family-wise bar**, the p the list's best gene reaches by chance only
  5% of the time (min-p), in place of 0.05 over the number of genes.
- **A pooled test per pathway:** do its genes shift the spectrum the same
  way? It is the gene test run on the pathway's hits together, with a
  branch that hits two member genes counted once and every member's hit
  mutations removed from it.
- **An any-direction test per pathway:** the sum of the member genes'
  omnibus statistics, each scaled by its own null mean so that a gene with
  few events doesn't dominate. Genes pushing opposite ways cancel in the
  pooled test but not here.

On simulated trees, 300 with no effect gave p < 0.05 in 2% (count), 4%
(pooled) and 5% (any-direction). Two genes shifting the same way were found
by the pooled test in 97.5% of 40 trees. Two shifting opposite ways were
found any-direction in all 40 and pooled in 5%: they cancel there, as they
should. Those simulations first ran at 8% to 10% under no effect; the cause
was the simulation, not the test. It drew hit children uniformly where the
test assumes a size-biased draw, and gave two simulated genes clashing node
IDs. Both are fixed, and the context calibration above was remeasured after
the fix.

**The heaviest-branch check.** The weights n1 n2 / (n1 + n2) grow with
subtree size, so a hit on a deep, lineage-defining branch, compared with a
single sibling, can carry a whole test, and that is one lineage against
another. Every gene and pathway row therefore reports the share of its
weight on its two heaviest events, its shift without them, and how much of
the shift that keeps. The report marks a gene among those clearing the bar
or the leads, or a pathway with pooled p below 0.05, as HEAVY when under
half is kept. This automates the check `dnaE2` was given by hand above.
Capping the weights instead, for example by random-effects weighting, was
considered and not adopted: the current weights are efficient when an
effect is real, and the check shows when they are not doing what they
should.

*Worked example.* The repair and replication lists together, plus
`Rv3395c` (`imuA'`) and `Rv3394c` (`imuB`), the rest of `dnaE2`'s mutasome,
and `Rv3202c` (`adnA`) and `Rv3201c` (`adnB`), the AdnAB helicase-nuclease:
77 genes in the 13 pathways of `tb_repair_pathways.tsv`, 10,000 shared
permutations, in 14 minutes and 1.1 GB.

The per-gene results are those of the separate runs, and `dnaE2` is still
the only gene to clear its bar, now 0.05 / 77. Its two heaviest branches
carry 16% of its weight, and without them it keeps all of its shift
(`CT|GA` +1.11 points). The joint count:

| reading | effect set | genes at p < 0.05 | expected, joint null | p |
| --- | --- | --- | --- | --- |
| all hits | synonymous | 5 of 77 | 3.8 | 0.35 |
| all hits | nonsynonymous | 11 of 77 | 3.8 | 0.008 |
| internal hits | synonymous | 4 | 3.8 | 0.50 |
| internal hits | nonsynonymous | 7 | 3.8 | 0.13 |

The excess is real with every hit counted, but weaker than the binomial says
(12 genes, P = 0.0004), and it does not hold on internal hits alone. The
family-wise bar came out at 0.0007, essentially Bonferroni's 0.00065: the
genes behave almost independently, so min-p gains nothing here.

Four leads are HEAVY: `Rv3202c` (`adnA`, 87% of its weight in two branches,
none of its shift kept), `nrdE` (69%, 7% kept), `Rv3395c` (`imuA'`, 59%,
17% kept) and `lexA` (47%, 28% kept). The other seven leads are not.

The pathways, nonsynonymous:

| pathway | pooled p, all hits | largest shift | kept without top 2 | pooled p, internal hits | any-direction p, all / internal |
| --- | --- | --- | --- | --- | --- |
| recombination and SOS | 0.0017 | `CT\|GA` −2.19 pts | 16%, HEAVY | 0.0021 (16%, HEAVY) | 0.081 / 0.071 |
| Y-family translesion (`dinX`, `dinP`) | 0.022 | `CT\|GA` +0.91 pts | 70% | 0.057 (71%) | 0.036 / 0.14 |
| mutasome (`dnaE2`, `imuA'`, `imuB`) | 0.089 | `TC\|AG` +0.67 pts | 55% | 0.11 (37%) | 0.0003 / 0.028 |
| dNTP supply | 0.43 | | | 0.31 | 0.024 / 0.007 |
| single-strand annealing (`recBCD`) | 0.76 | | | 0.80 | 0.37 / 0.46 |

The other pathways have pooled p of 0.068 or more.

Recombination passes the 13-pathway bar (0.0038) on both readings, and it
is a lineage contrast. Two branches carry 46% of its weight, 56% on
internal hits, and nearly all of its shift. Both are `adnA` missense changes
at the tree's deepest binary splits, each compared with a single sibling:
`node_13967` (A748V) is lineages 2, 3 and 4 against lineage 1, and
`node_14005` (E902K) is lineage 2 against lineage 3. Lineage 1 against the
rest is the search's strongest split (step 8), so a hit on that branch
inherits the largest spectrum difference on the tree. The other 6,396
branches add under 0.1 point between them. Each of those two `adnA`
branches is one comparison of two lineages, and it carries every other
mutation that separates them, so the data can't attribute the difference
to `adnA`.

The grouping matters, and both versions belong in any report of this. The
first draft had RecBCD in recombination and no AdnAB, following *E. coli*.
In mycobacteria, AdnAB resects DNA ends for RecA-dependent recombination
and RecBCD works in RecA-independent single-strand annealing (Gupta et al.
2011, *Mol Microbiol* 79:316, in *M. smegmatis*), so the pathways were
regrouped. Under the first grouping, recombination's pooled p was 0.40; the
regrouping brought the two `adnA` branches in.

The Y-family polymerases' `CT|GA` rise is the one pathway lead that doesn't
rest on its heaviest branches, and it runs the same way as `dnaE2`'s. At
0.022 it does not clear the 13-pathway bar. The mutasome's any-direction p
of 0.0003 is mostly `dnaE2`'s own result. `imuA'` shares `dnaE2`'s `CG|GC`
fall, but on two branches (HEAVY above), so that is no support either.

The synonymous rows are not all clean. Replisome's pooled p is 0.045 on all
hits and 0.035 on internal hits, HEAVY on both (33% and 44% kept), and nucleotide
excision's any-direction p is 0.021 and 0.040. Across 13 pathways, two
statistics and two readings, a few nominal results are expected by chance,
but a nonsynonymous result for those two pathways should be read against
its own synonymous rows.

So on this tree the pathway view adds little to `dnaE2`. Its strongest
signal was a lineage contrast that the heaviest-branch check catches, and
its one robust lead points the same way as `dnaE2` but is only nominal.

**A figure.** `plot_gene_hits.py` shows one gene against the rest of its list
and against its own null. It needs the permutation draws, which
`gene_hit_contrasts.py` writes for the genes named in `--null_draws_genes`,
for both readings. The summaries are unchanged by asking for them. It needs
only matplotlib, and `--prefix` finds the run's files:

```
plot_gene_hits.py --prefix hit --gene dnaE2 --class 'CG|GC,CT|GA' \
    --label udgB,ogt,lexA --output dnaE2_figure.pdf
```

Each class gets a column. On top, every gene in the list, both effect sets,
by its shift in the class and its omnibus p, with the list's bar drawn in;
below, the gene's null histogram for each run, with the observed shift, its
p for the class and the effective event count.

**Where in the class: trinucleotide context.** A class shift can be spread
over the class's 16 contexts in proportion to how common each one is, or it
can sit in a few. Which one it is says more about mechanism than the class
total does. `--context_genes` repeats the test for the genes named in the
context classes, in both readings, with the same events, comparators and
controls: the 96 folded ones with `--fold`, the 192 forward-strand ones
without. It needs numpy. Give `--mask` the BEDs the tree was masked with, as
for `profile_tree.py`, and the reference must be the root's genome (step 5).
Name the genes worth a closer look, not the whole list: each costs about 25
seconds and 1 GB per reading. The report ends with a TRINUCLEOTIDE CONTEXT
section, the shape p values per gene, and the files are
`hit_summary_context.tsv` (and `_internal`), with
`hit_context.tsv` (and `_internal`) per context. The context summary
has every column of a `--fold` run, with the class shifts computed by summing
contexts (the self-test checks they equal the 6-class statistic for the same
events). Contexts that cannot be called are dropped, and the draws differ,
so the p values will not match the 6-class run digit for digit. It adds
`shape_p_<class>`: does the class's mix of contexts change below the hit? It
compares each context's share *of its class*, hit child against siblings, so
it asks where the class shift falls, not whether there is one. The test does
not penalise a shift spread in proportion across contexts: in the self-test, a
3× rise of a whole class gave shape p < 0.05 in none of 30 trees while its
class test fired in all of them. A version built on the overall shifts got
this wrong 40% of the time, because a class that really rises has noisier
contexts than any no-shift permutation. `hit_context.tsv` then says which
contexts: per test and context, the context's tree-wide share of its class,
its shift, the shift a proportional spread would give, the excess, the
within-class share shift and their p values. Read the shape p before any one
context's: a strong class shift pushes per-context within-class p < 0.05
above 5%, and with 96 of them a few small ones are expected regardless.

`--calibrate_context N` measures this on N simulated trees per scenario
instead of testing your tree. It needs numpy but not `bte`, and checks the
method, not your data:

```
gene_hit_contrasts.py --calibrate_context 200 --fold --permutations 200
```

| Scenario, 96 classes | shape, every class | shape, C>T | within-class, every context | C>T class test |
| --- | --- | --- | --- | --- |
| no shift | 4.5% | 3.5% | 5.0% | 6.5% |
| C>T 3× in proportion | 6.1% | 7.5% | 6.8% | 100% |
| one C>T context raised | 21% | 100% | 6.8% | 96% |

Each cell is the share of p < 0.05, 200 trees and 200 permutations each, in
40 seconds. The shape test holds near 5% with no shift and under a
proportional shift, and finds a one-context shift every time. Without
`--fold` the 192 classes behave the same: 5.5% shape p with no shift, 5.0% for
C>T's shape under the proportional shift, and 98% power for one context.
(These figures are from simulated trees that draw each hit child with the
size-biased weight the test assumes. An earlier version drew it uniformly,
and gave 5.4% shape p with no shift and 3.5% for C>T's shape under the
proportional shift.)

*Worked example*, `dnaE2`, 10,000 permutations, about 25 seconds and 1 GB per run. The
class result reproduces: `CG|GC` p 0.0001 and `CT|GA` p 0.0003 (internal
hits: 0.005 and 0.007), on the same 1,109 informative events worth 54
effective ones. The context breakdown does not localise the shift:

| | shape p, `CG\|GC` | shape p, `CT\|GA` | shape p, all six |
|---|---|---|---|
| nonsynonymous, all hits | 0.23 | 0.046 | 0.042 |
| nonsynonymous, internal hits | 0.49 | 0.16 | 0.081 |
| synonymous (2 effective events) | 0.62 | 0.096 | 0.47 |

The `CG|GC` fall is spread across its contexts as their frequencies predict.
`CT|GA`'s marginal shape p is not stable: the contexts that lead it with all
hits (`T[C>T]G`, `T[C>T]C`, `A[C>T]A` up; `C[C>T]A` down) are not the ones that
lead it on internal hits (`C[C>T]G` and `A[C>T]A` up; `G[C>T]C` down), and
nothing survives counting 16 contexts. So on this tree `dnaE2`'s shift is a
change in how much of each class there is, not in which contexts it falls.
That is worth knowing: a mechanism that predicts a strong context preference
-- most do -- would have to explain why none shows. It is also as far as 54
effective events go: a +0.93-point class shift split 16 ways leaves about
0.06 points per context, near the noise.

**What it cannot see.** The tree carries SNVs only, so frameshifts,
deletions and IS6110 insertions -- the commonest ways a gene is actually lost
-- are invisible, as are promoter changes. A negative result here does not
clear a gene.

If a late change forces a rerun from the mask on, as a re-mask prompted by a
locus found under a split would, match the new splits to the old ones by
membership: see [Testing the method](#testing-the-method).

## Step 13. Ask whether the splits' spectra show in their resistance mutations

How-to: [Step 13. Ask whether the splits' spectra show in their resistance mutations](HOWTO.md#step-13-ask-whether-the-splits-spectra-show-in-their-resistance-mutations).

The search never sees resistance mutations: step 3 masks their genes because
selection makes the same sites recur on thousands of branches. That
recurrence is also what makes the reverse question answerable. Each
acquisition of, say, `katG` S315T is an independent event on a known branch,
and the step-13 question is whether a split's spectrum shapes those events.
There are two ways it could:

- in **rate**, by making resistance mutations more or less available;
- in **composition**, by tilting which of the resistant variants arise.

### What the tree can show

`resistance_spectra.py` takes the graded variants of the chosen confidence
groups from the catalogue, and keeps those one single-nucleotide
substitution can make:

- missense, stop-gained, start-lost and stop-lost codon changes;
- upstream (`c.-N`) and rRNA (`n.N`) changes.

Two of the catalogue's rules reach past the variants it lists:

- a graded `GENE_LoF` row covers every stop-gained or start-lost
  substitution in the gene;
- a graded `RRDR` row covers every missense substitution in `rpoB` codons
  426 to 452.

A variant the catalogue grades explicitly keeps its own grade, even where a
rule would say otherwise. Frameshifts, in-frame indels and gene deletions
are invisible, which matters most for the drugs whose resistance comes
mostly from losing a gene: ethionamide (`ethA`), pyrazinamide (`pncA`),
streptomycin (`gid`).

Each substitution in a catalogue gene, or upstream of it, is named against
its parent's genome, reconstructed down the tree, as the catalogue would
name it (`rpoB_p.Ser450Leu`). A lineage whose codon differs from H37Rv's
then gets the name its isolates would get.

Naming against the parent needs one guard: a substitution that gives back
the reference's amino acid or base is never an event. Without it, `rpoB`
Leu450Ser on an S450L background, a back-mutation to susceptibility, is a
missense change in the RRDR and the rule counts it as rifampicin
resistance. Such back-mutations are counted and reported per drug, because
a lot of them would mean the tree places resistance mutations badly.

Splits are carried over to this tree by the most recent common ancestor of
their samples, since node IDs don't transfer. Branches are assigned to
residual groups, as in `spectra.tsv`.

### Rate

Each group's spectrum is its trinucleotide-context counts at unmasked
sites, shrunk towards the whole tree's by `--prior` pseudo-mutations. A
context's rate is its share of the spectrum over its opportunity, as in
step 11. A resistance substitution gets its own context's rate. A drug's
events are then shared out over the groups in two ways:

- `E_global`, by each group's mutations;
- `E_spectrum`, by mutations times the summed rate of the drug's
  resistance substitutions under that group's spectrum.

`spectrum_factor` = `E_spectrum` / `E_global` is how much the spectrum
alone changes the supply of resistance mutations. `O / E_spectrum` is what
is left for treatment and selection. The rate side has no p value: how
often, and with what, a lineage's patients were treated swamps mutation
supply, and nothing in the tree models that.

### Composition, and why its null permutes the spectrum

Among a split's events for one drug, two statistics compare its variants
with the rest of the tree's:

- `G` is the 2 x V likelihood-ratio statistic: the variants differ, for any
  reason.
- `llr` asks whether the split's spectrum predicts how they differ.
  - The rest's variant shares `q` are reweighted by `rho`, the ratio of
    each variant's rate under the split's spectrum to its rate under the
    rest's, and renormalised to `q'`.
  - `llr` is the sum over the split's events of `log(q'/q)`. It is positive
    when the events lean towards the variants the split's spectrum favours.
  - `llr_expected` = n KL(q' || q) is the llr the events would give if the
    spectrum were all that differed. Near zero, the spectrum can't move the
    composition enough to see.

`G` is tested by shuffling the drug's events between splits. The obvious
null for `llr` is the same shuffle, and it is wrong. On this tree every
large split's composition differs from the rest's (G q of about 0.01 for
nearly all of them), for reasons unrelated to mutation, such as treatment
history and local epidemics. Under a label shuffle, a composition that
differs that much and happens to lean the spectrum's way looks like a
spectrum effect.

The self-test checks this. It simulates splits whose compositions are
shifted at random, independently of their spectra. Label shuffling called
the spectrum responsible in 27 to 34 cases out of 100. Instead,
`resistance_spectra.py`:

- holds each split's counts as observed;
- permutes `rho` across the drug's variants, with one permutation shared by
  all splits.

The question becomes: does the composition lean towards the variants this
spectrum favours more than towards a random reassignment of the same
ratios? In the same simulation that called a spectrum effect in 2 to 4
cases out of 100. The per-drug total of llr is tested on the same
permutations.

A caution: llr is dominated by a drug's commonest variants. A split
whose events are mostly one variant (`katG` S315T, `rpsL` K43R) passes or
fails largely on where that variant's ratio falls among the others. A drug
with only a handful of variants (amikacin has 4) has too few distinct
permutations for its p values to get small.

### On the MTBC tree

The run used WHO catalogue 3rd edition groups 1 and 2, the round-4 masks,
`--fold`, `--max_branch_length 250` and 10,000 permutations.

**What is visible.** One substitution can make only part of what the
catalogue grades:

| Drug | SNV-type rows | Graded rows |
| --- | --- | --- |
| Ethionamide | 91 | 429 |
| Streptomycin | 58 | 231 |
| Pyrazinamide | 254 | 451 |
| Rifampicin | 110 | 149, plus the RRDR rule |

Frameshifts account for most of the rest.

**Events.** The tree has 44,919 resistance substitutions, 57,075 counting
each drug they confer. Between 50% and 70% are on terminal branches. There
were 574 back-mutations of rifampicin substitutions, against 10,309
rifampicin events, and 445 of isoniazid ones, against 9,718 events: about
5%.

**Rate.** Across the splits, `spectrum_factor` for all drugs pooled runs
from 0.93 to 1.11, so no split's spectrum changes its supply of resistance
mutations by more than about 10%. The observed rates vary far more: among
splits with at least 100 events, from 0.3 times `E_spectrum` (lineage 1
splits `node_4965`, `node_4975`) to 9.3 times (`node_21970`, lineage
2.2.1). The lineage 2.2.1 splits are the richest in resistance, with
`O/E_sp` of 2 to 9, and their spectra predict slightly fewer resistance
mutations, not more (factor about 0.95). Their resistance comes from
treatment and transmission, not from mutation.

**Composition.** `G` is significant for nearly every split with a few
hundred events. The top variants follow lineage: `rpsL` K43R in lineage
2.2.1, `katG` S315T in lineage 3 and lineage 4 splits, and `inhA` c.-777C>T
(`fabG1` c.-15C>T) in several lineage 4 splits.

`llr`, summed over the splits for all drugs pooled, is positive:

| Reading | llr_total | p |
| --- | --- | --- |
| All branches | 320 | 3e-4 |
| Internal branches | 129 | 1e-4 |

So the splits' resistance variants lean, on the whole, towards what their
spectra favour. Few splits carry it on their own at q < 0.05:

- **All branches:** lineage 4.5 (`node_69238`), lineage 3 (`node_36254`),
  and lineage 2.2.1 (`node_14336`, `node_21897`).
- **Internal branches:** `node_36254`, `node_69238`, `node_14331` (lineage
  2.2.1), `node_285` (*M. africanum* La1) and `node_34074` (lineage
  2.2.1).
- **Both readings:** only `node_36254` and `node_69238`.

The clearest single drug is streptomycin:

| Reading | llr_total | p |
| --- | --- | --- |
| All branches | 210 | 1e-4 |
| Internal branches | 82 | 4e-4 |

On internal branches nine splits pass at q < 0.05. Lineage 2.2.1's
`node_14336` leads, with llr 45.7 against an expected 8.7, through `rpsL`
K43R.

An llr five times what the spectrum alone could produce says the
composition leans the spectrum's way, but more strongly than mutation
supply explains. A fitness cost of K43R that differs by lineage would do
the same, and this test can't tell the two apart.

Some drugs go the other way:

| Drug | llr_total, all branches |
| --- | --- |
| Kanamycin | -101 |
| Amikacin | -27 |
| Levofloxacin | -81 |

Their variants lean against the spectra (kanamycin's through the `eis`
promoter variants of lineage 2.2.1 splits), though the permutations do
not call any of these significant.

The conclusion: the splits' spectra do not make resistance more or less
available, and they tilt which resistance variants arise only slightly,
detectably when pooled and in a few lineages, most clearly for
streptomycin. Lineage and treatment history, which `G` and `O/E_spectrum`
show, matter far more.

## Record: the four rounds of tree modification

The tree the worked examples come from was not built in one pass. It went
through a divergence prune, a laboratory-strain prune and four rounds of
masking, and several were run in a different order from the one
[HOWTO.md](HOWTO.md) now gives. This section records what each change was,
what prompted it and what it bought, because the worked examples compare
against the earlier versions, and because a new tree will need rounds of
its own.

### The prunes

The divergence prune at `--threshold 250` was step 2 from the start. The
laboratory-strain prune, of the 336 BCG genomes, came later. The first
masking rounds and searches were run with BCG still in the tree, and BCG
was pruned from the round 3 tree before it was searched again. The run in
the order [HOWTO.md](HOWTO.md) gives makes both prunes at step 2, before any
masking. Pruning and masking commute (below), so the order changes nothing
about the mask. [Testing the method](#testing-the-method) compares the
searches before and after the BCG prune.

### Masking comes in rounds

[HOWTO step 3](HOWTO.md#step-3-build-and-apply-a-region-mask)'s commands build the finished mask, but nobody could have written it
in one go. Masking is iterative: build a mask, apply it, re-profile (step 4),
see which sites still recur most, and add what they point to. Each
build-apply-re-profile cycle is a *round*. The MTBC tree took four, the
last prompted not by the profile but by step 7's null:

| round | adds | mask | prompted by |
| --- | --- | --- | --- |
| 1 | `tb_resistance_genes.txt` with a 200 bp upstream flank, and the `^(PE\|PPE\|esx)` name pattern | 348 kb, 7.9% | Prior knowledge: loci under selection, and the best-known paralogous families |
| 2 | The first section of `tb_alignment_loci.txt` (mobile elements, prophages, repeats), the transposase / integrase / insertion element product pattern, and `--mask_pad 100` | 469 kb | Round 1's surviving top sites, many of them just outside a feature |
| 3 | The round 3 section of `tb_alignment_loci.txt`, `^maturase$`, and `--mask_pad 150` | 504 kb, 11.4% | Round 2's surviving top sites: synonymous-only hotspots, clusters such as `Rv2019`'s just past the 100 bp pad, the `Rv0071` maturase site, and `Rv2082` under a clade-defining branch |
| 4 | Two sites in `lppN`, 2544936 and 2544938 (added by hand, as HOWTO step 3 shows) | 2 bp more | Step 7's hotspot check: the two sites, flipping back and forth inside one 1,593-tip clade, set the `sites` null's maximum |

HOWTO step 3's commands are round 3, and its hand-added sites are round
4. Following [HOWTO.md](HOWTO.md) you can run them directly; you don't
need to repeat the history. The earlier rounds are
here because the worked examples compare against them, and because on a
different tree you will need rounds of your own: you can't write an artifact
list until you have seen which sites recur. Step 4's "Four rounds, and when
to stop" has what each round bought and how to tell when you are done. The
rounds were run on a tree with BCG still in it; pruning it in step 2 instead
changes nothing here, because pruning and masking commute (step 4).

*Worked example.* In round 1, the resistance list resolved to 41 intervals
over 81 kb and the PE/PPE pattern to 151 intervals over 266 kb, 348 kb of the
4,411,532 bp reference together. Round 3's commands, as HOWTO step 3 gives them, cover
506,763 bp, 11.5%, and round 4 makes it 506,765. That is 3,017 bp more than the original round 3, because
`annotate_sites.py` then dropped features typed `pseudogene`. The
`^(PE|PPE|esx)` pattern reaches five of NC_000962.3's: `PE21`, `PE_PGRS36`,
`PPE47`, `PPE48` and `esxM`, which `tb_alignment_loci.txt` also names.

The walkthrough's mask is 1,220 bp larger than the current lists give
(505,545 bp with round 4). `tb_resistance_genes.txt` then listed `ribD`,
which in NC_000962.3 is both Rv2671, the para-aminosalicylic acid gene, and a
synonym of `ribG` (Rv1409), so the mask covered Rv1409 and its 200 bp flank
too. `check_gene_names.py` caught it, and the list now names Rv2671. The
extra 0.2% only removes ordinary sites and adds no artifact, so the
walkthrough was not re-masked.

Bases removed is not the number that matters, though; the next step has it.

### What each round bought

Step 3 lists what each round added. On the MTBC tree the returns fell off
sharply:

| | raw | round 1 | round 2 | round 3 | round 4 |
| --- | --- | --- | --- | --- | --- |
| usable SNVs | 3,249,554 | 2,972,123 | 2,886,050 | 2,845,156 | 2,843,878 |
| worst single site | 6,102 | 1,989 | 318 | 225 | 225 |
| top 20 sites' share | 1.24% | 0.50% | 0.16% | 0.12% | 0.12% |
| `CT\|GA` | 0.5011 | 0.5058 | 0.5067 | 0.5080 | 0.5080 |
| `TC\|AG` | 0.1605 | 0.1581 | 0.1578 | 0.1575 | 0.1575 |

Rounds 1 to 3 are the figures from the original runs, with BCG still in the
tree; round 4 is the walkthrough tree, the same as [step 4's table](#step-4-re-profile-and-check-what-the-mask-bought). The
difference between the two trees' round 3 and round 4 is accounted for
exactly:

- That run used an unanchored `maturase`, which also masked `Rv0693` (see
  below): 1,014 mutations. With the anchored pattern it would have counted
  2,846,170.
- Pruning BCG took 535 mutations, giving 2,845,635.
- The pseudogenes added to the mask in step 3 held 1,520, giving
  2,844,115.
- Round 4's two `lppN` sites held 226. The tree also lost 56 internal
  nodes, and 11 mutations at 11 other sites, one each. Most likely those
  nodes existed only to carry `lppN` changes, and collapsing them merged
  identical mutations on branches that had been split by them. That gives
  2,843,878.

The BCG prune and the pseudogenes were measured in opposite orders. The earlier tree was masked
and then pruned of BCG; this one was pruned and then masked. Counting the
pseudogene positions' mutations on the earlier tree and subtracting gives
this tree's total to the mutation, so pruning and masking commute here, as
they should.

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

### The pad: from 100 to 150, and the maturase

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

### Round 4: the `lppN` sites

*Worked example.* On the walkthrough tree as steps 2 and 3 first built it,
with 20 replicates, the `sites` null's maximum was 45.4, and the summary recommended
`--min_chi 91`. The check flagged a lineage 2.2.1 clade of 1,593 tips,
`node_21012`. With its child `node_21013`, it held the two largest maxima,
45.4 and 31.6, and no other clade's exceeded 24.4. Two sites 2 bp apart in
`lppN` (`Rv2270`), 2544936 and 2544938, carry 224 of its 12,828 mutations,
and 224 of the 226 those sites have in the whole tree:

| mutation | branches | | mutation | branches |
| --- | --- | --- | --- | --- |
| `A2544938G` | 72 | | `T2544936G` | 68 |
| `G2544938A` | 43 | | `G2544936T` | 41 |

Each site keeps flipping back, and 94 branches carry both. Two adjacent sites
changing together, forward and back, a hundred times inside one clade look
like a genotyping or alignment artifact, not independent mutation. On the
earlier tree with BCG still in it the same clade, as `node_21060` and
`node_21059`, held the three largest maxima, 32.0, 29.2 and 26.6. So it also
set that run's `--min_chi 65`, and the tier B split `node_21060` found there
lost 18% of its chi-squared without site 2544938 alone. The check then
counted only how often a node was the maximum, and as two IDs with 2
replicates each the clade passed. The new tree's 45.4 is the same artifact
drawing a rarer class.

So mask the two sites. Add them to step 3's build (its "Sites added by hand"
lines), and mask them on the current tree:

```
printf 'N2544936N\nN2544938N\n' > lppN_mask.txt
matUtils mask -i tree.pruned.masked.pb.gz -m lppN_mask.txt \
    -o tree.pruned.masked.lppN.pb.gz
mv tree.pruned.masked.lppN.pb.gz tree.pruned.masked.pb.gz
```

Then rerun step 4's profile, so that later steps read the final tree, and
`permutation_null.py`. On the MTBC tree the check then comes back clean, over
100 replicates: the two largest maxima come from different clades, and no
clade is the maximum in more than 5 of them.

## Testing the method

[HOWTO.md](HOWTO.md) is what to run on a new tree. This section is what was
run to test the workflow itself, and what to reach for when a result has to
bear more weight than the tiers can carry.

### Replicate the search on a changed tree, and match by membership

The null (step 7) asks whether chance reaches a split's chi-squared, and the
bootstrap (step 8) resamples alignment positions. Neither asks whether the
same splits come back when the set of genomes, and so the search's carving
of the tree, changes. A second search over nearly the same data does, and
because node IDs do not survive the change, the splits are matched by sample
membership.

It is not one of the twelve steps because it is expensive: a search
with bootstrap and step 10, and a new cascade null if the change touched the
mask, hours on a large tree.
Once the method has been tested, the tiers carry most of what it would tell
you: on the larger MTBC replications below, the splits that failed to come
back were tier C or, with one exception, explained by the change itself, and
on a prune alone every split came back but the one the `--min_mutations`
floor removed. It is worth running when a claim rests on a marginal
split, and it is nearly free when a late change has forced a rerun anyway:

- a locus found under a split that has to be masked, which is how this
  tree's round 3 and round 4 masks came about (steps 4 and 7);
- a read check that turns up a systematic misalignment, such as one under
  the `polA` pair on lineage 2's branch (step 8);
- a quality gradient in one clade (step 9), or a lab stock or other
  non-natural population found in a split (step 11);
- a newer release of the tree.

With none of those pending, make a deliberate change. Dropping a random
tenth of the samples tests sensitivity to which genomes were sampled:

```
matUtils summary -i tree.pruned.masked.pb.gz -s samples_masked.txt
tail -n +2 samples_masked.txt | cut -f1 | shuf -n $(( $(wc -l < samples_masked.txt) / 10 )) \
    --random-source=<(yes) > drop.txt
matUtils extract -i tree.pruned.masked.pb.gz -s drop.txt --prune \
    -o tree.pruned.masked.sub.pb.gz
```

Then rerun step 8 on the new tree with the first run's parameters,
`--min_chi` included, since the point is to change only the tree, and
step 10. Step 7 is needed again only if the change touched the mask. A prune
leaves the recurrent sites that drive the sites null as they were, so the
first run's `null_sites_cascade.tsv` still gives the right ceilings; a new
run of it would only be another independent draw from the same
distribution. A re-mask changes those sites and needs a new cascade null.
The example below used a mask round plus a BCG prune, because that is the
change the MTBC tree went through; the `_nobcg` file names are that run's.

Run step 8 in a directory of its own, since its bootstrap files would
overwrite the first run's. Then build the new run's step 10 inputs with the
same loops as before, on the new tree and the new run's files, and score it:

```
for f in bootstrap_*_splits_output.tsv; do n=${f#bootstrap_}; n=${n%%_*};
  cut -f1 "$f" | grep ^node_ | sed "s/^/$n\t/"; done > bootstrap_by_replicate_nobcg.tsv

for n in $(cut -f1 spectra_nobcg.tsv | grep ^node_); do
  matUtils extract -i tree.pruned.masked.nobcg.pb.gz -I $n -u split_samples.txt >& /dev/null
  sed "s/^/$n\t/" split_samples.txt
done > split_clade_samples_nobcg.tsv

# --null: a new cascade null instead, if the change touched the mask
split_confidence.py --splits_log splits_nobcg.log --spectra spectra_nobcg.tsv \
    --null null_sites_cascade.tsv \
    --bootstrap_replicates bootstrap_by_replicate_nobcg.tsv \
    --clade_samples split_clade_samples_nobcg.tsv \
    --output split_confidence_nobcg.tsv
```

The two `node<TAB>sample` tables, the first run's `split_clade_samples.tsv`
and the new one, are what the matching needs. `match_splits.py` pairs each
old split with the new split sharing the most samples, by Jaccard index
(shared samples over samples in either):

```
match_splits.py --old split_clade_samples.tsv \
    --new split_clade_samples_nobcg.tsv --output split_matches.tsv > split_matches.txt
```

Both sides are first cut down to the samples present in both trees, read from
their `node_1` rows, so a sample the change removed does not count against
the clade that held it. It reports how many old splits match exactly, how
many nearly (Jaccard 0.95 or more, `--near`) and how many have no
counterpart (below 0.5, `--lost`), then one row per old split with its best
match, the samples only on each side, and whether the match is mutual. New
splits that are no old split's mutual match are listed last: they are the
change's own findings, not failures to replicate.

*Worked example: a prune alone.* The round 3 tree was searched again with
its 336 BCG genomes pruned, at the same settings, and scored against the
first run's cascade null. On the 128,918 samples present in both trees:

```
  68 old splits, 68 new
  Jaccard 1.000: 67
  no counterpart (best Jaccard < 0.5): 1: node_10
  new splits that are no old split's mutual match: 1: node_286 (best node_6, 0.876)
```

All 67 matches are mutual and exact, to the sample. Every tier C split came
back. The one loss and the one new split are the same event, below.

### The loss was the `--min_mutations` floor, not the biology

`node_10` was tier A at 7.5x margin and 100% support, and it is gone. The
reason is arithmetic. Its residual held 10,077 mutations, and pruning BCG
removed 535 of them, all necessarily inside it since BCG sits nowhere else.
The 577 samples left carry 9,542 mutations, below `--min_mutations 10000`.
The clade did not stop being different; it stopped being *testable*.

Watch what the search did instead, because it is the reason this is worth
documenting rather than filing as a loss. It carved out the **complement**:
a new 7,356-tip split (`node_286`, tier A, 99% support) inside the same
parent, leaving a 768-tip residual that contains all 577 surviving `node_10`
samples and has `node_10`'s spectrum (`CA|GT` 0.0791 and `TC|AG` 0.2635,
against the old 0.0777 and 0.2688; the parent's old residual had 0.2456).
`node_286` shares not one sample with `node_10`.

So the signal is fully intact and has moved from a split to a residual. Two
things follow. A clade dropping off the split list is not evidence it stopped
being distinctive, and you cannot tell the difference by reading node IDs or
split lists -- only by asking where its samples went. And a `--min_mutations`
floor does not merely hide small clades, it can *relocate* a finding from the
split you would report into a residual you might not think to look at.

### Tier changes between runs are mostly noise

Of the 67 matched splits, 61 kept their tier, and the second run's counts
were 35/29/4 against 35/27/6. All six that moved had unchanged margins and
moved on support alone, across the 70% line: `node_73757` from 59% to 70%,
`node_74367` from 72% to 62%. With 100 replicates, support near 70% varies
by about 9 points between runs, so a tier border is soft on that axis as
well as on margin.

The null is no steadier. The pruned tree's own cascade null was also run,
for comparison. Scored against it instead of the first run's, margins
changed by a factor of 0.64 to 1.37 (median 1.08) and three splits changed
tier, two up and one down. That is the scatter of two independent draws of
the same noisy ceiling, which is why reusing the first run's null after a
prune is sound.

### A larger change: round 2 to round 3

The round 2 tree's splits, matched against this same round 3 tree with BCG
pruned, test a re-mask as well as the prune. That gave 59 splits against the
previous 62, and **58 of the 62 matched at Jaccard 1.000** -- 56 of them
exactly, two off by a handful of tips. Tier counts went 36/20/6 to 33/22/4.
That is the result that says the search is finding structure rather than
fitting the artifact landscape, and no amount of null calibration or
bootstrapping could have established it.

What failed to replicate is what the tiering predicted. Three of the four
lost splits -- `node_74277`, `node_73308`, `node_73473` -- were tier C, and
were the three *lowest* margins in the entire old set (2.37x, 2.41x, 2.83x).
The three tier C splits that did survive came back tier C. The fourth loss
was round 2's `node_10`, lost to the floor exactly as above. So under a
change large enough to move the recurrent sites, the tier was a prediction
about reproducibility, and it held. Under a prune alone, nothing but the
floor case failed.

### The walkthrough itself: prune at step 2, and mask `lppN`

The run in the order HOWTO.md gives is a replication too. Against the round 3
run with BCG, its tree differs by the BCG prune, now made at step 2, and by
round 4's two `lppN` sites, and it was scored against its own cascade null,
as a re-mask requires. Matched by membership:

```
match_splits.py --old old/split_clade_samples.tsv \
    --new split_clade_samples.tsv --output split_matches.tsv > split_matches.txt
```

```
SPLIT MATCHING by sample membership, on the 128,918 samples present in both trees
  68 old splits, 66 new
  Jaccard 1.000: 63
  Jaccard >= 0.95 but < 1: 0
  in between: 2
  no counterpart (best Jaccard < 0.5): 3: node_45190 node_21060 node_10
  ...
  new splits that are no old split's mutual match: 1: node_285 (best node_6, 0.876)
```

`node_10` and `node_285` are the prune, exactly as in the prune-alone run
above. The other two losses are worth tracing, because each failed for a
different reason, and neither reason is visible in a list of split IDs.

`node_21060` (lineage 2.2.1, 1,590 samples, tier B at 2.8x) was the split
step 8's `split_drivers.py` found leaning hardest on one site: 18% of its
chi-squared came from 2,544,938, one of the two `lppN` sites. Masking them
took 241 mutations from its residual. Recomputed on the new tree, the same
clade (`node_21013`) scores 29.3 against its parent's row, down from 44.7,
and `--min_chi 30` no longer accepts it. This is the loss the re-mask was
meant to cause: the split rested partly on an artifact.

`node_45190` (lineage 4.1.2.1, 9,881 samples, tier C at 2.45x and 67%
support) lost to a neighbour. Its own residual barely changed: 150,527
mutations against 150,584. But it is tested against the rest of its parent's
row, and that row grew from 66,371 mutations to 71,372. The 5,001 extra came
from `node_51974` (lineage 4.1.1.1), which the new search found one node
lower, as `node_51867`, with 342 fewer samples (Jaccard 0.67). Their 5,033
mutations stayed in the parent's row and diluted the contrast, and the clade
(`node_45082` now) scores 27.8. This is the carving at work (step 10): a
split's chi-squared depends on what else the search has already removed
from the table, so a marginal split can fall below the bar without anything
inside it changing.

The other partial match, `node_69997` to `node_69899` (605 samples to 464,
Jaccard 0.77), was tier C on both runs.

Of the 65 matched splits, 57 kept their tier, and the counts went 35/27/6 to
33/28/5. The eight that moved were all within about 15% of the 5x margin
line or the 70% support line on one run or the other. Margins moved by a
factor of 0.69 to 1.30, median 0.91. That is wider than the prune alone
gave, as expected, since the two runs were scored against different nulls,
and the second null takes the 0.99 quantile of each size window rather than
its maximum (step 7).
