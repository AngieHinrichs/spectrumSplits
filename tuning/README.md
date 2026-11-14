# Parameter tuning helpers

Three read-only scripts for choosing `spectrumSplits.py` parameters on a tree
other than the one the defaults were calibrated on. None of them modifies the
tree, and all three need only `bte` and the Python standard library.

| Script | What it answers |
| --- | --- |
| `profile_tree.py` | What is in this tree, and what do its numbers imply for `--min_mutations`, `--max_branch_length` and whether to fold? |
| `permutation_null.py` | How large a chi-squared does chance alone reach on this tree, and so where should `--min_chi` sit? |
| `annotate_sites.py` | What are the recurrent sites, and which regions should be masked before any of the above is trusted? |
| `branch_class_spectra.py` | Is a split a real mutational difference, or a sequencing-quality gradient on terminal branches? |

The steps depend on one another, so the order matters more than any single
setting: masking changes the spectrum, the spectrum sets the floor under
`--min_mutations`, and both reshape the null that `--min_chi` is calibrated
against. [Suggested order](#suggested-order) at the end is the part to work
from; the sections in between exist to explain why each of its steps is there.

## Why the spectrumSplits.py defaults need calibration

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
  and against `scipy.stats.chi2_contingency`, and exits. Worth running once in
  your environment to confirm the statistic matches `spectrumSplits.py`.

The statistic is reimplemented in pure Python rather than called through scipy,
so all candidate nodes can be scored in one vectorless pass without scipy's
per-call overhead. It matches `chi2_contingency` exactly for tables larger than
2 x 2 (no continuity correction); columns whose total is zero are skipped, where
scipy would raise.

## 3. `annotate_sites.py` — identify recurrent sites, and build mask regions

`profile_tree.py` tells you *which* sites mutate repeatedly; this tells you
*what they are*, so you can decide what to mask. It does two jobs.

**Identify sites.** Reads positions straight out of a profile report:

```
annotate_sites.py --annotation NC_000962.3.gff3.gz --reference NC_000962.3.fna \
    --from_profile mtb_profile.txt --output annotated_sites.tsv
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
annotate_sites.py --annotation NC_000962.3.gff3.gz \
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

### A caution about the gene list

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

### The other gene list: `tb_repair_genes.txt`

`tb_repair_genes.txt` lists DNA repair, replication-fidelity and dNTP-pool
genes, and exists for the opposite purpose. **Do not mask these.** A
repair-gene mutation is the *explanation* you want when a clade turns out to
have a distinct spectrum: a lineage that lost or altered a repair function is
what a real mutational-spectrum split looks like. Masking them destroys the
signal.

Use it with `--highlight`, which flags matching sites in a `highlight` column
and counts them, changing nothing else:

```
annotate_sites.py --annotation NC_000962.3.gff3.gz --reference NC_000962.3.fna \
    --positions <mutations on the clade's defining branch> \
    --highlight tb_repair_genes.txt
```

This list carries gene names and **no locus tags at all** — the Rv numbers are
in comments as hints to check. A name missing from your annotation is reported
as NOT FOUND, which is harmless, whereas a misremembered locus tag would resolve
to the wrong gene and mislead you silently. Verify any you rely on against
Mycobrowser or the NCBI H37Rv annotation.

### The third gene list: `tb_alignment_loci.txt`

`tb_alignment_loci.txt` is also masked, but for the opposite reason to
`tb_resistance_genes.txt`. Resistance loci recur because the same mutation is
selected repeatedly — the variants are real. These loci recur because short
reads cannot be placed unambiguously: paralogous families (PE/PPE, `esx`, the
`Rv1318c`-`Rv1320c` adenylate cyclases, `lldD1`/`lldD2`), repeats (IS6110, the
13E12 family) and prophages (phiRv1, phiRv2). There the variants are largely
not real at all.

Keeping the two lists apart keeps each auditable, and they want different
widening, so resolve them separately and concatenate:

```
annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions tb_resistance_genes.txt --promoter_flank 200 \
    --mask_bed selection_loci.bed

annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions tb_alignment_loci.txt --mask_pattern '^(PE|PPE|esx)' \
    --mask_product_pattern '(?i)(transposase|integrase|insertion element)' \
    --mask_pad 100 --mask_bed alignment_loci.bed
```

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

## Node IDs are not stable across trees

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
    --threshold 25 --output_tree tree.pruned.pb.gz
```

Choose `--threshold` from the MUTATIONS:TIPS section rather than accepting the
auto-threshold.

*Worked example.* The report shows mean 18.6, median 7, p90 54, p99 114,
p99.9 192, max 909, and an auto-threshold (mean + 2 SD) of 71.7 with 10,930
nodes at or above it. Taking the auto-threshold would prune 10,930 nodes — 5%
of the tree — on a species where a high ratio is often real divergence. A
threshold well above p99.9 removed only 60 samples here and still cut the
worst terminal branch from 909 mutations to 234, which is the whole point: the
absurd genomes go, the diverged ones stay.

Read the distribution before choosing. A median of 7 with a p90 of 54 is not a
tight distribution with outliers, it is a broad one, and mean + 2 SD on a broad
distribution cuts into the bulk.

### 3. Build and apply a region mask

```
annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions resistance_genes.txt --promoter_flank 200 \
    --mask_bed selection_loci.bed --mask_positions selection_pos.txt

annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions alignment_loci.txt --mask_pattern '^(PE|PPE|esx)' \
    --mask_product_pattern '(?i)(transposase|integrase|insertion element)' \
    --mask_pad 100 \
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

*Worked example.* The first pass resolved 41 intervals over 81 kb from the
resistance list; the PE/PPE pattern gave 151 intervals over 266 kb — 348 kb in
all, 7.9% of the 4,411,532 bp reference. The alignment list and `--mask_pad`
came later, after the first re-profile showed what had survived, which is the
normal way round: you cannot write that list until you have seen which sites
recur.

Bases removed is not the number that matters, though; the next step has it.

### 4. Re-profile, and check what the mask bought

```
profile_tree.py --input_tree tree.pruned.masked.pb.gz > profile_masked.txt

annotate_sites.py --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --from_profile profile_masked.txt --output surviving_sites.tsv
```

*Worked example.* Judge the mask by recurrence, not by total count:

| | raw | pruned + masked |
| --- | --- | --- |
| usable SNVs | 3,249,554 | 2,972,123 (−8.5%) |
| distinct mutated sites | 1,988,164 | 1,865,580 |
| mutations per site (mean) | 1.634 | 1.593 |
| worst single site | 6,102 | 1,989 |
| top 20 sites' share | 1.24% | 0.50% |

Masking 7.9% of the genome removed 8.5% of the mutations, so by total count it
looks like it bought nothing. It bought a great deal: 13 of the raw top 20
sites are gone, including all seven largest, and the worst hotspot fell by
two thirds. Those are the sites that break the independence assumption.

The reason the total barely moved is worth knowing, because it is easy to
misread. Only 122,584 of the ~348,000 masked positions had ever mutated — 35%,
against 45% genome-wide — because the upstream variant caller had already
excluded most of the PE/PPE families, which were 76% of the mask. So the
masked regions were below average in how many positions varied at all, and far
above average in how often those positions recurred.

If the new top sites in `surviving_sites.tsv` still annotate to resistance
loci, widen `resistance_genes.txt` instead of masking those sites one by one.
Watch in particular for a site where every possible substitution is
synonymous, or where a third codon position recurs hundreds of times:
selection cannot produce that, so it is a mapping artifact, and paralogous
gene families are the usual source.

Node IDs changed when the tree was masked, so compare clades between the two
profiles by sample membership, not by ID.

### 5. Decide whether to fold

*Worked example.* The 12-type spectrum gives C>T 0.2538 against G>A 0.2521,
C>A 0.0604 against G>T 0.0604, A>G 0.0794 against T>C 0.0788, and so on: every
complementary pair agrees to within 1.2%. There is no strand-asymmetric
process to preserve, so fold — the folded classes double every cell count and
take the degrees of freedom from 11 to 5. If one pair disagreed two-fold,
folding would average away a real asymmetry and should be skipped.

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
    --replicates 20 --include_observed --seed 1 --nthreads 8 \
    --output null_sites.tsv
```

Use `--output` for the per-split TSV; the summary goes to stdout, so
redirecting stdout saves the summary and silently discards the table. Pass the
same `--fold`, `--min_mutations` and `--max_branch_length` you settled above —
a threshold calibrated under different settings does not transfer — and the
same `--seed` in both, so the replicates are paired and the only difference
between the runs is the null model. `--include_observed` matters: it scores the
real tree through the identical code path, so the comparison is not against a
number from a different pass.

*Worked example*, the same masked tree as the steps above, folded, at
`--min_mutations 10000` and `--max_branch_length 250`:

| | min | median | max |
|---|---|---|---|
| `types` | 15.2 | 19.7 | 23.9 |
| `sites` | 17.0 | 26.3 | 62.7 |

Observed: **18,283.3**, at `node_5045`, splitting 592,967 mutations below from
2,287,262 above.

Three things to read off that.

**Which null to believe.** `sites` tops out 2.6x higher than `types`, and 11 of
its 20 replicates exceed the *entire* `types` maximum. A threshold taken from
the `types` null would be exceeded by chance more than half the time on this
tree. Calibrate against `sites`.

**Whether to keep masking.** The 38.8-unit gap between the nulls is what
residual recurrent-site clumping costs. That is the number to weigh against
another round of masking — and here it is immaterial, because the observed
maximum is 291x the `sites` maximum. Masking further could move the threshold
from 63 to perhaps 40, which changes nothing. Stop masking when the gap
between the nulls is small compared with the observed value, not when the
recurrence table looks tidy.

**The threshold.** A 2x margin on the `sites` maximum gives `--min_chi 126`,
and the observed maximum is still 146x that. When the headroom is this large
the exact choice does not matter; when it is not, prefer the larger margin.

One caveat on transferring this number: it was calibrated with the default
`--max_splits 1`, so it is the threshold the *first* split must beat. Deeper
splits in the greedy cascade are scored on smaller subtrees, where both the
observed chi-squared and the null are smaller. Re-run with `--max_splits 25`
before trusting one threshold to govern the whole cascade.

#### The threshold should scale with subtree size

That caveat turns out to matter more than a caveat usually does, so do this
run rather than skipping it. Score the cascade null with a permissive
threshold, since here `--min_chi` exists only to stop the search and setting
it high would truncate the null at exactly the depth being measured:

```
permutation_null.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_mutations 10000 --max_branch_length 250 --mode sites \
    --max_splits 25 --min_chi 1 --replicates 20 --include_observed \
    --seed 1 --nthreads 8 --output null_sites_cascade.tsv
```

The printed summary will not answer the question: it reports the maximum
per replicate, which with `--max_splits 25` is still each replicate's first
split, so it comes back close to the `--max_splits 1` run. The answer is in
the TSV. Group the null splits by the size of the subtree each was carved
from (`mutations_below + mutations_above`) and compare each observed split
against nulls from subtrees within, say, 2x its size.

*Worked example*, same tree and parameters as above, 500 null splits:

| parent subtree | null max chi-squared |
|---|---|
| 2.88M (whole tree) | 62.7 |
| 470-590k | 18.0-25.8 |
| 120-180k | 14.0-17.5 |
| 57-92k | 12.4-15.5 |

The ceiling falls by roughly 4x from the root to a 90k subtree, because
chi-squared scales with the counts in the table. A single `--min_chi`
calibrated on the first split is therefore far too strict everywhere below
it: 126 against a real ceiling of about 12 is an order of magnitude.

That cost real splits on this tree. Because `--min_chi 1` scores everything,
the cascade run also reveals what the production run was rejecting — six
splits between 54.7 and 111.1 that each clear their size-matched null by
3.1x to 8.9x, the largest being `node_10` at 111.1 against a ceiling of
12.4. Only one scored split genuinely failed, at 8.7 against 15.5.

Two ways to use this. The simple one is to set `--min_chi` from the
*smallest* subtree you intend to split rather than from the root, which on
this tree means something nearer 40 than 126. The better one is for the
threshold to scale with the candidate subtree's mutation count, which is a
change to `spectrumSplits.py` rather than a tuning choice; until it exists,
a fixed threshold is a known limitation and worth stating alongside results.

Two caveats on the size-matching itself. The 2x windows overlap, so the
ceilings for nearby sizes are not independent estimates. And the null
subtrees are the permuted tree's own, which do not correspond one-to-one
with the observed subtrees being compared against them — the matching is on
size, not identity.

If one node is the maximum in most `sites` replicates, a single hotspot is
setting your threshold: find and mask it rather than paying for it forever.
Here no node repeats more than twice across the 20 replicates, which is what
an absence of such a hotspot looks like.

### 8. Run the search and find where it saturates

```
spectrumSplits.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_chi 126 --min_mutations 10000 --max_branch_length 250 \
    --output_spectrum spectra.tsv 2> splits.log
```

*Worked example (hypothetical — this run has not been done).* Suppose
`splits.log` records each accepted split's chi-squared as 18283, 4100, 1900,
640, 310, 180, 140, 129. The run is still finding splits just above the
threshold, so the curve has not bent yet — either raise `--min_chi` and keep
the splits that clear the higher bar, or keep the rest and label them
marginal. A run that instead stopped at 18283, 4100, 1900 with nothing under
1000 has found its own ceiling, and the threshold is not what is limiting it.

Note that `--min_chi` governs every split, while step 7 calibrated it against
the first one. The deeper splits in a list like that are scored on smaller
subtrees, so a threshold that is strict for the root can be loose by the
fourth split. The `--max_splits 25` null run settles that, and on this tree
the answer was no: see "The threshold should scale with subtree size" above.
A run that stops with its last split a fraction above the threshold has not
saturated — it has hit the bar, and splits below it are going unexamined.

### 9. Confirm, then interpret

```
spectrumSplits.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_chi 126 --min_mutations 10000 --max_branch_length 250 \
    --bootstrap_splits 100 --nthreads 8 --output_spectrum spectra.tsv
```

Then read each split's per-class proportions out of `spectra.tsv`. A split
driven almost entirely by one mutation class is usually an artifact or a single
over-mutated site; one spread across several classes that survives the
bootstrap is the thing worth explaining. Running `annotate_sites.py --highlight
repair_genes.txt` on the mutations of the split's defining branch is the first
place to look for a mechanism.

**Check each split's lineage before reading a mechanism into it.** If the tree
carries clade annotations, this is a few commands and it changes the
interpretation completely:

```
matUtils summary -i tree.pruned.masked.pb.gz --sample-clades clades.txt
matUtils extract -i tree.pruned.masked.pb.gz -I node_N -u node_N_samples.txt
grep -Fwf node_N_samples.txt clades.txt | cut -f2 | sort | uniq -c | sort -nr
```

*Worked example.* The three largest splits on the tree above each turned out
to be a recognised, deeply diverged group, and the tip counts matched the
split sizes exactly:

| split | clades | tips | share |
|---|---|---|---|
| `node_5045` | lineage 1 and sublineages | 11,468 | 8.9% |
| `node_14473` | lineage 2.2.1 (Beijing) | 35,044 | 27.1% |
| `node_6` | La1/La2/La3, the animal-adapted members | 8,729 | 6.8% |

Two things follow, in opposite directions.

The repair-gene reading collapses for `node_14473`. One of only five
mutations on its defining branch is `mutT4` R48G, a 12x enrichment — and it
is lineage 2.2.1, so that is the documented Beijing `mutT4` polymorphism. The
variant and the spectrum are both properties of lineage 2, and a single
branch cannot separate them. Any deep branch carries hundreds of
lineage-defining mutations, so co-occurrence is guaranteed. What would be
evidence is independent recurrence: the same gene disrupted on an unrelated
branch with the same shift.

The quality reading collapses too, more thoroughly than step 10 could show.
The clade that broke the divergence correlation — high mutations per tip, but
a spectrum matching its low-divergence parent — is `node_3616`, with 269
tips. La3 also has exactly 269 samples, so that clade *is* La3, a distinct
animal-adapted species. Its long branches are sparse sampling of a diverged
species, and its spectrum matches its parent because it belongs to the
animal-adapted group. The anomaly has a phylogenetic explanation, which is
the strongest single sign that these spectra track lineage.

One caveat on `node_6`: 336 of its 8,729 tips are `La1.2.BCG`, a
laboratory-passaged vaccine strain. That is 3.9%, so it will not dominate,
but serial passage is not natural evolution and it is worth excluding before
making a quantitative claim about the animal-adapted spectrum.

Extract the split's samples before doing anything else to the tree. Node IDs
are assigned at read time, so the first mask, prune or collapse renumbers
every node and the ID in `splits.log` stops meaning anything (see "Node IDs
are not stable across trees" above). The sample set is the only durable
handle:

```
matUtils extract -i tree.pruned.masked.pb.gz -I node_5045 -u split_samples.txt
```

`-I` selects by internal node ID, which is what `splits.log` reports. `-c` is
for annotated clades — the lineage labels stored as node annotations in a MAT,
such as Pango names — and will not find a `node_N` ID.

### 10. Check the splits against sequencing quality

A split can be real biology or it can be a quality gradient, and the two look
identical in `spectra.tsv`. On a sparsely sampled tree most mutations sit on
terminal branches — 63% on the tree described here — and a terminal branch
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
    --output branch_class.tsv
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

*Worked example*, the 30 splits from the run above. The correlations with
divergence per tip are the same on both branch classes:

| class | internal r | terminal r | internal rho | terminal rho |
|---|---|---|---|---|
| `CA\|GT` | 0.429 | 0.461 | 0.450 | 0.347 |
| `CG\|GC` | -0.327 | -0.378 | -0.253 | -0.274 |
| `CT\|GA` | **0.644** | **0.643** | 0.646 | 0.636 |
| `TA\|AT` | -0.040 | -0.132 | -0.252 | -0.289 |
| `TC\|AG` | **-0.523** | **-0.524** | -0.644 | -0.648 |
| `TG\|AC` | **-0.691** | **-0.712** | -0.674 | -0.709 |

Within each clade, internal and terminal spectra differ by a median of 0.1 to
0.4 percentage points per class, with a single maximum of 2.28. So the
spectrum is a property of the clade and not of the branch class: terminal
branches are not carrying a distinct artifact spectrum, even though they hold
63% of the mutations. The quality explanation is ruled out for this tree.

The between-clade separation therefore stands on internal branches alone,
where per-genome error cannot reach: `TC|AG` runs 0.131 to 0.273 across
clades, a 2.08x spread on counts of 4,000 to 63,000 internal mutations each.
And the clade that broke the divergence correlation holds up — at a
divergence ratio of 46.3, close to the 54.3 of the most extreme clade, its
internal-only `TC|AG` of 0.230 matches its low-divergence parent's 0.246
rather than that clade's 0.142.

Note what this step does *not* establish. It rules out per-genome artifacts,
not shared ones: a repeat that misaligns the same way in every genome of a
lineage places its mutations on internal branches too. Masking is what
addresses those, which is why it comes first. A clean result here means the
splits are not a sequencing-quality gradient; it does not by itself make them
biology.
