# Parameter tuning helpers

Three read-only scripts for choosing `spectrumSplits.py` parameters on a tree
other than the one the defaults were calibrated on. None of them modifies the
tree, and all three need only `bte` and the Python standard library.

| Script | What it answers |
| --- | --- |
| `profile_tree.py` | What is in this tree, and what do its numbers imply for `--min_mutations`, `--max_branch_length` and whether to fold? |
| `permutation_null.py` | How large a chi-squared does chance alone reach on this tree, and so where should `--min_chi` sit? |
| `annotate_sites.py` | What are the recurrent sites, and which regions should be masked before any of the above is trusted? |

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

Names that don't resolve against your annotation are reported as NOT FOUND
rather than silently dropped.

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
annotation and sequence, `resistance_genes.txt` a list of loci to mask. The
numbers in the worked examples are illustrative, not recommendations.

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

*Worked example.* The report shows mean 4.4, p99 15, p99.9 22, max 310, and an
auto-threshold (mean + 2 SD) of 11.1 with 2,400 nodes at or above it. Pruning
2,400 nodes for exceeding 11 would discard a lot of real divergence on a
sparsely sampled tree. A threshold of 25 sits above p99.9, catches the max-310
outlier and its neighbours, and leaves the rest alone.

### 3. Build and apply a region mask

```
annotate_sites.py --annotation NC_012345.1.gff3.gz \
    --regions resistance_genes.txt --promoter_flank 200 \
    --mask_pattern '^(PE|PPE)' \
    --mask_bed mask.bed --mask_matutils matUtils_mask.txt
```

`--mask_matutils` writes the mask in the form `matUtils mask -m` expects — one
line per masked position, `N<position>N` — so it needs no conversion:

```
matUtils mask -i tree.pruned.pb.gz -m matUtils_mask.txt \
    -o tree.pruned.masked.pb.gz
```

Expect this to take a while: `matUtils mask` traverses the tree once per masked
position.

*Worked example.* The resistance loci resolve to 41 intervals over 81 kb and
the PE/PPE families to 151 intervals over 266 kb — 348 kb in all, about 8% of a
4.4 Mb genome. Bases removed is not the number that matters, though; the next
step has it.

### 4. Re-profile, and check what the mask bought

```
profile_tree.py --input_tree tree.pruned.masked.pb.gz > profile_masked.txt

annotate_sites.py --annotation NC_012345.1.gff3.gz --reference NC_012345.1.fna \
    --from_profile profile_masked.txt --output surviving_sites.tsv
```

*Worked example.* Usable SNVs fall from 4.1M to 3.3M — 20% of the mutations
gone for 8% of the genome, a 2.5x enrichment that says the mask was well aimed.
Had it taken 8% of the mutations as well, the mask was spectrum-neutral and the
recurrence problem lives somewhere you have not looked. If the new top sites in
`surviving_sites.tsv` still annotate to resistance loci, widen
`resistance_genes.txt` instead of masking those sites one by one.

Node IDs changed when the tree was masked, so compare clades between the two
profiles by sample membership, not by ID.

### 5. Decide whether to fold

*Worked example.* The 12-type spectrum gives C>T 0.184 against G>A 0.202, C>A
0.055 against G>T 0.064, and so on down the six pairs: every complementary pair
agrees to within about 4%. There is no strand-asymmetric process to preserve,
so fold — the folded classes double every cell count and take the degrees of
freedom from 11 to 5. If one pair disagreed two-fold, folding would average
away a real asymmetry and should be skipped.

Settle this now: it changes both the floor in step 6 and the null in step 7.

### 6. Pick `--min_mutations` and `--max_branch_length`

*`--min_mutations`.* The masked profile reports a rarest 12-type class at
p = 0.0085 (floor 590) and a rarest folded class at p = 0.0175 (floor 286).
Folding, so the binding floor is 286. SPLIT CANDIDATE AVAILABILITY shows 2,500
leaving 140 eligible nodes with `smallest_eligible_subtree_tips` 410, and 5,000
leaving 60 nodes. Take 2,500: it clears the floor nearly nine times over, keeps
enough candidates for a real search, and no 410-tip clade can be a single bad
genome.

*`--max_branch_length`.* Branch counts run median 3, p99 16, p99.9 29, max 950,
and the cap table shows a cap of 50 excluding 3 branches and 0.4% of mutations
while a cap of 10 excludes 66 branches and 37%. Cap at 50: it removes the few
genuinely implausible branches at almost no cost. A cap of 10 would discard a
third of the data, which on a diverged tree is real divergence rather than
noise.

### 7. Calibrate `--min_chi`

```
permutation_null.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_mutations 2500 --mode types --replicates 20 \
    --include_observed --nthreads 8 > null_types.tsv

permutation_null.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_mutations 2500 --mode sites --replicates 20 \
    --include_observed --nthreads 8 > null_sites.tsv
```

Pass the same `--fold`, `--min_mutations` and `--max_branch_length` you settled
above; a threshold calibrated under different settings does not transfer.

*Worked example.* `types` tops out at 48, `sites` at 394, and the observed tree
scores 20,568. The 8x gap between the two nulls is the residual recurrent-site
clumping. Set `--min_chi` above the worse null with margin, so 800 — and note
the observed maximum is 52x that, which is not a close call. If one node were
the maximum in most `sites` replicates, a single hotspot would be setting your
threshold: find and mask it rather than paying for it forever.

### 8. Run the search and find where it saturates

```
spectrumSplits.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_chi 800 --min_mutations 2500 \
    --output_spectrum spectra.tsv 2> splits.log
```

*Worked example.* `splits.log` records each accepted split's chi-squared:
20568, 9100, 4400, 2050, 1180, 940, 860, 815. The run is still finding splits
just above the threshold, so the curve has not bent yet — either raise
`--min_chi` to 2000 and keep the four splits that clear it, or keep the rest
and label them marginal.

### 9. Confirm, then interpret

```
spectrumSplits.py --input_tree tree.pruned.masked.pb.gz --fold \
    --min_chi 2000 --min_mutations 2500 --bootstrap_splits 100 \
    --nthreads 8 --output_spectrum spectra.tsv
```

Then read each split's per-class proportions out of `spectra.tsv`. A split
driven almost entirely by one mutation class is usually an artifact or a single
over-mutated site; one spread across several classes that survives the
bootstrap is the thing worth explaining. Running `annotate_sites.py --highlight
repair_genes.txt` on the mutations of the split's defining branch is the first
place to look for a mechanism.
