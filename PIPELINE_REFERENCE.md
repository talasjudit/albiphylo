# Pipeline Reference — Step-by-Step

---

## Step 1: Trimming with fastp

**Type:** SLURM array job (one task per sample)

**Parameters:**
- `--detect_adapter_for_pe` — auto-detect PE adapters
- `--qualified_quality_phred 20` — bases below Q20 are low quality
- `--length_required 50` — discard reads shorter than 50 bp after trimming

---

## MultiQC

**Type:** Single SLURM job.

**Purpose:** QC checkpoint after trimming. Review the aggregated HTML report before
proceeding to mapping. To exclude a failed sample from downstream analysis, remove
its row from `samples.tsv` (or blank the R1/R2 paths).

---

## Step 2: Mapping with bwa-mem

**Type:** SLURM array job (one task per sample)

**Parameters:**
- Read groups parsed from FASTQ header (flowcell, lane, barcode)
- `@RG` tags: ID=flowcell.lane, SM=sample_id, PL=ILLUMINA, LB=sample_id, PU=flowcell.lane.barcode

**Gotchas:**
- Read group parsing assumes standard Illumina FASTQ header format. Non-Illumina reads would need manual RG specification.
- Intermediate sorted BAM is deleted after MarkDuplicates — only the markdup BAM is kept.

---

## Step 3: HaplotypeCaller

**Type:** SLURM array job (one task per sample)

**Parameters:**
- `-ERC GVCF` — emit reference confidence (needed for joint calling)
- `--sample-ploidy 2` — diploid
- `--native-pair-hmm-threads` — uses all allocated CPUs

**GVCF is the durable checkpoint.** New samples only need steps 1-3; steps 4+ re-run on the combined GVCF set. runjob.sh auto-skips samples with existing GVCFs unless `status=rerun`.

**Post-processing for public cohort:** After step 3 completes for public samples, move GVCFs to `public_data/gvcfs/` so they survive `results/` cleanup and can be shared across projects.

---

## Step 4: Joint Genotyping

**Type:** Single SLURM job

**Parameters:**
- `--sample-ploidy 2` — diploid
- `--intervals` — interval list from setup_reference.sh (restricts to non-excluded chromosomes)
- `--reader-threads` — uses all allocated CPUs

**Behaviour:**
- First run: builds GenomicsDB from scratch
- Subsequent runs: auto-detects and appends only new samples, then re-genotypes all
- GenotypeGVCFs always re-runs (joint calling requires all samples)

**How it works:**
1. **GenomicsDBImport** — builds a database (matrix of genomic positions × samples, storing genotype likelihoods). The alternative was CombineGVCFs which re-reads all GVCFs and produces a combined GVCF every time before genotyping.
2. **GenotypeGVCFs** — walks through every position in the database, makes final genotype calls considering all samples simultaneously. This step re-runs every time new samples are added because allele frequencies change.

---

## Step 5: Filtering (two-stage)

**Type:** Single SLURM job

Because we skip BQSR, we use a two-stage filtering approach: Adamu Bukari et al. 2025
site-level hard filters followed by Anderson et al. 2023 per-genotype filters.

### Stage 1 — Adamu Bukari site-level hard filters

Standard GATK recommendations; catch poor-quality variant calls regardless of BQSR status.

| Filter | Threshold | What it catches |
|--------|-----------|-----------------|
| QD | < 2.0 | Low-confidence calls relative to depth |
| FS | > 60.0 | Strand bias (Fisher's exact test) |
| MQ | < 30.0 | Low mapping quality |
| ReadPosRankSum | < -8.0 | Alt alleles cluster at read ends (artefact) |

Then exclude variants in repetitive/problematic regions via the exclusion BED:

- Subtelomeres: 15 kb from each chromosome end
- Centromeres: from CGD GFF annotations
- MRS (major repeat sequences): HOK, RPS, RB2, Ca3 from CGD GFF
- chrM: excluded entirely (mitochondrial)
- chrR: NOT excluded (carries real LOH/CNV signal)

### Stage 2 — Anderson per-genotype filters (compensates for no BQSR)

Anderson demonstrated these filters produce high-quality variant sets (112,136 SNVs
across 431 samples) without BQSR. They're stricter than the Adamu Bukari filters but
apply per-genotype rather than per-site.

| Filter | Threshold | Effect |
|--------|-----------|--------|
| GQ (per-genotype) | < 20 (Phred) | Set individual genotype calls with prob < 0.99 to no-call (`./.`) |
| MQRankSum (site-level) | != 0.0 | Remove sites with any ref/alt mapping quality imbalance (subsumes Adamu Bukari's MQRS < -12.5) |
| Missing genotype fraction | > 5% | Remove sites where more than 5% of samples are no-call after the GQ filter |

**Gotchas:**
- Config.conf thresholds **must be floats** (e.g. `2.0` not `2`). GATK's JEXL engine infers the comparison type from the literal — integer literals cause `NumberFormatException` when comparing to float VCF fields like FS=8.451.
- Sites without an `MQRankSum` annotation (e.g. all-hom-ref loci) are kept by the stage-2 filter — only sites with `MQRankSum` present AND ≠ 0 are removed.

---

## Step 6: Global Phylogeny (FastTree)

**Type:** Single SLURM job

Builds a maximum-likelihood tree from all 190 samples (8 clinical + 182 Ropars backbone).
Used for clade assignment and global context.

**Sub-steps:**

1. **SelectVariants** — extract biallelic SNPs only (`--select-type-to-include SNP --restrict-alleles-to BIALLELIC`). Multi-allelic sites violate the assumption of the GTR substitution model, and indels can't be represented in a single-character FASTA alignment.
2. **vcf2phylip** — convert VCF to FASTA alignment
   - `--min-samples-locus 4` — a site needs genotype data in at least 4/190 samples
   - Heterozygous genotypes → IUPAC ambiguity codes (R, Y, S, W, K, M)
   - Missing genotypes (./.) → N
3. **FastTree** — maximum likelihood tree with `-gtr -gamma -nt`
   - GTR+gamma model (matches Adamu Bukari et al. 2025)
   - Double precision is default in FastTree 2.2.0
   - SH-like branch support values (not bootstrap)

**Gotchas:**
- vcf2phylip names output files as `<prefix>.min<N>.fasta` based on the --min-samples-locus value. The script renames this to `all_samples_snps.fasta`.

**Rooting:** The tree is unrooted when produced. In iTOL, root on the C. africana +
clade 13 outgroup — click the branch leading to that cluster → "Re-root the tree here".
Follows Adamu Bukari et al. 2025 ("C. albicans clade 13, i.e., C. africana").

---

## Step 7: Clinical Phylogeny (RAxML-NG)

**Type:** Single SLURM job

Higher-resolution tree with bootstrap support for clinical isolates only. Subsets the
filtered VCF to clinical samples, re-converts to FASTA, then runs RAxML-NG.

**Parameters:**
- `--all` — combined ML search + bootstrap convergence test
- `--model GTR+G` — general time-reversible + gamma (matches Adamu Bukari 2025)
- `--bs-trees 100` — 100 bootstrap replicates (sufficient for 8 isolates; increase for publication)
- `--seed 12345` — reproducibility

**Rationale:** Adamu Bukari et al. 2025 used FastTree for the global tree and RAxML for
within-group trees, citing the higher branch support rigour of RAxML bootstrap values
over FastTree's SH-like supports.

---

## Step 8: Pairwise SNP Distances

**Type:** Single SLURM job

Computes pairwise SNP distances between clinical isolates following Anderson et al. 2023.

**Method:**
1. Subset the filtered VCF to clinical samples (GATK SelectVariants)
2. Extract genotype table (GATK VariantsToTable, `-GF GT`)
3. Python helper (`pairwise_snp_distance.py`) computes:
   - **Dosage coding:** 0/0 → 0, 0/1 → 1, 1/1 → 2, ./. → missing
   - **Per-site distance:** |dosage_i - dosage_j| / 2
   - **Pair total:** sum of per-site distances / number of non-missing sites compared
   - **Raw count:** number of sites where genotypes differed at all

**Outputs:**
- `clinical_distance_matrix.tsv` — normalised distance (Anderson formula)
- `clinical_snp_counts.tsv` — raw count of differing SNPs per pair

No participant groupings required — the matrix itself shows relatedness; the client can
overlay longitudinal/body-site metadata when visualising.

---

## Step 9: Summary Reporting

**Type:** Single SLURM job

Collates per-isolate results into final client-facing tables and updates iTOL annotations.

**Sub-steps:**

1. **Clade assignment** (`assign_clades.py`) — for each clinical isolate, computes pairwise
   distances to all 182 backbone samples, assigns the clade of the nearest neighbor, and
   reports top-3 consensus. Flagged if top 3 neighbors disagree on clade (indicates a
   sample sitting near a clade boundary).

2. **Summary statistics** — per-isolate SNP count in the filtered VCF, mean sequencing depth
   from the markdup BAM, assigned clade.

3. **iTOL annotation regeneration** (`generate_itol_annotations.py`) — with clade assignments
   now available, clinical isolates get their assigned clade colour in the colorstrip.

**Outputs:**
- `clinical_clade_assignments.tsv`
- `clinical_summary_stats.tsv`
- Updated `results/06_phylogeny/itol/*.txt` (clinical isolates now coloured by clade)

---

## LOH/CNV (manual — YMAP)

Manual BAM upload to [YMAP](http://lovelace.cs.umn.edu/Ymap/). Appropriate for current
small sample sets.

**Requires client input:**
- At what sample size should we switch to a scripted local solution? YMAP becomes impractical above ~30 samples.

---

## Where our methods differ from the reference papers

### vs Adamu Bukari et al. 2025

| Step | Adamu Bukari | This pipeline | Reason |
|------|-------------|---------------|--------|
| BQSR | Used BQSR with CGD known polymorphisms VCF | **Skipped** | CGD VCF is SC5314's own het sites, not a population truth set. For non-SC5314 strains, real variants would be treated as errors. Anderson et al. 2023 also skipped BQSR with good results. |
| Filtering | Site-level hard filters only (QD, FS, MQ, MQRS, RPRS) | **Two-stage: Adamu Bukari site-level + Anderson per-genotype** | Because we skip BQSR, we add Anderson's stricter per-genotype filtering (GQ < 20 → no-call, MQRS = 0.0) to compensate. |
| VCF→FASTA | Not stated in paper | vcf2phylip (Ortiz 2019) | Paper doesn't specify the tool. snp-sites cannot parse GATK diploid VCF. vcf2phylip is the standard tool for diploid GATK output. |
| Exclusion BED source | Coordinates in Table S3 | Derived from CGD GFF | Same approach (CGD annotations), but generated programmatically from the GFF. |
| Joint genotyping | CombineGVCFs (implied) | **GenomicsDBImport** | More efficient at ~190 samples; supports incremental append of new samples. Same GenotypeGVCFs step afterward. |
| Phylogeny backbone | Ropars et al. 2018 (182 isolates) | Same | Identical backbone dataset. |
| Within-group trees | RAxML (per-participant) | **RAxML on all clinical isolates together** | We don't have participant groupings. The clinical RAxML tree + pairwise SNP matrix are sufficient for strain-level comparison. |

### vs Anderson et al. 2023

| Step | Anderson | This pipeline | Reason |
|------|----------|---------------|--------|
| Sequencing | TELL-Seq linked reads | Illumina short reads | Different input data type. Our pipeline maps with BWA-MEM (same as Anderson). |
| BQSR | Skipped | Skipped | Same approach. |
| Genotype filters | Per-genotype: GQ >= 0.99, MQRS = 0.0; site-level missing ≤ 5% | **Same, applied as stage 2 after Adamu Bukari hard filters** | We combine Anderson's per-genotype filters with Adamu Bukari's site-level hard filters for best of both. |
| Pairwise distances | Dosage coding, |d_i - d_j|/2, normalised by non-missing sites | **Same** | Directly adopted Anderson's method for step 8. |
| Tree method | Neighbour-joining on custom distance matrix | **FastTree ML (global) + RAxML-NG (clinical)** | ML is more statistically rigorous. Anderson used NJ for distance-based comparisons; we use ML for phylogenetic resolution. |
| Structural variants | TELL-Seq assemblies + MiniMap2 + CHEF karyotyping | Not included (optional ONT track available) | Different scope. Could add if hybrid assemblies are needed. |

---

## Methods notes for publication

- **BQSR omission** justified by Anderson et al. 2023 (successful without BQSR) and the non-applicability of the CGD known variants VCF as a population-level truth set.
- **Two-stage filtering** combines Adamu Bukari (standard GATK site-level filters, 2025) with Anderson (per-genotype GQ + MQRS = 0.0, 2023) — cite both.
- **Clade assignment** is programmatic: nearest-neighbor to Ropars 2018 backbone by Anderson-style pairwise dosage distance. Ropars clade metadata joined with SRA run info (PRJNA432884) on CEC strain IDs.
- **Outgroup rooting** follows Adamu Bukari 2025 using C. africana (= clade 13 sensu Ropars).
- **vcf2phylip** should be cited as Ortiz 2019 (doi:10.1111/1755-0998.13115) for VCF→FASTA conversion.
- **GenomicsDBImport vs CombineGVCFs** is an implementation choice that does not affect the variant calls — the GenotypeGVCFs step is identical.

---

## Summary of parameters to confirm with client

| Parameter | Current value | Source | Question |
|-----------|--------------|--------|----------|
| Subtelomere buffer | 15 kb | Adamu Bukari 2025 | Not justified in paper — appropriate for this dataset? |
| Bootstrap replicates (step 7) | 100 | Default | Sufficient for publication, or increase to 1000? |
| `--min-samples-locus` | 4 | vcf2phylip default | Appropriate threshold? |
| LOH/CNV method | YMAP (manual) | Adamu Bukari 2025 | Switch to scripted at what sample size? |
| Pixy (nucleotide diversity) | Not run | — | Needs all-sites VCF — add before re-running step 4? |
