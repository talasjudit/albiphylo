# C. albicans Clinical Isolate Genomic Pipeline

Reference-based variant calling, phylogenetic placement, clade assignment, and pairwise
SNP comparison for *C. albicans* clinical isolates. GVCF checkpoint architecture —
new samples only need per-sample steps; downstream analyses re-run on the full set.

- Reference: SC5314 Assembly 22 (A22-s08-m01-r36)
- Public backbone: Ropars et al. 2018, 182 isolates (PRJNA432884)

---

## Samplesheet setup

The pipeline reads clinical sample information from `metadata/samples.tsv` (tab-separated).

### Columns

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | Yes | Unique identifier — used in all output filenames and VCF headers |
| `status` | No | Leave blank for normal processing; set to `rerun` to force GVCF regeneration |
| `illumina_irida_project` | No | IRIDA project ID (traceability only) |
| `illumina_irida_sample` | No | IRIDA sample ID (traceability only) |
| `illumina_r1` | Yes | Absolute path to Illumina R1 FASTQ |
| `illumina_r2` | Yes | Absolute path to Illumina R2 FASTQ |
| `ont_irida_project` | No | IRIDA project ID for ONT (not used by core pipeline) |
| `ont_irida_sample` | No | IRIDA sample ID for ONT (not used by core pipeline) |
| `ont_reads` | No | Absolute path to ONT FASTQ (not used by core pipeline) |
| `assembly` | No | Absolute path to assembly FASTA (not used by core pipeline) |

### Status field

The status column has one meaningful value: `rerun`. Any other value (including blank) is
treated as "process normally".

| Status | Behaviour |
|--------|-----------|
| *(blank)* | Process normally. Step 3 skips if a GVCF already exists. |
| `rerun` | Force step 3 to regenerate the GVCF even if one exists. |

**To exclude a sample:** either remove the row entirely, or leave `illumina_r1`/`illumina_r2`
blank — the sample will be skipped in steps 1-3 automatically. Steps 7-9 use all rows
with a `sample_id`, so removing the row is the cleanest way to exclude.

### Example

```
sample_id	status	illumina_irida_project	illumina_irida_sample	illumina_r1	illumina_r2
NCYC_1470		2570	PID-2045-1470	/path/to/R1.fastq.gz	/path/to/R2.fastq.gz
NCYC_1471	rerun	3027	PID-2462-1471	/path/to/R1.fastq.gz	/path/to/R2.fastq.gz
```

- Paths must be **absolute**
- Save as tab-separated (not comma-separated). If editing in Excel, watch for CRLF line
  endings — the pipeline strips `\r` but best to avoid

---

## Running the pipeline

All commands assume you are in the `scripts/` directory:

```bash
cd scripts/
```

### One-time setup

Copy the config template and edit it to match your environment — replace
`<PROJECT_ROOT>` with the absolute path to your project directory, and set
`SLURM_PARTITION` to your HPC partition name.

```bash
cp config.conf.example config.conf
# edit config.conf (e.g. with vim, nano, or your editor of choice)

bash install.sh          # pull Singularity images
bash setup_reference.sh  # index reference, build exclusion BED + interval list
```

The local `config.conf` is gitignored so environment-specific paths and
notification settings stay private.

### Process the public backbone (once)

The 182 Ropars et al. isolates form the phylogenetic backbone:

```bash
./runjob.sh 01_trim.slurm               --cohort public
./runjob.sh 02_map.slurm                --cohort public
./runjob.sh 03_haplotypecaller.slurm    --cohort public
```

Then move public GVCFs to their permanent location so they survive any `results/` cleanup:

```bash
mkdir -p ../public_data/gvcfs
mv ../results/03_gvcfs/SRR*.g.vcf.gz* ../public_data/gvcfs/
```

### Process clinical samples

```bash
./runjob.sh 01_trim.slurm
./runjob.sh multiqc.slurm               # optional QC checkpoint after trimming
./runjob.sh 02_map.slurm
./runjob.sh 03_haplotypecaller.slurm
```

### Joint genotyping, filtering, phylogeny, and reporting

Run once all GVCFs are ready (clinical in `results/03_gvcfs/`, public in `public_data/gvcfs/`):

```bash
./runjob.sh 04_joint_genotyping.slurm    # GenomicsDBImport + GenotypeGVCFs
./runjob.sh 05_filter.slurm              # Adamu Bukari hard filters + Anderson per-genotype filters
./runjob.sh 06_phylogeny_fasttree.slurm  # global tree (all 190 samples)
./runjob.sh 07_phylogeny_raxml.slurm     # clinical-only RAxML tree with bootstraps
./runjob.sh 08_pairwise_snps.slurm       # pairwise SNP distance matrix (clinical)
./runjob.sh 09_summary.slurm             # clade assignment + summary stats + iTOL annotations
```

### Tree visualisation

1. Upload `results/06_phylogeny/global_tree.nwk` to [iTOL](https://itol.embl.de/)
2. Drag-and-drop the files in `results/06_phylogeny/itol/` onto the tree:
   - `itol_clade_colorstrip.txt` — clade colours
   - `itol_clinical_labels.txt` — bold red labels on clinical isolates
   - `itol_popup_info.txt` — hover popups
3. Root the tree on the C. africana + clade 13 outgroup (black + magenta in the colorstrip):
   click the branch leading to that cluster → "Re-root the tree here"

---

## Adding new samples

1. Add rows to `metadata/samples.tsv` with Illumina read paths (status blank)
2. Run steps 1-3 — existing samples with GVCFs are auto-skipped
3. Re-run steps 4-9 — GenomicsDBImport auto-appends new samples; all downstream analyses rebuild with the full set

Set `status=rerun` on any sample whose GVCF you want regenerated (e.g. after a reference change).

---

## Pipeline steps at a glance

| Step | Script | Type | Tool(s) | Input | Output |
|------|--------|------|---------|-------|--------|
| 1 | `01_trim.slurm` | Array | fastp 1.1.0 | Raw FASTQs | Trimmed FASTQs + QC reports |
| QC | `multiqc.slurm` | Single | MultiQC v1.33 | fastp JSON reports | Aggregated QC HTML |
| 2 | `02_map.slurm` | Array | BWA-MEM, samtools, GATK MarkDuplicates | Trimmed FASTQs | `<sample>.markdup.bam` |
| 3 | `03_haplotypecaller.slurm` | Array | GATK HaplotypeCaller | markdup BAM | `<sample>.g.vcf.gz` (GVCF) |
| 4 | `04_joint_genotyping.slurm` | Single | GATK GenomicsDBImport + GenotypeGVCFs | All GVCFs | `all_samples.vcf.gz` |
| 5 | `05_filter.slurm` | Single | GATK VariantFiltration + bedtools | Joint VCF | `all_samples.filtered.vcf` |
| 6 | `06_phylogeny_fasttree.slurm` | Single | GATK SelectVariants, vcf2phylip, FastTree | Filtered VCF | `global_tree.nwk` |
| 7 | `07_phylogeny_raxml.slurm` | Single | GATK SelectVariants, vcf2phylip, RAxML-NG | Filtered VCF (clinical subset) | `clinical.raxml.support` (ML tree with bootstrap values) |
| 8 | `08_pairwise_snps.slurm` | Single | GATK SelectVariants + VariantsToTable, Python | Filtered VCF (clinical subset) | `clinical_distance_matrix.tsv`, `clinical_snp_counts.tsv` |
| 9 | `09_summary.slurm` | Single | GATK VariantsToTable, Python | Filtered VCF + pairwise outputs | `clinical_clade_assignments.tsv`, `clinical_summary_stats.tsv`, updated iTOL annotations |

Steps 1-3 are per-sample array jobs submitted via `runjob.sh`.
Steps 4-9 are single jobs that operate on all samples (or the clinical subset) together.

---

## Outputs

After a full pipeline run, the key client-facing outputs are:

| Output | File | Description |
|--------|------|-------------|
| Global phylogeny | `results/06_phylogeny/global_tree.nwk` + iTOL annotations | 190-sample tree for clade context |
| Clinical phylogeny | `results/07_raxml/clinical.raxml.support` | Higher-resolution tree with bootstrap support |
| Pairwise SNP matrix | `results/08_pairwise/clinical_snp_counts.tsv` | Raw SNP differences between every clinical pair |
| Normalised distances | `results/08_pairwise/clinical_distance_matrix.tsv` | Anderson-style dosage distances |
| Clade assignments | `results/09_summary/clinical_clade_assignments.tsv` | Per-isolate clade + nearest backbone neighbor + top-3 consensus |
| Summary stats | `results/09_summary/clinical_summary_stats.tsv` | SNP counts, mean depth, assigned clade per isolate |

See `PIPELINE_REFERENCE.md` for detailed methodology and parameter rationale.
