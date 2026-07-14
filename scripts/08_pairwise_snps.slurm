#!/bin/bash
#SBATCH --job-name=08_pairwise
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-02:00
#SBATCH --partition=qib-compute

# Step 8: Pairwise SNP distance matrix for clinical isolates
# Single job. Rerunnable — add new isolates to samples.tsv, re-run steps 4-5, then this.
# Computes per-pair distances following Anderson et al. 2023:
#   dosage coding (0/0=0, 0/1=1, 1/1=2), distance = |d_i - d_j| / 2, normalised
#   by number of non-missing sites.
# No participant groupings needed — the matrix itself shows relatedness;
# downstream visualisation/interpretation left to the client.
# Input:  results/05_filtered/all_samples.filtered.vcf + metadata/samples.tsv
# Output: results/08_pairwise/clinical_distance_matrix.tsv  (normalised distances)
#         results/08_pairwise/clinical_snp_counts.tsv        (raw differing SNP counts)

# ===== SETUP =====

set -euo pipefail

source "${SCRIPTS_DIR}/config.conf"

OUT_DIR="${RESULTS_DIR}/08_pairwise"
mkdir -p "$OUT_DIR"

FILTERED_VCF="${RESULTS_DIR}/05_filtered/all_samples.filtered.vcf"

# ===== Extract clinical sample IDs =====

SAMPLE_FILE="${OUT_DIR}/clinical_samples.txt"
awk -F'\t' 'NR>1 && $1!="" { print $1 }' "$SAMPLES_TSV" > "$SAMPLE_FILE"

N_SAMPLES=$(wc -l < "$SAMPLE_FILE")
echo "=== $N_SAMPLES clinical samples ==="
cat "$SAMPLE_FILE"

# ===== Subset VCF to clinical samples =====

CLINICAL_VCF="${OUT_DIR}/clinical_snps.vcf"

echo "=== Subsetting VCF to clinical samples ==="

# GATK SelectVariants has no --sample-file; build repeated -sn flags instead.
SN_ARGS=()
while IFS= read -r s; do
    [[ -n "$s" ]] && SN_ARGS+=(-sn "$s")
done < "$SAMPLE_FILE"

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk SelectVariants \
    -R "$REF" \
    -V "$FILTERED_VCF" \
    -O "$CLINICAL_VCF" \
    "${SN_ARGS[@]}" \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    --remove-unused-alternates \
    --exclude-non-variants

N_SNPS=$(grep -cv '^#' "$CLINICAL_VCF")
echo "=== $N_SNPS biallelic SNPs ==="

# ===== Extract genotype table =====
# GATK VariantsToTable outputs a flat TSV with one row per variant,
# columns: CHROM, POS, <sample>.GT for each sample.

GT_TABLE="${OUT_DIR}/clinical_genotypes.tsv"

echo "=== Extracting genotype table ==="

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk VariantsToTable \
    -R "$REF" \
    -V "$CLINICAL_VCF" \
    -F CHROM -F POS -F REF -F ALT \
    -GF GT \
    -O "$GT_TABLE"

echo "=== $(tail -n +2 "$GT_TABLE" | wc -l) variants in genotype table ==="

# ===== Compute pairwise distances =====

echo "=== Computing pairwise SNP distances (Anderson et al. 2023 method) ==="

singularity exec "${SIF_DIR}/vcf2phylip-2.9.sif" python3 \
    "${SCRIPTS_DIR}/pairwise_snp_distance.py" \
    "$GT_TABLE" \
    "${OUT_DIR}/clinical"

echo "=== Done ==="
echo "Distance matrix: ${OUT_DIR}/clinical_distance_matrix.tsv"
echo "SNP counts:      ${OUT_DIR}/clinical_snp_counts.tsv"
