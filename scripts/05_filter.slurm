#!/bin/bash
#SBATCH --job-name=05_filter
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0-04:00
#SBATCH --partition=qib-compute

# Step 5: Variant filtering — site-level hard filters + Anderson per-genotype filters
# Single job. Two-stage approach because we skip BQSR:
#
# Stage 1 (Adamu Bukari et al. 2025): site-level hard filters
#   QD < 2, FS > 60, MQ < 30, RPRS < -8
#   + exclusion BED (subtelomeres, centromeres, MRS, chrM)
#
# Stage 2 (Anderson et al. 2023): per-genotype quality + strict site filters
#   Compensates for no BQSR — Anderson used this approach on 431 samples without BQSR.
#   - Set genotypes with GQ < 20 (Phred; prob < 0.99) to no-call
#   - Remove sites where MQRankSum ≠ 0.0 (only keep perfect mapping quality balance)
#   - Remove sites with > 5% missing genotypes (after GQ filtering)
#
# Input:  results/04_joint/all_samples.vcf.gz
# Output: results/05_filtered/all_samples.filtered.vcf

# ===== SETUP =====

source "${SCRIPTS_DIR}/config.conf"

IN_VCF="${RESULTS_DIR}/04_joint/all_samples.vcf.gz"
OUT_DIR="${RESULTS_DIR}/05_filtered"
mkdir -p "$OUT_DIR"

# Intermediate files
FLAGGED_VCF="${OUT_DIR}/all_samples.flagged.vcf.gz"
PASS_VCF="${OUT_DIR}/all_samples.pass.vcf.gz"
BED_FILTERED_VCF="${OUT_DIR}/all_samples.bed_filtered.vcf"
GQ_FLAGGED_VCF="${OUT_DIR}/all_samples.gq_flagged.vcf.gz"
GQ_NOCALL_VCF="${OUT_DIR}/all_samples.gq_nocall.vcf.gz"
MQRS_FLAGGED_VCF="${OUT_DIR}/all_samples.mqrs_flagged.vcf.gz"
MQRS_FILTERED_VCF="${OUT_DIR}/all_samples.mqrs_filtered.vcf.gz"
FILTERED_VCF="${OUT_DIR}/all_samples.filtered.vcf"

# ===================================================================
# STAGE 1: Site-level hard filters (Adamu Bukari et al. 2025)
# ===================================================================
# Standard GATK recommendations for site-level quality.
# These catch poor-quality variant calls regardless of BQSR status.

echo "=== Stage 1: Site-level hard filters ==="

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk VariantFiltration \
    -R "$REF" \
    -V "$IN_VCF" \
    -O "$FLAGGED_VCF" \
    --filter-expression "QD < ${FILTER_QD}"    --filter-name "QD_filter" \
    --filter-expression "FS > ${FILTER_FS}"    --filter-name "FS_filter" \
    --filter-expression "MQ < ${FILTER_MQ}"    --filter-name "MQ_filter" \
    --filter-expression "ReadPosRankSum < ${FILTER_RPRS}" --filter-name "RPRS_filter"

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk SelectVariants \
    -R "$REF" \
    -V "$FLAGGED_VCF" \
    -O "$PASS_VCF" \
    --exclude-filtered

BEFORE=$(singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk CountVariants -V "$IN_VCF" 2>/dev/null | tail -1)
AFTER_HARD=$(singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk CountVariants -V "$PASS_VCF" 2>/dev/null | tail -1)
echo "=== Hard filters: $BEFORE → $AFTER_HARD variants ==="

# ===== Apply exclusion BED =====
# Remove variants in subtelomeres (15kb), centromeres, MRS, chrM

echo "=== Exclusion BED ==="

singularity exec "${SIF_DIR}/bedtools_2.31.1.sif" bedtools intersect \
    -v \
    -a "$PASS_VCF" \
    -b "$EXCL_BED" \
    -header \
    > "$BED_FILTERED_VCF"

AFTER_BED=$(grep -cv '^#' "$BED_FILTERED_VCF")
echo "=== After BED exclusion: $AFTER_BED variants ==="

# ===================================================================
# STAGE 2: Anderson et al. 2023 per-genotype filters
# ===================================================================
# Since we skip BQSR, we compensate with stricter per-genotype filtering.
# Anderson demonstrated this produces high-quality variant sets (112,136 SNVs
# across 431 samples) without BQSR.

echo "=== Stage 2: Anderson per-genotype filters ==="

# ----- 2a: Per-genotype GQ filter -----
# Set individual genotype calls with GQ < 20 (prob < 0.99) to no-call (./.)
# This removes low-confidence calls at the genotype level, not the site level.

echo "=== GQ filter: setting GQ < ${FILTER_GQ_MIN} genotypes to no-call ==="

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk VariantFiltration \
    -R "$REF" \
    -V "$BED_FILTERED_VCF" \
    -O "$GQ_FLAGGED_VCF" \
    --genotype-filter-expression "GQ < ${FILTER_GQ_MIN}" \
    --genotype-filter-name "lowGQ"

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk SelectVariants \
    -R "$REF" \
    -V "$GQ_FLAGGED_VCF" \
    -O "$GQ_NOCALL_VCF" \
    --set-filtered-genotype-to-no-call

# ----- 2b: Site-level MQRankSum = 0.0 filter -----
# Anderson required MQRankSum to be exactly 0.0 ("the perfect score") —
# no mapping quality difference between ref and alt supporting reads.
# This is much stricter than Adamu Bukari's MQRS < -12.5 and subsumes it.
# Sites lacking MQRankSum annotation (e.g. all-hom-ref) are kept.

echo "=== MQRankSum == ${FILTER_MQRS_EXACT} filter ==="

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk VariantFiltration \
    -R "$REF" \
    -V "$GQ_NOCALL_VCF" \
    -O "$MQRS_FLAGGED_VCF" \
    --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum != ${FILTER_MQRS_EXACT}" \
    --filter-name "MQRS_ne_zero"

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk SelectVariants \
    -R "$REF" \
    -V "$MQRS_FLAGGED_VCF" \
    -O "$MQRS_FILTERED_VCF" \
    --exclude-filtered

AFTER_MQRS=$(singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk CountVariants -V "$MQRS_FILTERED_VCF" 2>/dev/null | tail -1)
echo "=== After MQRankSum filter: $AFTER_MQRS variants ==="

# ----- 2c: Missing genotype filter -----
# Remove sites where > 5% of samples have no-call genotypes.
# This catches sites degraded by the GQ filter above.

echo "=== Missing genotype filter: max ${FILTER_MAX_NOCALL} no-call fraction ==="

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk SelectVariants \
    -R "$REF" \
    -V "$MQRS_FILTERED_VCF" \
    -O "$FILTERED_VCF" \
    --max-nocall-fraction "$FILTER_MAX_NOCALL"

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk IndexFeatureFile -I "$FILTERED_VCF"

FINAL=$(grep -cv '^#' "$FILTERED_VCF")
echo ""
echo "=== Filtering summary ==="
echo "  Input variants:           $BEFORE"
echo "  After hard filters:       $AFTER_HARD"
echo "  After BED exclusion:      $AFTER_BED"
echo "  After MQRankSum == 0:     $AFTER_MQRS"
echo "  After missing genotypes:  $FINAL (FINAL)"
echo "=== Output: $FILTERED_VCF ==="

# ===== Clean up intermediates =====
rm -f "$FLAGGED_VCF" "${FLAGGED_VCF}.tbi"
rm -f "$PASS_VCF" "${PASS_VCF}.tbi"
rm -f "$GQ_FLAGGED_VCF" "${GQ_FLAGGED_VCF}.tbi"
rm -f "$GQ_NOCALL_VCF" "${GQ_NOCALL_VCF}.tbi"
rm -f "$MQRS_FLAGGED_VCF" "${MQRS_FLAGGED_VCF}.tbi"
rm -f "$MQRS_FILTERED_VCF" "${MQRS_FILTERED_VCF}.tbi"

echo "=== Done ==="
