#!/bin/bash
#SBATCH --job-name=03_haplotypecaller
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-12:00
#SBATCH --partition=qib-compute

# Step 3: GATK HaplotypeCaller — per-sample GVCF
# Per-sample array job. GVCF is the durable checkpoint:
#   new samples only need steps 1-3; steps 4+ re-run on the combined GVCF set.
# Skipping logic handled by runjob.sh (checks for existing GVCF).
# Input:  results/02_mapped/<sample>.markdup.bam
# Output: results/03_gvcfs/<sample>.g.vcf.gz

# ===== SETUP =====

source "${SCRIPTS_DIR}/config.conf"

IFS=',' read -r SAMPLE_ID R1 R2 < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST_FILE")

BAM="${RESULTS_DIR}/02_mapped/${SAMPLE_ID}.markdup.bam"
OUT_DIR="${RESULTS_DIR}/03_gvcfs"
TMP_DIR="${RESULTS_DIR}/tmp/03_haplotypecaller_${SAMPLE_ID}"
mkdir -p "$OUT_DIR" "$TMP_DIR"
GVCF="${OUT_DIR}/${SAMPLE_ID}.g.vcf.gz"

echo "=== HaplotypeCaller: $SAMPLE_ID ==="

# ===== GATK HaplotypeCaller =====

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BAM" \
    -O "$GVCF" \
    -ERC GVCF \
    --sample-ploidy "$PLOIDY" \
    --native-pair-hmm-threads "$SLURM_CPUS_PER_TASK" \
    --tmp-dir "$TMP_DIR"

rm -rf "$TMP_DIR"

echo "=== Done: $SAMPLE_ID — GVCF: $GVCF ==="
