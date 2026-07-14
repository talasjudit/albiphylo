#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-01:00
#SBATCH --partition=qib-compute

# MultiQC — aggregate fastp reports
# Single job. Run after step 1 (trimming) to review QC before proceeding.
# Output: results/multiqc/multiqc_report.html

# ===== SETUP =====

set -euo pipefail

source "${SCRIPTS_DIR}/config.conf"

OUT_DIR="${RESULTS_DIR}/multiqc"
mkdir -p "$OUT_DIR"

echo "=== MultiQC ==="

singularity exec "${SIF_DIR}/multiqc_v1.33.sif" multiqc \
    "${RESULTS_DIR}/01_trimmed/" \
    --outdir "$OUT_DIR" \
    --filename "multiqc_report" \
    --force

echo "=== Done: ${OUT_DIR}/multiqc_report.html ==="
