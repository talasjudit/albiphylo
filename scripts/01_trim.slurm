#!/bin/bash
#SBATCH --job-name=01_trim
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-04:00
#SBATCH --partition=qib-compute

# Step 1: Trimming + QC with fastp
# Per-sample array job. Reads R1/R2 paths from SAMPLE_LIST_FILE.
# Output: trimmed FASTQs + fastp HTML/JSON reports - results/01_trimmed/

# ===== SETUP =====

set -euo pipefail

source "${SCRIPTS_DIR}/config.conf"

# Read this array task's sample from the list
# Format: sample_id,R1_path,R2_path
IFS=',' read -r SAMPLE_ID R1 R2 < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST_FILE")

OUT_DIR="${RESULTS_DIR}/01_trimmed/${SAMPLE_ID}"
mkdir -p "$OUT_DIR"

echo "=== fastp trim: $SAMPLE_ID ==="
echo "R1: $R1"
echo "R2: $R2"

# ===== FASTP =====

singularity exec "${SIF_DIR}/fastp_1.1.0--heae3180_0.sif" fastp \
    --in1 "$R1" \
    --in2 "$R2" \
    --out1 "${OUT_DIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz" \
    --out2 "${OUT_DIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz" \
    --html "${OUT_DIR}/${SAMPLE_ID}_fastp.html" \
    --json "${OUT_DIR}/${SAMPLE_ID}_fastp.json" \
    --thread "$SLURM_CPUS_PER_TASK" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50

echo "=== Done: $SAMPLE_ID ==="
