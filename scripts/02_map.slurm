#!/bin/bash
#SBATCH --job-name=02_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-08:00
#SBATCH --partition=qib-compute

# Step 2: BWA-MEM mapping → samtools sort → GATK MarkDuplicates
# Per-sample array job.
# Input:  results/01_trimmed/<sample>/
# Output: results/02_mapped/<sample>.markdup.bam + .bai

# ===== SETUP =====

source "${SCRIPTS_DIR}/config.conf"

IFS=',' read -r SAMPLE_ID R1 R2 < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST_FILE")

TRIM_DIR="${RESULTS_DIR}/01_trimmed/${SAMPLE_ID}"
OUT_DIR="${RESULTS_DIR}/02_mapped"
TMP_DIR="${RESULTS_DIR}/tmp/02_map_${SAMPLE_ID}"
mkdir -p "$OUT_DIR" "$TMP_DIR"

TRIM_R1="${TRIM_DIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz"
TRIM_R2="${TRIM_DIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz"
BAM_SORTED="${OUT_DIR}/${SAMPLE_ID}.sorted.bam"
BAM_MARKDUP="${OUT_DIR}/${SAMPLE_ID}.markdup.bam"
METRICS="${OUT_DIR}/${SAMPLE_ID}.markdup_metrics.txt"

echo "=== BWA-MEM map: $SAMPLE_ID ==="

# ===== Parse read group info from FASTQ header =====
# Illumina header: @instrument:run:flowcell:lane:... read:filter:flag:barcode
# ID  = flowcell.lane        — unique per sequencing lane
# PU  = flowcell.lane.barcode — unique per library per lane (used by BQSR if ever enabled)
# SM  = sample name          — appears as VCF column header in GATK output
# PL  = platform
# LB  = library              — used by MarkDuplicates to group duplicates correctly

HEADER=$(gunzip -c "$TRIM_R1" | head -1)
FLOWCELL=$(echo "$HEADER" | cut -d: -f3)
LANE=$(echo "$HEADER" | cut -d: -f4)
BARCODE=$(echo "$HEADER" | cut -d' ' -f2 | cut -d: -f4)
RG_ID="${FLOWCELL}.${LANE}"
RG_PU="${FLOWCELL}.${LANE}.${BARCODE}"

echo "Read group: ID=${RG_ID} PU=${RG_PU} SM=${SAMPLE_ID}"

# ===== BWA-MEM + samtools sort =====

singularity exec "${SIF_DIR}/bwa_0.7.19.sif" bwa mem \
    -t "$SLURM_CPUS_PER_TASK" \
    -R "@RG\tID:${RG_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:${SAMPLE_ID}\tPU:${RG_PU}" \
    "$REF" \
    "$TRIM_R1" "$TRIM_R2" \
    | singularity exec "${SIF_DIR}/samtools_1.21.sif" samtools sort \
        -@ "$SLURM_CPUS_PER_TASK" \
        -o "$BAM_SORTED"

singularity exec "${SIF_DIR}/samtools_1.21.sif" samtools index "$BAM_SORTED"

# ===== GATK MarkDuplicates =====

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk MarkDuplicates \
    -I "$BAM_SORTED" \
    -O "$BAM_MARKDUP" \
    -M "$METRICS" \
    --TMP_DIR "$TMP_DIR"

singularity exec "${SIF_DIR}/samtools_1.21.sif" samtools index "$BAM_MARKDUP"

# Clean up intermediate sorted BAM and tmp
rm "$BAM_SORTED" "${BAM_SORTED}.bai"
rm -rf "$TMP_DIR"

echo "=== Done: $SAMPLE_ID ==="
