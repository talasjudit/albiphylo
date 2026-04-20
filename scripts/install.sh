#!/bin/bash
# install.sh — pull/build all Singularity images for the C. albicans pipeline
# Run once from the project root before submitting any jobs.
# Usage: bash scripts/install.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SIF_DIR="$PROJECT_DIR/singularity/sif"
DEF_DIR="$PROJECT_DIR/singularity/def"

mkdir -p "$SIF_DIR"
cd "$SIF_DIR"

echo "=== C. albicans pipeline — Singularity image installer ==="
echo "Images will be written to: $SIF_DIR"
echo ""

# Helpers 

pull_image() {
    local name="$1"
    local source="$2"
    local sif="${name}.sif"
    if [[ -f "$sif" ]]; then
        echo "[SKIP]  $sif already exists"
    else
        echo "[PULL]  $source"
        singularity pull "$sif" "$source"
        echo "[DONE]  $sif"
    fi
}

# Pre-built Docker images

pull_image "gatk-4.6.2.0"      "docker://broadinstitute/gatk:4.6.2.0"
pull_image "fastp_1.1.0--heae3180_0"  "docker://biocontainers/fastp:1.1.0--heae3180_0"
pull_image "bwa_0.7.19"               "docker://staphb/bwa:0.7.19"
pull_image "samtools_1.21"            "docker://staphb/samtools:1.21"
pull_image "bedtools_2.31.1"          "docker://staphb/bedtools:2.31.1"
pull_image "multiqc_v1.33"            "docker://multiqc/multiqc:v1.33"
pull_image "fasttree_2.2.0"           "docker://staphb/fasttree:2.2.0"
pull_image "vcf2phylip-2.9"          "oras://ghcr.io/talasjudit/calbicancs/vcf2phylip:2.9-1"

echo ""
echo "=== Done ==="
ls -lh "$SIF_DIR"
