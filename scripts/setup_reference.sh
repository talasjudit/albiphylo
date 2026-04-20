#!/bin/bash
# setup_reference.sh — one-time reference preparation for C. albicans pipeline
# Run interactively or on any node with samtools, bwa, gatk (via Singularity), and bedtools available.
# Idempotent — skips steps whose outputs already exist.
# Usage: bash scripts/setup_reference.sh

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPTS_DIR/config.conf"

REF_DIR="$(dirname "$REF")"
GFF="${REF_DIR}/$(basename "$REF" _chromosomes.fasta)_features.gff"
FAI="${REF}.fai"
DICT="${REF%.fasta}.dict"
EXCL_BED_TMP_DIR="${REF_DIR}/tmp_bed"

echo "=== C. albicans pipeline — reference setup ==="
echo "Reference: $REF"
echo "GFF:       $GFF"
echo "Output BED: $EXCL_BED"
echo ""

# ── Validate inputs ───────────────────────────────────────────────────────────

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF"
    exit 1
fi

if [[ ! -f "$GFF" ]]; then
    echo "ERROR: GFF not found: $GFF"
    echo "Expected: $GFF"
    echo "Download from: https://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/"
    exit 1
fi

# ── samtools faidx ────────────────────────────────────────────────────────────

if [[ -f "$FAI" ]]; then
    echo "[SKIP] samtools faidx — $FAI exists"
else
    echo "[RUN]  samtools faidx"
    singularity exec "${SIF_DIR}/samtools_1.21.sif" samtools faidx "$REF"
    echo "[DONE] $FAI"
fi

# ── BWA index ─────────────────────────────────────────────────────────────────

if [[ -f "${REF}.bwt" ]]; then
    echo "[SKIP] bwa index — ${REF}.bwt exists"
else
    echo "[RUN]  bwa index (this takes ~30 seconds)"
    singularity exec "${SIF_DIR}/bwa_0.7.19.sif" bwa index "$REF"
    echo "[DONE] BWA index"
fi

# ── GATK CreateSequenceDictionary ─────────────────────────────────────────────

if [[ -f "$DICT" ]]; then
    echo "[SKIP] gatk CreateSequenceDictionary — $DICT exists"
else
    echo "[RUN]  gatk CreateSequenceDictionary"
    singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk CreateSequenceDictionary \
        -R "$REF"
    echo "[DONE] $DICT"
fi

# ── Interval list (one line per chromosome, for GenomicsDBImport) ─────────────

if [[ -f "$INTERVALS" ]]; then
    echo "[SKIP] interval list — $INTERVALS exists"
else
    echo "[RUN]  building interval list from $FAI"
    awk '{print $1}' "$FAI" > "$INTERVALS"
    echo "[DONE] $INTERVALS ($(wc -l < "$INTERVALS") intervals)"
fi

# ── Exclusion BED ─────────────────────────────────────────────────────────────

if [[ -f "$EXCL_BED" ]]; then
    echo "[SKIP] exclusion BED — $EXCL_BED exists"
else
    echo "[RUN]  building exclusion BED"

    mkdir -p "$EXCL_BED_TMP_DIR"

    # 1. Subtelomeres: 15 kb from each end of every sequence
    awk -v T=15000 'BEGIN{OFS="\t"} {
        print $1, 0, T;
        print $1, $2-T, $2
    }' "$FAI" > "${EXCL_BED_TMP_DIR}/subtelomeres.bed"

    # 2. Centromeres (GFF feature type: centromere; convert 1-based GFF to 0-based BED)
    awk 'BEGIN{OFS="\t"} $3=="centromere" {print $1, $4-1, $5}' "$GFF" \
        > "${EXCL_BED_TMP_DIR}/centromeres.bed"

    # 3. MRS regions (GFF feature type: repeat_region — all are MRS components in this assembly)
    awk 'BEGIN{OFS="\t"} $3=="repeat_region" {print $1, $4-1, $5}' "$GFF" \
        > "${EXCL_BED_TMP_DIR}/mrs.bed"

    # 4. chrM entirely (mitochondrial — not relevant for nuclear variant calling)
    grep "^Ca22chrM" "$FAI" | awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' \
        > "${EXCL_BED_TMP_DIR}/chrM.bed"

    # 5. Combine, sort, merge
    cat \
        "${EXCL_BED_TMP_DIR}/subtelomeres.bed" \
        "${EXCL_BED_TMP_DIR}/centromeres.bed" \
        "${EXCL_BED_TMP_DIR}/mrs.bed" \
        "${EXCL_BED_TMP_DIR}/chrM.bed" \
        | singularity exec "${SIF_DIR}/bedtools_2.31.1.sif" bedtools sort -i - \
        | singularity exec "${SIF_DIR}/bedtools_2.31.1.sif" bedtools merge -i - \
        > "$EXCL_BED"

    rm -rf "$EXCL_BED_TMP_DIR"

    echo "[DONE] $EXCL_BED ($(wc -l < "$EXCL_BED") regions)"
fi

echo ""
echo "=== Reference setup complete ==="
echo "  FAI:       $FAI"
echo "  BWA:       ${REF}.bwt (+ .amb .ann .pac .sa)"
echo "  Dict:      $DICT"
echo "  Intervals: $INTERVALS"
echo "  Excl BED:  $EXCL_BED"
