#!/bin/bash
#SBATCH --job-name=04_joint_genotyping
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00
#SBATCH --partition=qib-compute

# Step 4: Joint genotyping — GenomicsDBImport + GenotypeGVCFs
# Single job. Run after all per-sample GVCFs are ready (clinical + public).
# Input:  results/03_gvcfs/*.g.vcf.gz  +  public_data/gvcfs/*.g.vcf.gz
# Output: results/04_joint/all_samples.vcf.gz (raw, unfiltered)
#
# GenomicsDBImport behaviour:
#   - First run:    builds DB from all GVCFs found in GVCF dirs
#   - Subsequent:   appends only samples not already in the DB, then re-genotypes

# ===== SETUP =====

source "${SCRIPTS_DIR}/config.conf"

GVCF_DIR="${RESULTS_DIR}/03_gvcfs"
PUBLIC_GVCF_DIR="${PUBLIC_DATA_DIR}/gvcfs"
OUT_DIR="${RESULTS_DIR}/04_joint"
DB="${OUT_DIR}/genomicsdb"
TMP_DIR="${RESULTS_DIR}/tmp/04_joint_genotyping"
GENOTYPED_VCF="${OUT_DIR}/all_samples.vcf.gz"
SAMPLE_MAP="${OUT_DIR}/sample_map.tsv"

mkdir -p "$OUT_DIR" "$TMP_DIR"

# ===== BUILD SAMPLE MAP =====
# Tab-separated: sample_name <TAB> /path/to/sample.g.vcf.gz
# GenomicsDBImport uses this instead of repeated -V flags.

> "$SAMPLE_MAP"
for gvcf in "${GVCF_DIR}"/*.g.vcf.gz "${PUBLIC_GVCF_DIR}"/*.g.vcf.gz; do
    [[ -f "$gvcf" ]] || continue
    sample=$(basename "$gvcf" .g.vcf.gz)
    printf '%s\t%s\n' "$sample" "$gvcf" >> "$SAMPLE_MAP"
done

echo "=== Sample map: $(wc -l < "$SAMPLE_MAP") samples ==="

# ===== GenomicsDBImport =====
# If the DB already exists, identify samples not yet imported and append only those.
# If the DB does not exist, build it from scratch.

if [[ -d "$DB" ]]; then
    echo "=== GenomicsDBImport — appending new samples to existing DB ==="

    # Samples already in the DB
    EXISTING_SAMPLES=$(singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" \
        gatk GenomicsDBImport \
        --genomicsdb-workspace-path "$DB" \
        --print-interval-list 2>/dev/null | grep -v "^#" || true)

    # Build a map of only new samples
    NEW_SAMPLE_MAP="${OUT_DIR}/new_sample_map.tsv"
    > "$NEW_SAMPLE_MAP"
    while IFS=$'\t' read -r sample path; do
        if ! grep -qxF "$sample" <<< "$EXISTING_SAMPLES"; then
            printf '%s\t%s\n' "$sample" "$path" >> "$NEW_SAMPLE_MAP"
        fi
    done < "$SAMPLE_MAP"

    NEW_COUNT=$(wc -l < "$NEW_SAMPLE_MAP")
    if [[ "$NEW_COUNT" -eq 0 ]]; then
        echo "=== No new samples to append — skipping GenomicsDBImport ==="
    else
        echo "=== Appending $NEW_COUNT new sample(s) ==="
        singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk GenomicsDBImport \
            --genomicsdb-update-workspace-path "$DB" \
            --sample-name-map "$NEW_SAMPLE_MAP" \
            --intervals "$INTERVALS" \
            --reader-threads "$SLURM_CPUS_PER_TASK" \
            --tmp-dir "$TMP_DIR"
    fi

    rm -f "$NEW_SAMPLE_MAP"

else
    echo "=== GenomicsDBImport — building DB from scratch ($(wc -l < "$SAMPLE_MAP") samples) ==="
    singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk GenomicsDBImport \
        --genomicsdb-workspace-path "$DB" \
        --sample-name-map "$SAMPLE_MAP" \
        --intervals "$INTERVALS" \
        --reader-threads "$SLURM_CPUS_PER_TASK" \
        --tmp-dir "$TMP_DIR"
fi

# ===== GenotypeGVCFs =====
# Always re-runs — joint calling requires all samples to be re-genotyped together.

echo "=== GenotypeGVCFs ==="

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://${DB}" \
    -O "$GENOTYPED_VCF" \
    --sample-ploidy "$PLOIDY" \
    --tmp-dir "$TMP_DIR"

rm -rf "$TMP_DIR"

echo "=== Done: $GENOTYPED_VCF ==="