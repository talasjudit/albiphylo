#!/bin/bash
#SBATCH --job-name=09_summary
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0-04:00
#SBATCH --partition=qib-compute

# Step 9: Summary reporting — clade assignment, statistics, pairwise heatmap
# Single job. Produces the final output tables:
#   - Per-isolate clade assignment (nearest-neighbor from backbone)
#   - Summary statistics (coverage, SNP counts per sample)
#   - Pairwise SNP heatmap (R plot)
# Input:  results/05_filtered/all_samples.filtered.vcf
#         results/08_pairwise/clinical_snp_counts.tsv
#         public_data/ropars2018_clade_metadata.tsv
#         metadata/samples.tsv
# Output: results/09_summary/

# ===== SETUP =====

source "${SCRIPTS_DIR}/config.conf"

OUT_DIR="${RESULTS_DIR}/09_summary"
mkdir -p "$OUT_DIR"

FILTERED_VCF="${RESULTS_DIR}/05_filtered/all_samples.filtered.vcf"
CLADE_METADATA="${SCRIPTS_DIR}/../public_data/ropars2018_clade_metadata.tsv"
PAIRWISE_COUNTS="${RESULTS_DIR}/08_pairwise/clinical_snp_counts.tsv"

# ===== Extract clinical sample IDs =====

SAMPLE_FILE="${OUT_DIR}/clinical_samples.txt"
awk -F'\t' 'NR>1 && $1!="" { print $1 }' "$SAMPLES_TSV" > "$SAMPLE_FILE"

N_SAMPLES=$(wc -l < "$SAMPLE_FILE")
echo "=== $N_SAMPLES clinical samples ==="

# =====================================================================
# 1. Clade assignment — nearest-neighbor from full genotype table
# =====================================================================
# Uses ALL samples (clinical + backbone) to compute distances between
# each clinical isolate and every backbone sample, then assigns the
# clade of the nearest neighbor.

echo "=== Extracting genotype table (all samples) ==="

GT_TABLE="${OUT_DIR}/all_samples_genotypes.tsv"

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk VariantsToTable \
    -R "$REF" \
    -V "$FILTERED_VCF" \
    -F CHROM -F POS \
    -GF GT \
    --select-type-to-include SNP \
    -O "$GT_TABLE"

echo "=== Assigning clades ==="

singularity exec "${SIF_DIR}/vcf2phylip-2.9.sif" python3 \
    "${SCRIPTS_DIR}/assign_clades.py" \
    "$GT_TABLE" \
    "$CLADE_METADATA" \
    "$SAMPLE_FILE" \
    "${OUT_DIR}/clinical"

# =====================================================================
# 2. Summary statistics per clinical sample
# =====================================================================
# Collect: mapped reads, coverage, variant counts from existing outputs

echo "=== Collecting summary statistics ==="

STATS_FILE="${OUT_DIR}/clinical_summary_stats.tsv"
echo -e "sample_id\ttotal_snps_in_vcf\tmean_depth\tclade_assigned" > "$STATS_FILE"

# Get clade assignments we just computed
declare -A CLADES
while IFS=$'\t' read -r sid clade rest; do
    CLADES["$sid"]="$clade"
done < <(tail -n +2 "${OUT_DIR}/clinical_clade_assignments.tsv")

while read -r sample_id; do
    # Count variants for this sample (non-missing, non-ref genotypes)
    sample_snps=$(awk -F'\t' -v s="${sample_id}.GT" '
        NR==1 { for(i=1;i<=NF;i++) if($i==s) col=i; next }
        col && $col != "./." && $col != "0/0" && $col != "0|0" { count++ }
        END { print count+0 }
    ' "$GT_TABLE")

    # Get mean depth from BAM if available
    bam="${RESULTS_DIR}/02_mapped/${sample_id}.bam"
    if [[ -f "$bam" ]]; then
        mean_depth=$(singularity exec "${SIF_DIR}/samtools_1.21.sif" samtools depth -a "$bam" \
            | awk '{sum+=$3; n++} END {if(n>0) printf "%.1f", sum/n; else print "NA"}')
    else
        mean_depth="NA"
    fi

    clade="${CLADES[$sample_id]:-NA}"

    echo -e "${sample_id}\t${sample_snps}\t${mean_depth}\t${clade}" >> "$STATS_FILE"
done < "$SAMPLE_FILE"

echo "=== Summary stats: $STATS_FILE ==="

# =====================================================================
# 3. Regenerate iTOL annotations with clade assignments
# =====================================================================
# Now that we have clade assignments, regenerate the iTOL annotations
# so clinical isolates get their clade colors in the colorstrip.

echo "=== Regenerating iTOL annotations with clade assignments ==="

singularity exec "${SIF_DIR}/vcf2phylip-2.9.sif" python3 \
    "${SCRIPTS_DIR}/generate_itol_annotations.py" \
    "$CLADE_METADATA" \
    "$SAMPLES_TSV" \
    "${RESULTS_DIR}/06_phylogeny/itol" \
    "${OUT_DIR}/clinical_clade_assignments.tsv"

# =====================================================================
# Done
# =====================================================================

echo ""
echo "=== Summary outputs ==="
echo "  Clade assignments: ${OUT_DIR}/clinical_clade_assignments.tsv"
echo "  Summary stats:     ${OUT_DIR}/clinical_summary_stats.tsv"
echo "  iTOL annotations:  ${RESULTS_DIR}/06_phylogeny/itol/ (updated with clade colors)"
echo "=== Done ==="
