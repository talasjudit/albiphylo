#!/bin/bash
#SBATCH --job-name=07_raxml
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00
#SBATCH --partition=qib-compute

# Step 7: Clinical isolate phylogeny — RAxML-NG (GTR+G, bootstrap)
# Single job. Builds ML tree with bootstrap support for clinical isolates only.
# Provides higher resolution than the global FastTree (step 6) for strain comparison.
# Following Adamu Bukari et al. 2025: RAxML for within-group trees.
# Input:  results/05_filtered/all_samples.filtered.vcf + metadata/samples.tsv
# Output: results/07_raxml/clinical.raxml.bestTree
#         results/07_raxml/clinical.raxml.support (tree with bootstrap values)

# ===== SETUP =====

set -euo pipefail

source "${SCRIPTS_DIR}/config.conf"

OUT_DIR="${RESULTS_DIR}/07_raxml"
mkdir -p "$OUT_DIR"

FILTERED_VCF="${RESULTS_DIR}/05_filtered/all_samples.filtered.vcf"

# ===== Extract clinical sample IDs from samples.tsv =====
# Only samples with status "ready" or "rerun" are included.
# As new isolates are added to samples.tsv, re-running this step includes them automatically.

SAMPLE_FILE="${OUT_DIR}/clinical_samples.txt"
awk -F'\t' 'NR>1 && $1!="" { print $1 }' "$SAMPLES_TSV" > "$SAMPLE_FILE"

N_SAMPLES=$(wc -l < "$SAMPLE_FILE")
echo "=== $N_SAMPLES clinical samples ==="
cat "$SAMPLE_FILE"

if [[ "$N_SAMPLES" -lt 4 ]]; then
    echo "ERROR: Need at least 4 samples for meaningful RAxML tree (have $N_SAMPLES)" >&2
    exit 1
fi

# ===== Subset VCF to clinical samples =====
# --remove-unused-alternates: drop alt alleles absent in clinical subset
# --exclude-non-variants: drop sites monomorphic after subsetting

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
echo "=== $N_SNPS biallelic SNPs in clinical samples ==="

# ===== VCF → FASTA alignment (vcf2phylip) =====
# IUPAC ambiguity codes for hets (e.g. R=A/G) — preserves diploid signal.

echo "=== vcf2phylip: VCF -> FASTA alignment ==="

singularity exec "${SIF_DIR}/vcf2phylip-2.9.sif" vcf2phylip \
    --input "$CLINICAL_VCF" \
    --fasta \
    --min-samples-locus 4 \
    --output-prefix "${OUT_DIR}/clinical"

ALIGNMENT="${OUT_DIR}/clinical.min4.fasta"

N_SEQS=$(grep -c '^>' "$ALIGNMENT" || true)
echo "=== Alignment: $N_SEQS sequences ==="

if [[ "$N_SEQS" -lt 4 ]]; then
    echo "ERROR: Fewer than 4 sequences in alignment" >&2
    exit 1
fi

# ===== RAxML-NG =====
# --all: combined ML search + bootstrap convergence test
# GTR+G: general time-reversible + gamma rate heterogeneity (Adamu Bukari 2025)
# 100 bootstrap replicates (sufficient for small sample sets;
#   Adamu Bukari used 1000 for publication-quality trees — increase if needed)

echo "=== RAxML-NG: ML tree + 1000 bootstrap replicates ==="

singularity exec "${SIF_DIR}/raxml-ng-2.0.0--h2105a86_1.sif" raxml-ng \
    --all \
    --msa "$ALIGNMENT" \
    --model GTR+G \
    --bs-trees 1000 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --seed 12345 \
    --prefix "${OUT_DIR}/clinical"

echo "=== Done ==="
echo "Best ML tree:       ${OUT_DIR}/clinical.raxml.bestTree"
echo "Bootstrap support:  ${OUT_DIR}/clinical.raxml.support"
echo "All bootstrap trees: ${OUT_DIR}/clinical.raxml.bootstraps"
