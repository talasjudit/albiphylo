#!/bin/bash
#SBATCH --job-name=06_fasttree
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00
#SBATCH --partition=qib-compute

# Step 6: Global phylogeny — FastTree (GTR+gamma)
# Single job. Builds global ML tree from SNP alignment of all samples
# (clinical isolates + Ropars et al. 2018 backbone).
# Input:  results/05_filtered/all_samples.filtered.vcf
# Output: results/06_phylogeny/global_tree.nwk

# ===== SETUP =====

set -euo pipefail

source "${SCRIPTS_DIR}/config.conf"

OUT_DIR="${RESULTS_DIR}/06_phylogeny"
mkdir -p "$OUT_DIR"

FILTERED_VCF="${RESULTS_DIR}/05_filtered/all_samples.filtered.vcf"
ALIGNMENT="${OUT_DIR}/all_samples_snps.fasta"
TREE="${OUT_DIR}/global_tree.nwk"

# ===== VCF → SNP-only VCF =====
# Extract biallelic SNPs only — indels and multi-allelic sites are excluded.

SNP_VCF="${OUT_DIR}/all_samples_snps.vcf"

echo "=== Extracting biallelic SNPs ==="

singularity exec "${SIF_DIR}/gatk-4.6.2.0.sif" gatk SelectVariants \
    -R "$REF" \
    -V "$FILTERED_VCF" \
    -O "$SNP_VCF" \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC

echo "=== $(grep -cv '^#' "$SNP_VCF") biallelic SNPs ==="

# ===== VCF → FASTA alignment (vcf2phylip) =====
# vcf2phylip (Ortiz 2019, Mol Ecol Resour) converts multi-sample diploid VCF
# to FASTA alignment. Hets encoded as IUPAC ambiguity codes by default
# (e.g. R=A/G, Y=C/T) — preserves het signal for C. albicans phylogeny.
# snp-sites was originally chosen but cannot parse GATK's diploid VCF format.

echo "=== vcf2phylip: VCF -> FASTA alignment ==="

singularity exec "${SIF_DIR}/vcf2phylip-2.9.sif" vcf2phylip \
    --input "$SNP_VCF" \
    --fasta \
    --min-samples-locus 4 \
    --output-prefix "${OUT_DIR}/all_samples"

# vcf2phylip names output as <prefix>.min<N>.fasta
mv "${OUT_DIR}/all_samples.min4.fasta" "$ALIGNMENT"

N_SEQS=$(grep -c '^>' "$ALIGNMENT" || true)
echo "=== Alignment: $ALIGNMENT ($N_SEQS sequences) ==="

if [[ "$N_SEQS" -eq 0 ]]; then
    echo "ERROR: No sequences in alignment — check filtered VCF" >&2
    exit 1
fi

# ===== FastTree =====

echo "=== FastTree global tree ==="

singularity exec "${SIF_DIR}/fasttree_2.2.0.sif" FastTree \
    -gtr \
    -gamma \
    -nt \
    "$ALIGNMENT" \
    > "$TREE"

echo "=== Done: $TREE ==="
