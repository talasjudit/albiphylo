#!/bin/bash
#SBATCH --job-name=mlst_typing
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-00:30
#SBATCH --partition=qib-compute

# MLST typing of clinical hybrid assemblies using PubMLST C. albicans scheme.
# Single job. Compares MLST sequence types with WGS nearest-neighbor clade assignments.
# Input:  results/hybrid_assembly/<sample>/<sample>.fasta
# Output: results/mlst/mlst_results.tsv

# ===== SETUP =====

set -euo pipefail

source "${SCRIPTS_DIR}/config.conf"

OUT_DIR="${RESULTS_DIR}/mlst"
mkdir -p "$OUT_DIR"

MLST_DATADIR="${PROJECT_DIR}/databases/pubmlst"
MLST_SCHEME="${MLST_DATADIR}/calbicans"
BLAST_DIR="${MLST_DATADIR}/blast"
BLAST_DB="${BLAST_DIR}/mlst.fa"
LOCI=(AAT1a ACC1 ADP1 MPIb SYA1 VPS13 ZWF1b)

# ===== Build BLAST database if needed =====
# Built here (not in setup_mlst_db.sh) because singularity is only
# available on compute nodes, not login/software nodes.

if [[ ! -f "${BLAST_DB}.nin" ]]; then
    echo "=== Building MLST BLAST database ==="
    mkdir -p "$BLAST_DIR"
    > "$BLAST_DB"
    for locus in "${LOCI[@]}"; do
        awk -v locus="$locus" -v scheme="calbicans" '
            /^>/ {
                id = substr($0, 2)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", id)
                n = split(id, parts, "_")
                allele_num = parts[n]
                print ">" scheme "." locus "-" allele_num
                next
            }
            { print }
        ' "${MLST_SCHEME}/${locus}.tfa" >> "$BLAST_DB"
    done

    singularity exec "${SIF_DIR}/mlst_2.33.1.sif" makeblastdb \
        -in "$BLAST_DB" -dbtype nucl -title "C. albicans MLST" \
        -out "$BLAST_DB" > /dev/null

    echo "[DONE] BLAST DB: $(grep -c '^>' "$BLAST_DB") sequences"
fi

# ===== Collect assemblies from samplesheet =====

ASSEMBLIES=()
while IFS=$'\t' read -r sample_id assembly; do
    [[ -z "$sample_id" ]] && continue
    if [[ -z "$assembly" ]]; then
        echo "WARN: no assembly path for ${sample_id}, skipping" >&2
        continue
    fi
    if [[ ! -f "$assembly" ]]; then
        echo "WARN: assembly not found for ${sample_id}: ${assembly}" >&2
        continue
    fi
    ASSEMBLIES+=("$assembly")
done < <(awk -F'\t' 'NR>1 && $1!="" { print $1"\t"$10 }' "$SAMPLES_TSV")

if [[ ${#ASSEMBLIES[@]} -eq 0 ]]; then
    echo "ERROR: no assemblies found — populate the 'assembly' column in $SAMPLES_TSV" >&2
    exit 1
fi

echo "=== MLST typing: ${#ASSEMBLIES[@]} assemblies ==="

# ===== Run mlst =====

singularity exec "${SIF_DIR}/mlst_2.33.1.sif" mlst \
    --scheme calbicans \
    --datadir "$MLST_DATADIR" \
    --blastdb "${MLST_DATADIR}/blast/mlst.fa" \
    "${ASSEMBLIES[@]}" \
    > "${OUT_DIR}/mlst_results.tsv"

# ===== Display results =====

echo ""
column -t "${OUT_DIR}/mlst_results.tsv"

# ===== Assign MLST clades and compare with WGS =====
# Uses Alkhars et al. 2024 DST→clade mapping (3,215 DSTs, 21 clades)
# to independently assign MLST clades by nearest-DST allelic distance.

WGS_CLADES="${RESULTS_DIR}/09_summary/clinical_clade_assignments.tsv"
DST_CLADE="${MLST_DATADIR}/dst_to_clade.tsv"
PROFILES="${MLST_SCHEME}/calbicans.txt"
COMPARISON="${OUT_DIR}/mlst_vs_wgs_clades.tsv"

if [[ ! -f "$DST_CLADE" ]]; then
    echo ""
    echo "WARN: DST→clade lookup not found at $DST_CLADE"
    echo "      Download Alkhars et al. 2024 S4 Fig and extract with assign_mlst_clades.py"
elif [[ ! -f "$WGS_CLADES" ]]; then
    echo ""
    echo "WARN: WGS clade assignments not found at $WGS_CLADES — run step 9 first."
else
    echo ""
    echo "=== MLST vs WGS clade comparison (Alkhars et al. 2024 DST→clade lookup) ==="

    python3 "${SCRIPTS_DIR}/assign_mlst_clades.py" \
        "$PROFILES" "$DST_CLADE" "${OUT_DIR}/mlst_results.tsv" "$WGS_CLADES" "$COMPARISON"

    echo ""
    column -t "$COMPARISON"
    echo ""
    echo "=== Comparison: ${COMPARISON} ==="
fi

echo ""
echo "=== Done ==="
