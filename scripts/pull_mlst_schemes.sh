#!/bin/bash
# setup_mlst_db.sh — download C. albicans MLST scheme from PubMLST
#
# C. albicans is not in the standard mlst/mlstdb bacterial databases.
# This script pulls the scheme (7 loci, ~4700 STs) from the PubMLST
# REST API and formats it for Seemann's mlst tool (--datadir option).
#
# Downloads only — the BLAST database is built by mlst_typing.slurm on a
# compute node (singularity is required for makeblastdb but unavailable
# on login/software nodes).
#
# Idempotent — skips downloads if files already exist. Re-run with --force to refresh.
# Requires: curl, internet access (run on login node).
#
# Usage:
#   bash scripts/setup_mlst_db.sh [--force]

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPTS_DIR/config.conf"

MLST_DB="${PROJECT_DIR}/databases/pubmlst/calbicans"
API_BASE="https://rest.pubmlst.org/db/pubmlst_calbicans_seqdef"
SCHEME_ID=1
LOCI=(AAT1a ACC1 ADP1 MPIb SYA1 VPS13 ZWF1b)

FORCE=false
[[ "${1:-}" == "--force" ]] && FORCE=true

echo "=== C. albicans MLST database setup ==="
echo "Scheme: PubMLST C. albicans MLST (${#LOCI[@]} loci)"
echo "Output: $MLST_DB"
echo ""

mkdir -p "$MLST_DB"

# ── Download allele sequences ────────────────────────────────────────────────

for locus in "${LOCI[@]}"; do
    tfa="${MLST_DB}/${locus}.tfa"
    if [[ -f "$tfa" && "$FORCE" == false ]]; then
        echo "[SKIP] ${locus}.tfa exists ($(grep -c '^>' "$tfa") alleles)"
    else
        echo "[DOWN] ${locus} alleles..."
        curl -sf "${API_BASE}/loci/${locus}/alleles_fasta" > "$tfa"
        echo "[DONE] ${locus}: $(grep -c '^>' "$tfa") alleles"
    fi
done

# ── Download profiles (ST definitions) ───────────────────────────────────────

PROFILES="${MLST_DB}/calbicans.txt"
if [[ -f "$PROFILES" && "$FORCE" == false ]]; then
    echo "[SKIP] profiles exist ($(tail -n +2 "$PROFILES" | wc -l) STs)"
else
    echo "[DOWN] ST profiles..."
    curl -sf "${API_BASE}/schemes/${SCHEME_ID}/profiles_csv" > "$PROFILES"
    echo "[DONE] $(tail -n +2 "$PROFILES" | wc -l) STs"
fi

echo ""
echo "=== Downloads complete ==="
echo "  Alleles:  ${MLST_DB}/*.tfa (${#LOCI[@]} loci)"
echo "  Profiles: ${PROFILES}"
echo ""
echo "BLAST database will be built automatically by mlst_typing.slurm."
echo "Submit: bash scripts/runjob.sh mlst_typing.slurm"
