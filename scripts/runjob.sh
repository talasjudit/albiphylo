#!/bin/bash
# runjob.sh — submit C. albicans pipeline SLURM jobs
# Usage: runjob.sh <script.slurm> [--cohort clinical|public]
#
# --cohort clinical  (default) reads samples from metadata/samples.tsv
# --cohort public    reads Ropars et al. paired FASTQs from public_data/
#
# Per-sample array scripts (steps 1-3) run once per sample.
# Step 3 (HaplotypeCaller) skips samples with existing GVCFs unless status=rerun.
# Single-job scripts (steps 4+) are submitted as-is.

set -euo pipefail

# ===== ARGUMENTS =====

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <script.slurm> [--cohort clinical|public]"
    exit 1
fi

script_file="$1"
script_name=$(basename "$script_file")
cohort="clinical"

# Parse optional --cohort flag
shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cohort) cohort="$2"; shift 2 ;;
        *) echo "ERROR: Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ "$cohort" != "clinical" && "$cohort" != "public" ]]; then
    echo "ERROR: --cohort must be 'clinical' or 'public'"
    exit 1
fi

# ===== VALIDATION =====

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -f "$script_file" ]]; then
    echo "ERROR: Script file '$script_file' not found"
    exit 1
fi

if [[ ! -f "$SCRIPTS_DIR/config.conf" ]]; then
    echo "ERROR: config.conf not found in $SCRIPTS_DIR"
    exit 1
fi

if [[ ! -f "$SCRIPTS_DIR/sample_detection.sh" ]]; then
    echo "ERROR: sample_detection.sh not found in $SCRIPTS_DIR"
    exit 1
fi

source "$SCRIPTS_DIR/config.conf"
source "$SCRIPTS_DIR/sample_detection.sh"

# ===== LOGGING =====

log_dir="${RESULTS_DIR}/logs"
mkdir -p "$log_dir"

email_params=""
if [[ "${email_notifications:-false}" == "true" ]] && [[ -n "${email_address:-}" ]]; then
    email_params="--mail-type=${email_events:-END,FAIL} --mail-user=$email_address"
    echo "Email notifications: $email_address"
fi

# ===== JOB TYPE =====

array_scripts=(
    "01_trim.slurm"
    "02_map.slurm"
    "03_haplotypecaller.slurm"
)

needs_array=false
for s in "${array_scripts[@]}"; do
    [[ "$script_name" == "$s" ]] && needs_array=true && break
done

# ===== ARRAY JOB =====

if [[ "$needs_array" == true ]]; then

    mkdir -p "$log_dir/${script_name%.slurm}"
    sample_list_file="$log_dir/${script_name%.slurm}/${script_name%.slurm}_${cohort}_sample_list.txt"
    rm -f "$sample_list_file"

    echo "Building sample list (cohort: $cohort)..."

    case "$cohort" in
        clinical)
            detect_clinical_samples "$SAMPLES_TSV" > "$sample_list_file"
            ;;
        public)
            detect_public_samples "$PUBLIC_DATA_DIR" > "$sample_list_file"
            ;;
    esac

    # For step 3: remove samples that already have a GVCF (unless status=rerun)
    if [[ "$script_name" == "03_haplotypecaller.slurm" && "$cohort" == "clinical" ]]; then
        filtered_list="${sample_list_file%.txt}_filtered.txt"
        > "$filtered_list"
        while IFS=',' read -r sample_id r1 r2; do
            gvcf="${RESULTS_DIR}/03_gvcfs/${sample_id}.g.vcf.gz"
            # Check if rerun is set in TSV for this sample
            status=$(tr -d '\r' < "$SAMPLES_TSV" | awk -F'\t' -v s="$sample_id" '$1==s {print $2}')
            if [[ -f "$gvcf" && "$status" != "rerun" ]]; then
                echo "  [SKIP] $sample_id — GVCF exists (set status=rerun to force)"
            else
                echo "${sample_id},${r1},${r2}" >> "$filtered_list"
            fi
        done < "$sample_list_file"
        mv "$filtered_list" "$sample_list_file"
    fi

    sample_count=$(wc -l < "$sample_list_file")
    if [[ $sample_count -eq 0 ]]; then
        echo "ERROR: No samples to process"
        exit 1
    fi

    max_index=$((sample_count - 1))
    echo "Submitting array job: $sample_count samples (indices 0-$max_index)"
    echo "Sample list: $sample_list_file"

    mkdir -p "$log_dir/${script_name%.slurm}"

    sbatch \
        --array=0-${max_index} \
        --partition="${SLURM_PARTITION}" \
        --export=SAMPLE_LIST_FILE="$sample_list_file",SCRIPTS_DIR="$SCRIPTS_DIR" \
        --output="$log_dir/${script_name%.slurm}/${script_name%.slurm}-%A_%a.out" \
        --error="$log_dir/${script_name%.slurm}/${script_name%.slurm}-%A_%a.err" \
        $email_params \
        "$script_file"

# ===== SINGLE JOB =====

else
    echo "Submitting single job: $script_name"

    mkdir -p "$log_dir/${script_name%.slurm}"

    sbatch \
        --partition="${SLURM_PARTITION}" \
        --export=SCRIPTS_DIR="$SCRIPTS_DIR" \
        --output="$log_dir/${script_name%.slurm}/${script_name%.slurm}-%j.out" \
        --error="$log_dir/${script_name%.slurm}/${script_name%.slurm}-%j.err" \
        $email_params \
        "$script_file"
fi

echo "Job submitted successfully!"
