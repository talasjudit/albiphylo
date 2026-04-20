#!/bin/bash
# sample_detection.sh — sample list generation for C. albicans pipeline
# Sourced by runjob.sh. Provides two detection functions:
#   detect_clinical_samples  — reads metadata/samples.tsv, filters by status
#   detect_public_samples    — scans public_data/ for SRR* paired FASTQs

# ── Clinical samples (from samples.tsv) ───────────────────────────────────────
# Output per line: sample_id,R1_path,R2_path
# Includes all rows where illumina_r1 and illumina_r2 are populated.

detect_clinical_samples() {
    local tsv="$1"

    if [[ ! -f "$tsv" ]]; then
        echo "ERROR: samples.tsv not found: $tsv" >&2
        return 1
    fi

    # Use awk for TSV parsing — handles empty fields correctly (bash read does not).
    # tr -d '\r' strips Windows CRLF line endings (Excel-saved TSVs).
    # TSV columns: 1=sample_id 2=status 3=illu_proj 4=illu_sample
    #              5=illumina_r1 6=illumina_r2 7=ont_proj 8=ont_sample 9=ont_reads 10=assembly
    # Status: blank/"ready" = process normally; "rerun" = force reprocess (step 3 only)
    local output
    output=$(tr -d '\r' < "$tsv" | awk -F'\t' '
        NR == 1 { next }
        $1 == "" { next }
        $5 == "" || $6 == "" { next }
        { print $1 "," $5 "," $6 }
    ')

    if [[ -z "$output" ]]; then
        echo "ERROR: No samples with Illumina reads found in $tsv" >&2
        return 1
    fi

    echo "$output"
}

# ── Public samples (Ropars et al., public_data/) ──────────────────────────────
# Output per line: sample_id,R1_path,R2_path
# Expects paired files named SRR*_1.fastq.gz / SRR*_2.fastq.gz

detect_public_samples() {
    local data_dir="$1"

    if [[ ! -d "$data_dir" ]]; then
        echo "ERROR: Public data directory not found: $data_dir" >&2
        return 1
    fi

    declare -A r1_map r2_map

    while IFS= read -r -d '' filepath; do
        local filename
        filename=$(basename "$filepath")

        if [[ "$filename" =~ ^(SRR[0-9]+)_1\.(fastq|fq)\.gz$ ]]; then
            r1_map["${BASH_REMATCH[1]}"]="$filepath"
        elif [[ "$filename" =~ ^(SRR[0-9]+)_2\.(fastq|fq)\.gz$ ]]; then
            r2_map["${BASH_REMATCH[1]}"]="$filepath"
        fi
    done < <(find "$data_dir" -maxdepth 1 -name "SRR*.fastq.gz" -print0 | sort -z)

    local found=0
    for sample in $(printf '%s\n' "${!r1_map[@]}" | sort); do
        if [[ -n "${r1_map[$sample]}" && -n "${r2_map[$sample]:-}" ]]; then
            echo "${sample},${r1_map[$sample]},${r2_map[$sample]}"
            found=$(( found + 1 ))
        else
            echo "WARNING: $sample missing R2 — skipping" >&2
        fi
    done

    if [[ $found -eq 0 ]]; then
        echo "ERROR: No paired SRR FASTQ files found in $data_dir" >&2
        return 1
    fi
}
