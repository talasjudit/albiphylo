#!/usr/bin/env python3
"""
pairwise_snp_distance.py — Pairwise SNP distance matrix from GATK genotype table.

Computes pairwise distances between all sample pairs following Anderson et al. 2023:
  - Each genotype coded as dosage: 0/0=0, 0/1=1, 1/1=2
  - Per-site distance = |dosage_i - dosage_j| / 2
    D(ref/ref, ref/alt) = 0.5
    D(ref/ref, alt/alt) = 1.0
    D(ref/alt, alt/alt) = 0.5
  - Total distance = sum(per-site distances) / number of non-missing sites

Outputs:
  <prefix>_distance_matrix.tsv  — normalised pairwise distance (Anderson formula)
  <prefix>_snp_counts.tsv       — raw count of sites where genotypes differ

Usage:
  python3 pairwise_snp_distance.py <genotype_table.tsv> <output_prefix>

Input: TSV from GATK VariantsToTable with -GF GT
  Columns: CHROM  POS  sample1.GT  sample2.GT  ...
"""

import sys
import csv
from itertools import combinations


def parse_dosage(gt_str, ref, alt):
    """Nucleotide genotype → ALT dosage (0, 1, 2) or None for missing.

    GATK VariantsToTable with -GF GT emits nucleotide genotypes (e.g. 'C/T'),
    not numeric ('0/1'). Dosage = count of alleles matching ALT.
    """
    if gt_str in ('./.', '.|.', '.', '', 'NA'):
        return None
    alleles = gt_str.replace('|', '/').split('/')
    if len(alleles) != 2:
        return None
    if any(a == '.' for a in alleles):
        return None
    # Biallelic restriction upstream means alleles should be REF or ALT only
    dosage = 0
    for a in alleles:
        if a == alt:
            dosage += 1
        elif a != ref:
            return None  # unexpected allele — treat as missing
    return dosage


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <genotype_table.tsv> <output_prefix>",
              file=sys.stderr)
        sys.exit(1)

    gt_file = sys.argv[1]
    prefix = sys.argv[2]

    # Read header to identify sample columns
    with open(gt_file) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)

    # Locate REF/ALT columns (needed to convert nucleotide GT to dosage)
    try:
        ref_idx = header.index('REF')
        alt_idx = header.index('ALT')
    except ValueError:
        print("ERROR: genotype table must have REF and ALT columns "
              "(run VariantsToTable with -F REF -F ALT)", file=sys.stderr)
        sys.exit(1)

    # Sample columns end with .GT
    sample_cols = [(i, col.replace('.GT', ''))
                   for i, col in enumerate(header) if col.endswith('.GT')]
    samples = [name for _, name in sample_cols]
    n_samples = len(samples)

    print(f"Samples: {n_samples}", file=sys.stderr)
    for s in samples:
        print(f"  {s}", file=sys.stderr)

    # Initialise pairwise counters
    pair_diff = {}    # sum of |dosage_i - dosage_j| / 2
    pair_count = {}   # number of non-missing sites compared
    pair_raw = {}     # count of sites where genotypes differ at all

    for i, j in combinations(range(n_samples), 2):
        pair_diff[(i, j)] = 0.0
        pair_count[(i, j)] = 0
        pair_raw[(i, j)] = 0

    # Process variants
    n_variants = 0
    with open(gt_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip header
        for row in reader:
            n_variants += 1
            ref = row[ref_idx]
            alt = row[alt_idx]
            dosages = []
            for col_idx, _ in sample_cols:
                dosages.append(parse_dosage(row[col_idx], ref, alt))

            for i, j in combinations(range(n_samples), 2):
                if dosages[i] is not None and dosages[j] is not None:
                    diff = abs(dosages[i] - dosages[j]) / 2.0
                    pair_diff[(i, j)] += diff
                    pair_count[(i, j)] += 1
                    if dosages[i] != dosages[j]:
                        pair_raw[(i, j)] += 1

    print(f"Variants processed: {n_variants}", file=sys.stderr)

    # Write normalised distance matrix
    out_dist = f"{prefix}_distance_matrix.tsv"
    with open(out_dist, 'w') as f:
        f.write('\t' + '\t'.join(samples) + '\n')
        for i in range(n_samples):
            row = [samples[i]]
            for j in range(n_samples):
                if i == j:
                    row.append('0.000000')
                else:
                    a, b = (i, j) if i < j else (j, i)
                    n = pair_count[(a, b)]
                    row.append(f"{pair_diff[(a, b)] / n:.6f}" if n > 0 else 'NA')
            f.write('\t'.join(row) + '\n')

    # Write raw SNP count matrix
    out_snps = f"{prefix}_snp_counts.tsv"
    with open(out_snps, 'w') as f:
        f.write('\t' + '\t'.join(samples) + '\n')
        for i in range(n_samples):
            row = [samples[i]]
            for j in range(n_samples):
                if i == j:
                    row.append('0')
                else:
                    a, b = (i, j) if i < j else (j, i)
                    row.append(str(pair_raw[(a, b)]))
            f.write('\t'.join(row) + '\n')

    # Summary
    print(f"\nPairwise summary (raw differing SNPs):", file=sys.stderr)
    for i, j in combinations(range(n_samples), 2):
        n = pair_count[(i, j)]
        print(f"  {samples[i]} vs {samples[j]}: "
              f"{pair_raw[(i, j)]} SNPs differ "
              f"({n} sites compared, "
              f"dist={pair_diff[(i, j)] / n:.6f})" if n > 0 else "NA",
              file=sys.stderr)

    print(f"\nOutput: {out_dist}", file=sys.stderr)
    print(f"Output: {out_snps}", file=sys.stderr)


if __name__ == '__main__':
    main()
