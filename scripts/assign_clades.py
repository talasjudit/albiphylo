#!/usr/bin/env python3
"""
assign_clades.py — Assign clades to clinical isolates by nearest-neighbor distance.

For each clinical isolate, computes Anderson et al. 2023 pairwise distances to all
backbone samples, then assigns the clade of the nearest backbone neighbor.

Reports:
  - Assigned clade (clade of nearest backbone sample)
  - Distance to nearest neighbor
  - Top 3 nearest backbone samples (for confidence assessment)
  - Consensus: whether the top 3 neighbors all agree on clade

Input:
  - Genotype table (TSV from GATK VariantsToTable, all samples)
  - Clade metadata (public_data/ropars2018_clade_metadata.tsv)
  - Clinical sample IDs (one per line)

Output:
  - <prefix>_clade_assignments.tsv

Usage:
  python3 assign_clades.py <genotype_table.tsv> <clade_metadata.tsv> <clinical_samples.txt> <output_prefix>
"""

import sys
import csv
from collections import defaultdict


def parse_dosage(gt_str):
    """Convert GATK genotype string to dosage (0, 1, 2) or None for missing."""
    if gt_str in ('./.', '.|.', '.', ''):
        return None
    alleles = gt_str.replace('|', '/').split('/')
    if len(alleles) != 2:
        return None
    try:
        return int(alleles[0]) + int(alleles[1])
    except ValueError:
        return None


def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <genotype_table.tsv> <clade_metadata.tsv> "
              f"<clinical_samples.txt> <output_prefix>", file=sys.stderr)
        sys.exit(1)

    gt_file = sys.argv[1]
    clade_file = sys.argv[2]
    clinical_file = sys.argv[3]
    prefix = sys.argv[4]

    # Load clade metadata (SRR -> clade)
    srr_clade = {}
    with open(clade_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            srr_clade[row['srr_accession']] = row['clade_mlst']
    print(f"Loaded {len(srr_clade)} backbone clade assignments", file=sys.stderr)

    # Load clinical sample IDs
    with open(clinical_file) as f:
        clinical_ids = set(line.strip() for line in f if line.strip())
    print(f"Clinical samples: {len(clinical_ids)}", file=sys.stderr)

    # Read genotype table header to identify samples
    with open(gt_file) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)

    sample_cols = [(i, col.replace('.GT', ''))
                   for i, col in enumerate(header) if col.endswith('.GT')]
    all_samples = [name for _, name in sample_cols]

    # Classify columns as clinical or backbone
    clinical_indices = []
    backbone_indices = []
    for idx, (col_idx, name) in enumerate(sample_cols):
        if name in clinical_ids:
            clinical_indices.append(idx)
        elif name in srr_clade:
            backbone_indices.append(idx)

    print(f"Found in genotype table: {len(clinical_indices)} clinical, "
          f"{len(backbone_indices)} backbone", file=sys.stderr)

    if not clinical_indices or not backbone_indices:
        print("ERROR: Need both clinical and backbone samples in genotype table",
              file=sys.stderr)
        sys.exit(1)

    # Initialise pairwise distance accumulators (clinical x backbone)
    n_clinical = len(clinical_indices)
    n_backbone = len(backbone_indices)
    pair_diff = [[0.0] * n_backbone for _ in range(n_clinical)]
    pair_count = [[0] * n_backbone for _ in range(n_clinical)]

    # Process variants
    n_variants = 0
    with open(gt_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip header
        for row in reader:
            n_variants += 1
            # Parse dosages for all relevant samples
            dosages = []
            for col_idx, _ in sample_cols:
                dosages.append(parse_dosage(row[col_idx]))

            for ci, c_idx in enumerate(clinical_indices):
                if dosages[c_idx] is None:
                    continue
                for bi, b_idx in enumerate(backbone_indices):
                    if dosages[b_idx] is None:
                        continue
                    diff = abs(dosages[c_idx] - dosages[b_idx]) / 2.0
                    pair_diff[ci][bi] += diff
                    pair_count[ci][bi] += 1

    print(f"Variants processed: {n_variants}", file=sys.stderr)

    # Compute normalised distances and find nearest neighbors
    results = []
    for ci, c_idx in enumerate(clinical_indices):
        sample_name = all_samples[c_idx]
        distances = []
        for bi, b_idx in enumerate(backbone_indices):
            b_name = all_samples[b_idx]
            n = pair_count[ci][bi]
            if n > 0:
                dist = pair_diff[ci][bi] / n
                distances.append((dist, b_name, srr_clade.get(b_name, 'NA'), n))

        # Sort by distance
        distances.sort(key=lambda x: x[0])

        if distances:
            top3 = distances[:3]
            assigned_clade = top3[0][2]
            nearest_dist = top3[0][0]
            nearest_sample = top3[0][1]
            nearest_sites = top3[0][3]

            # Check consensus among top 3
            top3_clades = [t[2] for t in top3]
            consensus = "yes" if len(set(top3_clades)) == 1 else "no"

            top3_str = "; ".join(
                f"{t[1]} (clade {t[2]}, dist={t[0]:.6f})" for t in top3
            )
        else:
            assigned_clade = 'NA'
            nearest_dist = 'NA'
            nearest_sample = 'NA'
            nearest_sites = 0
            consensus = 'NA'
            top3_str = 'NA'

        results.append({
            'sample_id': sample_name,
            'assigned_clade': assigned_clade,
            'nearest_backbone': nearest_sample,
            'distance': nearest_dist,
            'sites_compared': nearest_sites,
            'top3_consensus': consensus,
            'top3_neighbors': top3_str,
        })

    # Write output
    out_file = f"{prefix}_clade_assignments.tsv"
    with open(out_file, 'w') as f:
        fields = ['sample_id', 'assigned_clade', 'nearest_backbone', 'distance',
                  'sites_compared', 'top3_consensus', 'top3_neighbors']
        f.write('\t'.join(fields) + '\n')
        for r in results:
            f.write('\t'.join(str(r[field]) for field in fields) + '\n')

    # Print summary
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Clade assignments:", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    for r in results:
        consensus_flag = " *" if r['top3_consensus'] == 'no' else ""
        print(f"  {r['sample_id']:20s}  →  clade {r['assigned_clade']}"
              f"  (nearest: {r['nearest_backbone']}, "
              f"dist={r['distance']:.6f}){consensus_flag}",
              file=sys.stderr)
    print(f"\n  * = top 3 neighbors disagree on clade (check tree manually)",
          file=sys.stderr)
    print(f"\nOutput: {out_file}", file=sys.stderr)


if __name__ == '__main__':
    main()
