#!/usr/bin/env python3
"""
assign_clades.py — Assign clades to clinical isolates by nearest-neighbor distance.

For each clinical isolate, computes Anderson et al. 2023 pairwise distances to all
backbone samples, then assigns the clade by **majority vote over the top 5 nearest
backbone neighbours**. Backbone samples with < 1000 shared sites are excluded as
unreliable (e.g. SRR6669970, SRR6669899 each have ~200 calls / 9000 SNPs after filtering).

Reports:
  - assigned_clade: majority vote over top 5 neighbours (tie-break: smallest mean distance)
  - nearest_clade: clade of the single closest neighbour (for comparison)
  - nearest_matches_majority: "no" if these differ → boundary sample, inspect on tree
  - top5_consensus: "yes" if all top 5 neighbours agree, else "no"
  - top5_votes: clade vote counts, e.g. "1:3, 9:2"
  - top5_neighbors: per-neighbour distance breakdown

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
    dosage = 0
    for a in alleles:
        if a == alt:
            dosage += 1
        elif a != ref:
            return None
    return dosage


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

    try:
        ref_idx = header.index('REF')
        alt_idx = header.index('ALT')
    except ValueError:
        print("ERROR: genotype table must have REF and ALT columns "
              "(run VariantsToTable with -F REF -F ALT)", file=sys.stderr)
        sys.exit(1)

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
            ref = row[ref_idx]
            alt = row[alt_idx]
            # Parse dosages for all relevant samples
            dosages = []
            for col_idx, _ in sample_cols:
                dosages.append(parse_dosage(row[col_idx], ref, alt))

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

        # Filter out backbone comparisons with too few shared sites — broken
        # samples (e.g. SRR6669970, SRR6669899 with ~200/9000 calls) produce
        # unreliable pairwise distances that inflate their apparent closeness.
        MIN_SITES = 1000
        distances = [d for d in distances if d[3] >= MIN_SITES]

        # Sort by distance
        distances.sort(key=lambda x: x[0])

        if distances:
            top5 = distances[:5]
            # Majority vote over top-5 (tie-break: smallest mean distance)
            from collections import Counter
            clade_votes = Counter(t[2] for t in top5)
            max_votes = max(clade_votes.values())
            tied_clades = [c for c, v in clade_votes.items() if v == max_votes]
            if len(tied_clades) == 1:
                assigned_clade = tied_clades[0]
            else:
                # Tie-break: clade whose top-5 representatives have smallest mean distance
                best = min(
                    tied_clades,
                    key=lambda c: sum(t[0] for t in top5 if t[2] == c)
                                  / sum(1 for t in top5 if t[2] == c))
                assigned_clade = best

            nearest_dist = top5[0][0]
            nearest_sample = top5[0][1]
            nearest_sites = top5[0][3]
            nearest_clade = top5[0][2]

            top5_clades = [t[2] for t in top5]
            consensus = "yes" if len(set(top5_clades)) == 1 else "no"

            # Note mismatch between nearest-neighbour and majority-vote assignment
            method_agree = "yes" if nearest_clade == assigned_clade else "no"

            top5_str = "; ".join(
                f"{t[1]} (clade {t[2]}, dist={t[0]:.6f})" for t in top5
            )
            vote_str = ", ".join(f"{c}:{v}" for c, v in clade_votes.most_common())
        else:
            assigned_clade = 'NA'
            nearest_dist = 'NA'
            nearest_sample = 'NA'
            nearest_sites = 0
            nearest_clade = 'NA'
            consensus = 'NA'
            method_agree = 'NA'
            top5_str = 'NA'
            vote_str = 'NA'

        results.append({
            'sample_id': sample_name,
            'assigned_clade': assigned_clade,
            'nearest_backbone': nearest_sample,
            'nearest_clade': nearest_clade,
            'distance': nearest_dist,
            'sites_compared': nearest_sites,
            'top5_consensus': consensus,
            'nearest_matches_majority': method_agree,
            'top5_votes': vote_str,
            'top5_neighbors': top5_str,
        })

    # Write output
    out_file = f"{prefix}_clade_assignments.tsv"
    with open(out_file, 'w') as f:
        fields = ['sample_id', 'assigned_clade', 'nearest_backbone', 'nearest_clade',
                  'distance', 'sites_compared', 'top5_consensus',
                  'nearest_matches_majority', 'top5_votes', 'top5_neighbors']
        f.write('\t'.join(fields) + '\n')
        for r in results:
            f.write('\t'.join(str(r[field]) for field in fields) + '\n')

    # Print summary
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Clade assignments (top-5 majority vote, min 1000 shared sites):", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    for r in results:
        flag = ""
        if r['top5_consensus'] == 'no':
            flag += " *"
        if r['nearest_matches_majority'] == 'no':
            flag += " BOUNDARY"
        print(f"  {r['sample_id']:20s}  →  clade {r['assigned_clade']}"
              f"  (votes: {r['top5_votes']}){flag}",
              file=sys.stderr)
    print(f"\n  * = top 5 neighbors span multiple clades",
          file=sys.stderr)
    print(f"  BOUNDARY = majority vote differs from single nearest neighbour",
          file=sys.stderr)
    print(f"\nOutput: {out_file}", file=sys.stderr)


if __name__ == '__main__':
    main()
