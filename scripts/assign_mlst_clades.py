#!/usr/bin/env python3
"""Assign MLST clades to clinical isolates by nearest-DST lookup.

Uses Alkhars et al. 2024 DST→clade mapping (3,215 DSTs, 21 clades)
and PubMLST allele profiles to find the nearest classified DST for
each isolate's allele profile, then assigns that DST's clade.

Usage:
    python3 assign_mlst_clades.py <profiles> <dst_clade> <mlst_results> <wgs_clades> <output>

Args:
    profiles:     PubMLST profiles file (calbicans.txt)
    dst_clade:    DST→clade lookup (dst_to_clade.tsv)
    mlst_results: mlst tool output (mlst_results.tsv)
    wgs_clades:   Step 9 clade assignments (clinical_clade_assignments.tsv)
    output:       Output comparison file
"""
import sys, re
from collections import defaultdict

profiles_file, dst_clade_file, mlst_file, wgs_file, out_file = sys.argv[1:6]

dst_clade = {}
with open(dst_clade_file) as f:
    f.readline()
    for line in f:
        dst, clade = line.strip().split('\t')
        dst_clade[int(dst)] = clade

sts = {}
with open(profiles_file) as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        sts[int(parts[0])] = tuple(parts[1:8])

labeled = {st: dst_clade[st] for st in sts if st in dst_clade}

wgs_clades = {}
with open(wgs_file) as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        wgs_clades[parts[0]] = parts[1]

with open(out_file, 'w') as fout:
    fout.write("sample_id\tmlst_clade\twgs_clade\tconcordant\tnearest_dst\tdist\ttop5_votes\tmlst_alleles\n")

    for line in open(mlst_file):
        parts = line.strip().split('\t')
        sample = parts[0].split('/')[-1].replace('.fasta', '')

        query = []
        for i in range(3, 10):
            m = re.search(r'\(~?(\d+)\??\)', parts[i])
            query.append(m.group(1) if m else None)
        allele_str = ' '.join(parts[3:10])

        matches = []
        for st, alleles in sts.items():
            if st not in labeled:
                continue
            shared = sum(1 for q, a in zip(query, alleles) if q is not None and q == a)
            matches.append((7 - shared, st, labeled[st]))
        matches.sort(key=lambda x: (x[0], x[1]))

        clade_votes = defaultdict(int)
        for d, st, cl in matches[:5]:
            clade_votes[cl] += 1
        top_votes = sorted(clade_votes.items(), key=lambda x: -x[1])

        mlst_clade = top_votes[0][0]
        wgs = wgs_clades.get(sample, 'NA')
        concordant = "yes" if mlst_clade == wgs else "no"
        votes_str = ', '.join(f'{cl}:{n}' for cl, n in top_votes[:3])

        fout.write(f"{sample}\t{mlst_clade}\t{wgs}\t{concordant}\t"
                   f"{matches[0][1]}\t{matches[0][0]}\t{votes_str}\t{allele_str}\n")

        print(f"{sample:<12s}  MLST={mlst_clade:>3s}  WGS={wgs:>3s}  "
              f"{'OK' if concordant == 'yes' else '**':>2s}  "
              f"DST-{matches[0][1]} ({matches[0][0]}/7)  {votes_str}")
