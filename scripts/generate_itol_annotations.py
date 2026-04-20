#!/usr/bin/env python3
"""
generate_itol_annotations.py — Create iTOL annotation files for the global tree.

Reads:
  - Ropars et al. 2018 clade metadata (public_data/ropars2018_clade_metadata.tsv)
  - Clinical sample list (metadata/samples.tsv)
  - Optionally: clade assignments for clinical isolates (from assign_clades.py)

Produces iTOL annotation files in the output directory:
  1. itol_clade_colorstrip.txt  — color strip: clade for all samples
  2. itol_clinical_labels.txt   — bold red labels for clinical isolates
  3. itol_popup_info.txt        — hover popups with strain details

Usage:
  python3 generate_itol_annotations.py <clade_metadata.tsv> <samples.tsv> <output_dir> [clade_assignments.tsv]

  The optional 4th argument is the output from assign_clades.py (step 9).
  If provided, clinical isolates get clade colors in the colorstrip.
  If omitted, clinical isolates are left blank (clade unknown).

Upload the global tree (.nwk) to iTOL, then drag-and-drop annotation files onto the tree.
To root: find the C. africana (black) + clade 13 (magenta) cluster → re-root on that branch.
"""

import sys
import os
import csv

# Clade color palette — visually distinct colors for each clade
CLADE_COLORS = {
    1:  '#e41a1c',   # red
    2:  '#377eb8',   # blue
    3:  '#4daf4a',   # green
    4:  '#984ea3',   # purple
    5:  '#ff7f00',   # orange
    7:  '#a65628',   # brown
    8:  '#f781bf',   # pink
    9:  '#999999',   # grey
    10: '#66c2a5',   # teal
    11: '#fc8d62',   # salmon
    12: '#8da0cb',   # light blue
    13: '#e78ac3',   # magenta
    16: '#a6d854',   # lime
    18: '#ffd92f',   # yellow
    'NC':          '#bdbdbd',  # light grey — not classified
    'C. africana': '#000000',  # black — outgroup
}
DEFAULT_COLOR = '#636363'


def get_color(clade):
    """Get color for a clade, handling string/int conversion."""
    if clade in CLADE_COLORS:
        return CLADE_COLORS[clade]
    try:
        return CLADE_COLORS.get(int(clade), DEFAULT_COLOR)
    except (ValueError, TypeError):
        return DEFAULT_COLOR


def main():
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print(f"Usage: {sys.argv[0]} <clade_metadata.tsv> <samples.tsv> <output_dir> "
              f"[clade_assignments.tsv]", file=sys.stderr)
        sys.exit(1)

    clade_file = sys.argv[1]
    samples_file = sys.argv[2]
    out_dir = sys.argv[3]
    assignments_file = sys.argv[4] if len(sys.argv) == 5 else None
    os.makedirs(out_dir, exist_ok=True)

    # ---- Load clade metadata (SRR -> clade) ----
    srr_info = {}
    with open(clade_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            srr_info[row['srr_accession']] = {
                'strain': row['strain'],
                'clade': row['clade_mlst'],
                'country': row.get('country', ''),
                'site': row.get('site', ''),
            }
    print(f"Loaded {len(srr_info)} backbone samples with clade info", file=sys.stderr)

    # ---- Load clinical samples ----
    clinical_ids = []
    with open(samples_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('sample_id', '').strip():
                clinical_ids.append(row['sample_id'].strip())
    print(f"Loaded {len(clinical_ids)} clinical samples", file=sys.stderr)

    # ---- Load clinical clade assignments (optional) ----
    clinical_clades = {}
    if assignments_file and os.path.exists(assignments_file):
        with open(assignments_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                clinical_clades[row['sample_id']] = row['assigned_clade']
        print(f"Loaded clade assignments for {len(clinical_clades)} clinical samples",
              file=sys.stderr)
    else:
        if assignments_file:
            print(f"WARNING: {assignments_file} not found — clinical clades will be blank",
                  file=sys.stderr)

    # ---- Collect unique clades for legend ----
    all_clades = set(info['clade'] for info in srr_info.values())
    all_clades.update(clinical_clades.values())
    clades_seen = sorted(all_clades,
                         key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))

    # ================================================================
    # 1. Clade color strip
    # ================================================================
    with open(os.path.join(out_dir, 'itol_clade_colorstrip.txt'), 'w') as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tClade (MLST)\n")
        f.write("COLOR\t#666666\n")
        f.write("STRIP_WIDTH\t25\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("\n")

        # Legend
        f.write("LEGEND_TITLE\tClade (MLST)\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(clades_seen)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(get_color(c) for c in clades_seen) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(f"Clade {c}" for c in clades_seen) + "\n")
        f.write("\n")

        f.write("DATA\n")
        # Backbone samples
        for srr, info in sorted(srr_info.items()):
            color = get_color(info['clade'])
            label = f"Clade {info['clade']}"
            f.write(f"{srr}\t{color}\t{label}\n")

        # Clinical samples — with clade color if assignment available, otherwise omitted
        for sample_id in clinical_ids:
            if sample_id in clinical_clades:
                clade = clinical_clades[sample_id]
                color = get_color(clade)
                f.write(f"{sample_id}\t{color}\tClade {clade}\n")

    # ================================================================
    # 2. Clinical isolate labels — bold red to stand out
    # ================================================================
    with open(os.path.join(out_dir, 'itol_clinical_labels.txt'), 'w') as f:
        f.write("DATASET_STYLE\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tClinical isolates\n")
        f.write("\n")
        f.write("DATA\n")
        for sample_id in clinical_ids:
            f.write(f"{sample_id}\tlabel\tnode\t#d62728\t2\tbold\n")

    # ================================================================
    # 3. Popup info (hover text with strain details)
    # ================================================================
    with open(os.path.join(out_dir, 'itol_popup_info.txt'), 'w') as f:
        f.write("POPUP_INFO\n")
        f.write("SEPARATOR TAB\n")
        f.write("\n")
        f.write("DATA\n")
        for srr, info in sorted(srr_info.items()):
            title = f"{srr} ({info['strain']})"
            content = (f"Strain: {info['strain']}<br>"
                       f"Clade (MLST): {info['clade']}<br>"
                       f"Country: {info['country']}<br>"
                       f"Site: {info['site']}")
            f.write(f"{srr}\t{title}\t{content}\n")
        for sample_id in clinical_ids:
            clade = clinical_clades.get(sample_id, 'not yet assigned')
            f.write(f"{sample_id}\t{sample_id} (clinical)\t"
                    f"Clinical isolate — VVC research<br>"
                    f"Assigned clade: {clade}\n")

    print(f"\niTOL annotation files written to {out_dir}/:", file=sys.stderr)
    print(f"  itol_clade_colorstrip.txt  — clade colors"
          f"{' (incl. clinical)' if clinical_clades else ' (backbone only)'}", file=sys.stderr)
    print(f"  itol_clinical_labels.txt   — bold red labels (clinical isolates)", file=sys.stderr)
    print(f"  itol_popup_info.txt        — hover popups with strain details", file=sys.stderr)
    print(f"\nUpload tree to iTOL, then drag-and-drop these files onto it.", file=sys.stderr)
    print(f"To root: find C. africana (black) + clade 13 (magenta) → re-root on that branch.",
          file=sys.stderr)


if __name__ == '__main__':
    main()
