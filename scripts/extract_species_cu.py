#!/usr/bin/env python3
"""
Extracts species names and their total root-to-tip branch lengths in coalescent units (cu)
from a species tree (default: store/63K.tre) and saves the results to a CSV file.
"""

import sys
import argparse
import csv
from pathlib import Path
import treeswift as ts

SCRIPTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPTS_DIR.parent

def calculate_cu_branch_lengths(tree_path: Path) -> list[tuple[str, float]]:
    """Calculates total root-to-tip coalescent unit (cu) distance for each leaf species."""
    tree = ts.read_tree_newick(tree_path)
    results = []
    for leaf in tree.traverse_leaves():
        dist = 0.0
        curr = leaf
        while curr:
            if curr.edge_length is not None:
                dist += curr.edge_length
            curr = curr.get_parent()
        results.append((leaf.label, dist))
    return results

def main():
    parser = argparse.ArgumentParser(description="Extract species and CU branch lengths from 63K.tre")
    parser.add_argument(
        "-i", "--input",
        type=Path,
        default=REPO_ROOT / "store" / "63K.tre",
        help="Input newick tree path (default: store/63K.tre)"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=REPO_ROOT / "species_cu.csv",
        help="Output CSV path (default: species_cu.csv)"
    )
    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Error: Input tree file not found at '{args.input}'")

    print(f"Reading tree from {args.input}...")
    species_cu = calculate_cu_branch_lengths(args.input)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "cu_branch_length"])
        for species, cu in species_cu:
            writer.writerow([species, cu])

    print(f"Successfully wrote {len(species_cu)} species CU branch lengths to '{args.output}'")

if __name__ == "__main__":
    main()
