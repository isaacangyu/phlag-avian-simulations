#!/usr/bin/env python3
"""
Checks if a FASTA file contains outgroup taxa species.
If unfiltered (contains any outgroup taxa), filters out those species and writes a new FASTA file.
If already filtered, reports that no changes are needed.
"""

import sys
import argparse
from pathlib import Path

# Standard 37 Palaeognathae and Galloanserae outgroup taxa
OUTGROUP_TAXA = {
    # Palaeognathae
    "Struthio_camelus", "Rhea_americana", "Rhea_pennata", "Casuarius_casuarius",
    "Dromaius_novaehollandiae", "Apteryx_haastii", "Apteryx_owenii", "Apteryx_australis",
    "Apteryx_rowi", "Eudromia_elegans", "Nothoprocta_perdicaria", "Nothoprocta_ornata",
    "Nothoprocta_pentlandii", "Nothocercus_julius", "Nothocercus_nigrocapillus",
    "Tinamus_guttatus", "Crypturellus_soui", "Crypturellus_cinnamomeus", "Crypturellus_undulatus",
    # Galloanserae
    "Chauna_torquata", "Anseranas_semipalmata", "Anser_cygnoid", "Cairina_moschata",
    "Asarcornis_scutulata", "Anas_platyrhynchos", "Anas_zonorhyncha", "Alectura_lathami",
    "Penelope_pileata", "Numida_meleagris", "Odontophorus_gujanensis", "Colinus_virginianus",
    "Callipepla_squamata", "Gallus_gallus", "Coturnix_japonica", "Phasianus_colchicus",
    "Meleagris_gallopavo", "Tympanuchus_cupido",
}


def read_fasta(filepath: Path) -> dict[str, str]:
    sequences = {}
    current_header = None
    current_seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)
    return sequences


def write_fasta(sequences: dict[str, str], filepath: Path):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n{seq}\n")


def filter_fasta_file(input_path: Path, output_path: Path = None, in_place: bool = False) -> bool:
    """
    Checks if input_path FASTA has outgroup taxa.
    If it contains outgroup taxa, filters them out and writes a FASTA file.
    If in_place is True or output_path == input_path, replaces the input file in-place.
    Returns True if filtered output was created/updated, False if already filtered.
    """
    if not input_path.exists():
        sys.exit(f"Error: Input file '{input_path}' does not exist.")

    fasta_data = read_fasta(input_path)
    found_outgroups = [h for h in fasta_data if h in OUTGROUP_TAXA]

    if not found_outgroups:
        print(f"[2/4 Check Neoaves Filter] '{input_path.name}' is already filtered ({len(fasta_data)} species).", file=sys.stderr)
        return False

    print(f"[2/4 Check Neoaves Filter] '{input_path.name}' is unfiltered (contains {len(found_outgroups)} outgroup taxa). Filtering...", file=sys.stderr)

    if in_place:
        target_output = input_path
    elif output_path is None:
        # Default output name: <stem>_filtered.fa
        target_output = input_path.parent / f"{input_path.stem}_filtered{input_path.suffix or '.fa'}"
    else:
        target_output = output_path

    filtered_data = {h: seq for h, seq in fasta_data.items() if h not in OUTGROUP_TAXA}
    write_fasta(filtered_data, target_output)

    print(f"[2/4 Check Neoaves Filter] Filter complete: retained {len(filtered_data)} / {len(fasta_data)} species in '{target_output}'.", file=sys.stderr)
    return True


def filter_tsv_file(input_path: Path, output_path: Path = None, in_place: bool = False) -> bool:
    """
    Checks if input_path TSV mapping file contains outgroup taxa in the first column.
    If it contains outgroup taxa, filters them out and writes a TSV file.
    If in_place is True or output_path == input_path, replaces the input file in-place.
    Returns True if filtered output was created/updated, False if already filtered.
    """
    if not input_path.exists():
        sys.exit(f"Error: Input file '{input_path}' does not exist.")

    lines = []
    found_outgroups = []
    with open(input_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            taxon = parts[0].strip()
            if taxon in OUTGROUP_TAXA:
                found_outgroups.append(taxon)
            else:
                lines.append(line)

    if not found_outgroups:
        print(f"[2/4 Check Neoaves Filter] TSV mapping '{input_path.name}' is already filtered ({len(lines)} taxa).", file=sys.stderr)
        return False

    print(f"[2/4 Check Neoaves Filter] TSV mapping '{input_path.name}' is unfiltered (contains {len(found_outgroups)} outgroup taxa). Filtering...", file=sys.stderr)

    if in_place:
        target_output = input_path
    elif output_path is None:
        target_output = input_path.parent / f"{input_path.stem}_filtered{input_path.suffix or '.tsv'}"
    else:
        target_output = output_path

    target_output.parent.mkdir(parents=True, exist_ok=True)
    with open(target_output, "w") as f:
        f.writelines(lines)

    print(f"[2/4 Check Neoaves Filter] Filter complete: retained {len(lines)} taxa in '{target_output}'.", file=sys.stderr)
    return True


def filter_file(input_path: Path, output_path: Path = None, in_place: bool = False) -> bool:
    """
    Filters outgroup taxa from a FASTA or TSV mapping file based on file extension.
    """
    if input_path.suffix.lower() == ".tsv":
        return filter_tsv_file(input_path, output_path=output_path, in_place=in_place)
    else:
        return filter_fasta_file(input_path, output_path=output_path, in_place=in_place)


def main():
    parser = argparse.ArgumentParser(description="Check and filter outgroup taxa from FASTA or TSV mapping file(s).")
    parser.add_argument("input_files", type=Path, nargs="+", help="Path to input FASTA or TSV mapping file(s).")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output file path (only valid if single input file specified).")
    parser.add_argument("-i", "--in-place", action="store_true", help="Overwrite input file(s) in-place.")
    args = parser.parse_args()

    if args.output and len(args.input_files) > 1:
        sys.exit("Error: -o/--output can only be specified when a single input file is provided.")

    for f in args.input_files:
        out_p = args.output if len(args.input_files) == 1 else None
        filter_file(f, output_path=out_p, in_place=args.in_place)


if __name__ == "__main__":
    main()
