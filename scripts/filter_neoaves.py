#!/usr/bin/env python3
"""
Filter FASTA files to retain only Neoaves species by removing
the 37 Palaeognathae and Galloanserae outgroup species.
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


def read_fasta(filepath: Path):
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


def write_fasta(sequences: dict, filepath: Path):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA files to retain only Neoaves species."
    )
    # Support positional input arguments (any number)
    parser.add_argument(
        "positional_inputs",
        type=Path,
        nargs="*",
        help="Input FASTA files to filter.",
    )
    # Support -i / --input / --input-fasta flag
    parser.add_argument(
        "-i",
        "--input",
        "--input-fasta",
        type=Path,
        dest="flag_input",
        action="append",
        help="Input FASTA file to filter.",
    )
    # Support -o / --output / --output-fasta flag
    parser.add_argument(
        "-o",
        "--output",
        "--output-fasta",
        type=Path,
        dest="output_fasta",
        default=None,
        help="Output FASTA file path (only used when a single input file is specified).",
    )

    args = parser.parse_args()

    inputs = []
    if args.flag_input:
        inputs.extend(args.flag_input)
    if args.positional_inputs:
        inputs.extend(args.positional_inputs)

    if not inputs:
        parser.print_help()
        sys.exit("\nError: Must specify at least one input file.")

    for input_path in inputs:
        if not input_path.exists():
            print(f"Error: Input file does not exist: {input_path}", file=sys.stderr)
            continue

        if len(inputs) == 1 and args.output_fasta is not None:
            output_path = args.output_fasta
        else:
            stem = input_path.stem
            if stem.startswith("default"):
                new_stem = stem.replace("default", "default-neoaves", 1)
            else:
                new_stem = f"{stem}-neoaves"
            output_path = input_path.parent / f"{new_stem}{input_path.suffix}"

        print(f"Reading FASTA: {input_path}", file=sys.stderr)
        fasta_data = read_fasta(input_path)
        initial_count = len(fasta_data)

        print("Filtering out 37 Palaeognathae and Galloanserae outgroup species...", file=sys.stderr)
        filtered_data = {h: seq for h, seq in fasta_data.items() if h not in OUTGROUP_TAXA}

        write_fasta(filtered_data, output_path)
        print(f"Filter complete. Retained {len(filtered_data)} / {initial_count} species.", file=sys.stderr)
        print(f"Saved filtered FASTA to: {output_path}\n", file=sys.stderr)


if __name__ == "__main__":
    main()
