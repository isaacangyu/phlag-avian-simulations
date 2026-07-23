import os
import sys
import argparse
from pathlib import Path


def get_output_dir() -> Path:
    phlag_dir = os.environ.get("PHLAG_DIR")
    if not phlag_dir:
        current = Path(__file__).resolve().parent
        for _ in range(3):
            env_path = current / ".env"
            if env_path.exists():
                try:
                    with open(env_path, "r") as f:
                        for line in f:
                            line = line.strip()
                            if line and not line.startswith("#") and "=" in line:
                                key, val = line.split("=", 1)
                                if key.strip() == "PHLAG_DIR":
                                    phlag_dir = val.strip().strip("'").strip('"')
                                    break
                except Exception:
                    pass
            if phlag_dir:
                break
            current = current.parent

    if phlag_dir:
        return Path(phlag_dir)
    else:
        repo_root = Path(__file__).resolve().parent.parent
        return repo_root / "data"


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


def mix_sequences(
    null_seq: str, alt_seq: str, pattern: str
) -> str:
    pieces = []
    for char in pattern:
        if char == "n":
            pieces.append(null_seq)
        elif char == "a":
            pieces.append(alt_seq)
        else:
            sys.exit(f"Error: Invalid character '{char}' in pattern. Only 'n' and 'a' (case-insensitive) are allowed.")
            
    return "".join(pieces)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate FASTA sequences by mixing null and alt/anomalous files using a pattern string (e.g. 'nanan')."
    )
    
    # Support positional arguments
    parser.add_argument(
        "positional_null",
        type=Path,
        nargs="?",
        help="Path to null FASTA file.",
    )
    parser.add_argument(
        "positional_alt",
        type=Path,
        nargs="?",
        help="Path to alternative/anomalous FASTA file.",
    )
    parser.add_argument(
        "positional_pattern",
        type=str,
        nargs="?",
        help="Pattern string specifying order of sequences (e.g., 'nanan').",
    )

    # Support flags
    parser.add_argument(
        "-n",
        "--null",
        type=Path,
        help="Path to null FASTA file.",
    )
    parser.add_argument(
        "-a",
        "--alt",
        "--alternative",
        "--anomalous",
        type=Path,
        dest="alt",
        help="Path to alternative/anomalous FASTA file.",
    )
    parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        help="Pattern string specifying order of sequences (e.g., 'nanan').",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Output FASTA path.",
    )

    args = parser.parse_args()

    # Resolve values preferring flags, falling back to positionals
    null_path = args.null or args.positional_null
    alt_path = args.alt or args.positional_alt
    pattern = args.pattern or args.positional_pattern

    if not null_path or not alt_path or not pattern:
        parser.print_help()
        sys.exit(
            "\nError: Must specify null file, alt file, and pattern (either as positional arguments or via flags)."
        )

    pattern = pattern.lower()
    for char in pattern:
        if char not in ("n", "a"):
            sys.exit(f"Error: Invalid character '{char}' in pattern. Only 'n' and 'a' (case-insensitive) are allowed.")

    if not null_path.exists():
        sys.exit(f"Error: Null FASTA file does not exist: {null_path}")
    if not alt_path.exists():
        sys.exit(f"Error: Alt/Anomalous FASTA file does not exist: {alt_path}")

    print(f"Reading null FASTA: {null_path}", file=sys.stderr)
    null_fasta = read_fasta(null_path)

    print(f"Reading alt FASTA: {alt_path}", file=sys.stderr)
    alt_fasta = read_fasta(alt_path)

    # Master list of headers from null fasta
    base_headers = list(null_fasta.keys())

    # Verify header consistency
    message = []
    null_missing_headers = set(base_headers) - set(alt_fasta.keys())
    if null_missing_headers:
        message.append(f"Headers present in null FASTA but missing in alt FASTA: {null_missing_headers}")

    
    alt_missing_headers = set(alt_fasta.keys()) - set(base_headers)
    if alt_missing_headers:
        message.append(f"Headers present in alt FASTA but missing in null FASTA: {alt_missing_headers}")

    if message:
        sys.exit("\n".join(message))

    output_dir = get_output_dir()
    (output_dir / "msa").mkdir(parents=True, exist_ok=True)

    if args.output is None:
        out_name = f"concat_{null_path.stem}_{alt_path.stem}_{pattern}.fa"
        output_path = output_dir / "msa" / out_name
    elif not args.output.is_absolute():
        output_path = output_dir / "msa" / args.output.name
    else:
        output_path = args.output

    mixed_fasta = {}
    print(f"Mixing sequences with pattern '{pattern}'...", file=sys.stderr)
    for header in base_headers:
        null_seq = null_fasta[header]
        alt_seq = alt_fasta[header]
        if len(null_seq) != len(alt_seq):
            sys.exit(
                f"Error: Length mismatch for header '{header}': null sequence length is {len(null_seq)}, alt sequence length is {len(alt_seq)}."
            )
        
        mixed_fasta[header] = mix_sequences(null_seq, alt_seq, pattern)

    print(f"Writing output FASTA: {output_path}", file=sys.stderr)
    with open(output_path, "w") as f:
        for header in base_headers:
            f.write(f">{header}\n")
            f.write(f"{mixed_fasta[header]}\n")

    print("Done.", file=sys.stderr)
