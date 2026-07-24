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
        sys.exit("Error: PHLAG_DIR environment variable is not set and could not be found in .env file.")


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


def resolve_fasta_file(f: Path, output_dir: Path) -> Path:
    if f.suffix == ".ts" or str(f).endswith(".ts"):
        simulated_fa = output_dir / "msa" / "simulated" / f"{f.stem}.fa"
        if simulated_fa.exists():
            return simulated_fa
        elif (f.parent / f"{f.stem}.fa").exists():
            return f.parent / f"{f.stem}.fa"
        elif f.with_suffix(".fa").exists():
            return f.with_suffix(".fa")
        elif f.exists():
            print(f"Simulating sequences for {f} -> {simulated_fa}...", file=sys.stderr)
            try:
                from simulate_sequences import simulate_ts_to_fasta
            except ImportError:
                sys.path.insert(0, str(Path(__file__).resolve().parent))
                from simulate_sequences import simulate_ts_to_fasta
            simulate_ts_to_fasta(f, simulated_fa)
            return simulated_fa
        else:
            sys.exit(f"Error: Tree sequence file does not exist: {f}")
    else:
        if f.exists():
            return f
        alt1 = output_dir / "msa" / "simulated" / f.name
        if alt1.exists():
            return alt1
        alt2 = output_dir / "msa" / "concat" / f.name
        if alt2.exists():
            return alt2
        sys.exit(f"Error: Sequence file does not exist: {f}")


def mix_sequences(null_seq: str, alt_seq: str, pattern: str) -> str:
    pieces = []
    for char in pattern:
        if char == "n":
            pieces.append(null_seq)
        elif char == "a":
            pieces.append(alt_seq)
        else:
            sys.exit(f"Error: Invalid character '{char}' in pattern. Only 'n' and 'a' (case-insensitive) are allowed.")
    return "".join(pieces)


import re


def parse_pattern_to_sequence(pattern: str) -> list[str]:
    pattern = pattern.lower()
    matches = re.findall(r"([na])(\d*)", pattern)
    if not matches or "".join(m[0] + m[1] for m in matches) != pattern:
        sys.exit(f"Error: Invalid pattern '{pattern}'. Pattern should consist of 'n' and 'a' optionally followed by counts (e.g., 'nanan' or 'n1a3a2n2').")
    seq = []
    for char, count_str in matches:
        count = int(count_str) if count_str else 1
        seq.extend([char] * count)
    return seq


def build_rle_pattern(file_types: list[str]) -> str:
    if not file_types:
        return ""
    result = []
    current_type = file_types[0]
    count = 1
    for t in file_types[1:]:
        if t == current_type:
            count += 1
        else:
            result.append(f"{current_type}{count}")
            current_type = t
            count = 1
    result.append(f"{current_type}{count}")
    return "".join(result)


def get_file_type(f: Path) -> str:
    s = f.stem.lower()
    if "anomaly" in s or "alt" in s or "anom" in s:
        return "a"
    return "n"


def strip_chunk_index(stem: str) -> str:
    """Removes trailing chunk index or embedded chunk index from filenames."""
    # Match patterns like -5-neoaves or -5
    return re.sub(r"-(?:neoaves-)?\d+(?:-[^/]+)?$", "", re.sub(r"-\d+$", "", stem))


def get_chunk_index(stem: str) -> str:
    """Extracts numeric chunk index from filename stem (e.g., '5' from default-500Kb-5 or default-neoaves-500Kb-5)."""
    m = re.search(r"-(\d+)(?:-[^/]+)?$", stem)
    if not m:
        m = re.search(r"-(\d+)-", stem)
    return m.group(1) if m else ""


def build_indexed_pattern(pos_files: list[Path]) -> str:
    parts = []
    for f in pos_files:
        ftype = get_file_type(f)
        idx = get_chunk_index(f.stem)
        parts.append(f"{ftype}{idx}")
    return "".join(parts)


def parse_indexed_pattern(pattern: str) -> list[tuple[str, str]]:
    """
    Parses pattern string such as 'a1n5a2a3n8n1' or 'n1a3a2n2' into tuples of (type, index_str).
    For example, 'a1n5' -> [('a', '1'), ('n', '5')].
    If no index digit is present (e.g. 'n'), index_str is ''.
    """
    pattern = pattern.lower()
    matches = re.findall(r"([na])(\d*)", pattern)
    if not matches or "".join(m[0] + m[1] for m in matches) != pattern:
        sys.exit(f"Error: Invalid pattern '{pattern}'. Pattern should consist of 'n' and 'a' followed by optional chunk numbers (e.g., 'a1n5a2a3n8n1' or 'nanan').")
    return [(char, idx_str) for char, idx_str in matches]


def construct_indexed_path(base_path: Path, idx_str: str) -> Path:
    if not idx_str:
        return base_path
    
    parent = base_path.parent
    suffix = base_path.suffix
    stem = base_path.stem if suffix else base_path.name

    # Check if stem contains a chunk index already, or insert chunk index at end
    if suffix:
        target = parent / f"{stem}-{idx_str}{suffix}"
    else:
        target = parent / f"{stem}-{idx_str}"
    
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate FASTA (or TS) sequences in order, or mix null/alt using a pattern."
    )

    parser.add_argument(
        "positional_args",
        type=str,
        nargs="*",
        help="Sequence files to concatenate in order, or (null_base, alt_base, pattern).",
    )

    parser.add_argument("-n", "--null", type=Path, help="Path to null FASTA file.")
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
        "-p", "--pattern", type=str, help="Pattern string specifying order of sequences (e.g., 'a1n5a2a3n8n1' or 'n1a3a2n2')."
    )
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output FASTA path.")

    args = parser.parse_args()

    output_dir = get_output_dir()
    output_subdir = output_dir / "msa" / "concat"
    output_subdir.mkdir(parents=True, exist_ok=True)

    pos = [Path(p) for p in args.positional_args]

    # Detect if pattern mode is requested
    is_pattern_mode = False
    if args.null and args.alt and args.pattern:
        is_pattern_mode = True
        null_path = args.null
        alt_path = args.alt
        pattern_raw = args.pattern
    elif len(pos) == 3 and not pos[2].exists():
        try:
            parse_indexed_pattern(pos[2].name)
            is_pattern_mode = True
            # Determine which positional argument is null and which is alt based on get_file_type
            p0_type = get_file_type(pos[0])
            p1_type = get_file_type(pos[1])
            if p0_type == "n" and p1_type == "a":
                null_path = pos[0]
                alt_path = pos[1]
            elif p0_type == "a" and p1_type == "n":
                null_path = pos[1]
                alt_path = pos[0]
            else:
                # Default to pos[0] as null, pos[1] as alt if types are identical/ambiguous
                null_path = pos[0]
                alt_path = pos[1]
            pattern_raw = pos[2].name
        except SystemExit:
            is_pattern_mode = False

    if is_pattern_mode:
        pattern_items = parse_indexed_pattern(pattern_raw)

        # Build list of specific files to resolve and read
        files_to_concat = []
        for ftype, idx_str in pattern_items:
            base_p = null_path if ftype == "n" else alt_path
            item_p = construct_indexed_path(base_p, idx_str)
            files_to_concat.append(item_p)

        resolved_files = [resolve_fasta_file(f, output_dir) for f in files_to_concat]
        fastas = []
        for rf in resolved_files:
            print(f"Reading FASTA: {rf}", file=sys.stderr)
            fastas.append(read_fasta(rf))

        base_headers = list(fastas[0].keys())
        for idx, fa in enumerate(fastas[1:], start=2):
            missing = set(base_headers) - set(fa.keys())
            extra = set(fa.keys()) - set(base_headers)
            if missing or extra:
                sys.exit(f"Error: Header mismatch in file #{idx} ({resolved_files[idx-1]}): missing={missing}, extra={extra}")

        clean_null_stem = strip_chunk_index(null_path.stem)
        clean_alt_stem = strip_chunk_index(alt_path.stem)

        if args.output is None:
            out_name = f"{clean_null_stem}_{clean_alt_stem}_{pattern_raw}.fa"
            output_path = output_subdir / out_name
        elif not args.output.is_absolute():
            output_path = output_subdir / args.output.name
        else:
            output_path = args.output

        print(f"Concatenating {len(resolved_files)} sequence files in order (pattern '{pattern_raw}')...", file=sys.stderr)
        with open(output_path, "w") as f:
            for header in base_headers:
                f.write(f">{header}\n")
                seq_concat = "".join(fa[header] for fa in fastas)
                f.write(f"{seq_concat}\n")

        print(f"Writing output FASTA: {output_path}", file=sys.stderr)
        print("Done.", file=sys.stderr)

    elif len(pos) >= 2:
        # Sequential concatenation mode
        file_types = [get_file_type(f) for f in pos]
        pattern_str = build_indexed_pattern(pos)

        resolved_files = [resolve_fasta_file(f, output_dir) for f in pos]
        fastas = []
        for rf in resolved_files:
            print(f"Reading FASTA: {rf}", file=sys.stderr)
            fastas.append(read_fasta(rf))

        base_headers = list(fastas[0].keys())
        for idx, fa in enumerate(fastas[1:], start=2):
            missing = set(base_headers) - set(fa.keys())
            extra = set(fa.keys()) - set(base_headers)
            if missing or extra:
                sys.exit(f"Error: Header mismatch in file #{idx} ({resolved_files[idx-1]}): missing={missing}, extra={extra}")

        first_null = next((rf for rf, ft in zip(resolved_files, file_types) if ft == "n"), None)
        first_alt = next((rf for rf, ft in zip(resolved_files, file_types) if ft == "a"), None)

        clean_null_stem = strip_chunk_index(first_null.stem) if first_null else None
        clean_alt_stem = strip_chunk_index(first_alt.stem) if first_alt else None

        if args.output is None:
            if clean_null_stem and clean_alt_stem:
                out_name = f"{clean_null_stem}_{clean_alt_stem}_{pattern_str}.fa"
            elif clean_null_stem:
                out_name = f"{clean_null_stem}_{pattern_str}.fa"
            else:
                out_name = f"{clean_alt_stem}_{pattern_str}.fa"
            output_path = output_subdir / out_name
        elif not args.output.is_absolute():
            output_path = output_subdir / args.output.name
        else:
            output_path = args.output

        print(f"Concatenating {len(resolved_files)} sequence files in order (pattern '{pattern_str}')...", file=sys.stderr)
        with open(output_path, "w") as f:
            for header in base_headers:
                f.write(f">{header}\n")
                seq_concat = "".join(fa[header] for fa in fastas)
                f.write(f"{seq_concat}\n")

        print(f"Writing output FASTA: {output_path}", file=sys.stderr)
        print("Done.", file=sys.stderr)

    else:
        parser.print_help()
        sys.exit("\nError: Must specify at least 2 files to concatenate, or specify (null, alt, pattern).")
