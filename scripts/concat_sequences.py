from __future__ import annotations

import os
import sys
import argparse
from typing import Union
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPTS_DIR.parent


def load_env_vars() -> dict[str, str]:
    env_vars = {}
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
                            env_vars[key.strip()] = val.strip().strip("'").strip('"')
            except Exception:
                pass
            break
        current = current.parent
    return env_vars


def get_large_files_dir() -> Path:
    val = os.environ.get("LARGE_FILES_DIR")
    if not val:
        env_vars = load_env_vars()
        val = env_vars.get("LARGE_FILES_DIR") or env_vars.get("PHLAG_DIR") or os.environ.get("PHLAG_DIR")
    if val:
        return Path(val)
    sys.exit("Error: LARGE_FILES_DIR environment variable is not set and could not be found in .env file.")


def get_output_dir() -> Path:
    return get_large_files_dir()


def get_phlag_dir() -> Path:
    val = os.environ.get("PHLAG_DIR")
    if not val:
        env_vars = load_env_vars()
        val = env_vars.get("PHLAG_DIR")
    if val and val != str(get_large_files_dir()):
        return Path(val)
    cand = Path("/home/iang/phlag")
    if cand.exists():
        return cand
    return Path(__file__).resolve().parent.parent


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


def get_base_name(path: str | Path | None) -> str:
    """Extracts base name from a file path or string, stripping .ts, .fa, or .fasta extensions if present."""
    if path is None:
        return ""
    name = Path(path).name
    for ext in (".ts", ".fa", ".fasta"):
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return Path(path).stem


def resolve_fasta_file(f: Path, output_dir: Path) -> Path:
    bname = get_base_name(f)

    # Standard candidate paths
    fasta_candidates = [
        f if f.suffix in (".fa", ".fasta") else f.with_suffix(".fa"),
        output_dir / "msa" / "simulated" / f"{bname}.fa",
        output_dir / "msa" / "concat" / f"{bname}.fa",
        output_dir / "msa" / "simulated" / f.name,
        output_dir / "msa" / "concat" / f.name,
        f.parent / f"{bname}.fa",
    ]

    for candidate in fasta_candidates:
        if candidate.exists() and candidate.is_file():
            return candidate

    sys.exit(f"Error: Sequence FASTA file does not exist for base name: '{bname}' ({f})")


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
    s = get_base_name(f).lower()
    if "anomaly" in s or "alt" in s or "anom" in s:
        return "a"
    return "n"


def strip_chunk_index(stem: str) -> str:
    """Removes trailing chunk index or embedded chunk index from base names."""
    clean_stem = get_base_name(Path(stem))
    return re.sub(r"-(?:neoaves-)?\d+(?:-[^/]+)?$", "", re.sub(r"-\d+$", "", clean_stem))


def get_chunk_index(stem: str) -> str:
    """Extracts numeric chunk index from base name (e.g., '5' from null-5 or legacy default-500Kb-5)."""
    m = re.search(r"-(\d+)(?:-[^/]+)?$", stem)
    if not m:
        m = re.search(r"-(\d+)-", stem)
    return m.group(1) if m else ""


def build_indexed_pattern(pos_files: list[Path]) -> str:
    parts = []
    for f in pos_files:
        ftype = get_file_type(f)
        idx = get_chunk_index(get_base_name(f))
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


def construct_indexed_path(base_name: str | Path, idx_str: str) -> Path:
    base_p = Path(base_name)
    clean_bname = strip_chunk_index(get_base_name(base_p))
    if not idx_str:
        return base_p.parent / clean_bname
    return base_p.parent / f"{clean_bname}-{idx_str}"


def find_existing_chunk_indices(base_path: str | Path, output_dir: Path) -> list[int]:
    """
    Finds all existing numeric chunk indices for a base name across ts, simulated, concat, and local directories.
    Returns a sorted list of unique integer chunk indices.
    """
    clean_base = strip_chunk_index(get_base_name(base_path)).lower()
    found_indices = set()

    search_dirs = [
        output_dir / "ts",
        REPO_ROOT / "store" / "ts",
        output_dir / "msa" / "simulated",
        output_dir / "msa" / "concat",
        Path(base_path).parent if Path(base_path).parent.exists() else Path("."),
        Path("."),
    ]

    for sdir in search_dirs:
        if not sdir.exists() or not sdir.is_dir():
            continue
        try:
            for p in sdir.iterdir():
                p_bname = get_base_name(p).lower()
                # If searching for null, match null- or default-
                matches_base = False
                if clean_base in ("null", "default") or clean_base.startswith("null") or clean_base.startswith("default"):
                    matches_base = p_bname.startswith("null") or p_bname.startswith("default")
                else:
                    p_clean = strip_chunk_index(p_bname)
                    matches_base = (clean_base in p_bname) or (clean_base in p_clean)

                if matches_base:
                    m = re.search(r"-(\d+)(?:-[^/]+)?$", p_bname)
                    if not m:
                        m = re.search(r"-(\d+)-", p_bname)
                    if m:
                        found_indices.add(int(m.group(1)))
        except Exception:
            pass

    return sorted(list(found_indices))


import random


def parse_pattern_to_12_chunks(pattern_raw: str | None) -> tuple[list[tuple[str, str]], str]:
    """
    Parses a pattern string, single count digit, or None into 12 chunk items.
    Returns (pattern_items, suffix_label).
      - If pattern_raw is None or empty: 3 alt chunks inserted randomly into 12, suffix_label = ""
      - If pattern_raw is a digit string (e.g., '3' or '5'): n alt chunks inserted randomly into 12, suffix_label = "{n}a"
      - If pattern_raw is a pattern string: validated to ensure exactly 12 chunks total, suffix_label = pattern_raw
    """
    if pattern_raw is None or pattern_raw == "":
        alt_positions = set(random.sample(range(12), 3))
        pattern_items = [("a" if i in alt_positions else "n", "") for i in range(12)]
        return pattern_items, ""

    pattern_raw_str = str(pattern_raw).strip()
    if pattern_raw_str.isdigit():
        n_alt = int(pattern_raw_str)
        if not (0 <= n_alt <= 12):
            sys.exit(f"Error: Number of alternative chunks must be between 0 and 12, got {n_alt}.")
        alt_positions = set(random.sample(range(12), n_alt))
        pattern_items = [("a" if i in alt_positions else "n", "") for i in range(12)]
        return pattern_items, f"{n_alt}a"

    # Parse indexed or RLE/character pattern string
    pattern_items = parse_indexed_pattern(pattern_raw_str)
    
    total_chunks = len(pattern_items)
    if total_chunks != 12:
        sys.exit(f"Error: Pattern must specify exactly 12 chunks total for a 6Mbp sequence (got {total_chunks} chunks in '{pattern_raw_str}').")

    return pattern_items, pattern_raw_str


def expand_pattern_to_target_files(
    null_path: Path, alt_path: Path, pattern_raw: str | None, output_dir: Path
) -> tuple[list[Path], str]:
    """
    Expands a pattern string into concrete target file paths for 12 chunks (6Mbp total).
    Returns (target_files, suffix_label).
    """
    pattern_items, suffix_label = parse_pattern_to_12_chunks(pattern_raw)
    null_indices = find_existing_chunk_indices(null_path, output_dir)
    alt_indices = find_existing_chunk_indices(alt_path, output_dir)

    null_ptr = 0
    alt_ptr = 0
    target_files = []

    for ftype, num_str in pattern_items:
        base_p = null_path if ftype == "n" else alt_path
        indices = null_indices if ftype == "n" else alt_indices
        ptr = null_ptr if ftype == "n" else alt_ptr

        if not num_str:
            if indices:
                idx_val = str(indices[ptr % len(indices)])
            else:
                idx_val = str(ptr + 1)
            ptr += 1
            target_files.append(construct_indexed_path(base_p, idx_val))
        else:
            val = int(num_str)
            if indices and val in indices:
                target_files.append(construct_indexed_path(base_p, num_str))
                ptr = (indices.index(val) + 1) % len(indices)
            else:
                if indices:
                    idx_val = str(indices[ptr % len(indices)])
                else:
                    idx_val = str(ptr + 1)
                ptr += 1
                target_files.append(construct_indexed_path(base_p, idx_val))

        if ftype == "n":
            null_ptr = ptr
        else:
            alt_ptr = ptr

    return target_files, suffix_label


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Concatenate FASTA (or TS) sequences into a 6Mbp sequence (12 chunks) using pattern or random alternative insertion."
    )

    parser.add_argument(
        "positional_args",
        type=str,
        nargs="*",
        help="Sequence files to concatenate: (null, alt, [pattern_or_count]) or 12 individual files.",
    )

    parser.add_argument("-n", "--null", type=Path, help="Base name or path to null sequence.")
    parser.add_argument(
        "-a",
        "--alt",
        "--alternative",
        "--anomalous",
        type=Path,
        dest="alt",
        help="Base name or path to alternative sequence.",
    )
    parser.add_argument(
        "-p", "--pattern", type=str, help="Pattern string (e.g. 'aanannnnaann'), single count n (e.g. '3'), or omitted for default 3 alt chunks."
    )
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output FASTA path.")

    args = parser.parse_args(argv)

    output_dir = get_output_dir()
    output_subdir = output_dir / "msa" / "concat"
    output_subdir.mkdir(parents=True, exist_ok=True)

    pos = [Path(p) for p in args.positional_args]

    null_path = None
    alt_path = None
    pattern_raw = None
    is_pair_mode = False

    if args.null and args.alt:
        is_pair_mode = True
        null_path = args.null
        alt_path = args.alt
        pattern_raw = args.pattern
    elif len(pos) >= 2 and len(pos) <= 3:
        is_pair_mode = True
        p0_type = get_file_type(pos[0])
        p1_type = get_file_type(pos[1])
        if p0_type == "n" and p1_type == "a":
            null_path = pos[0]
            alt_path = pos[1]
        elif p0_type == "a" and p1_type == "n":
            null_path = pos[1]
            alt_path = pos[0]
        else:
            null_path = pos[0]
            alt_path = pos[1]

        if len(pos) == 3:
            pattern_raw = pos[2].name
    elif len(pos) == 1 or (args.alt and not args.null):
        is_pair_mode = True
        alt_path = args.alt if args.alt else pos[0]
        pattern_raw = args.pattern if args.pattern else None
        
        if null_path is None:
            sys.exit("Error: Null sequence path (-n/--null) must be provided.")

    if is_pair_mode:
        files_to_concat, suffix_label = expand_pattern_to_target_files(null_path, alt_path, pattern_raw, output_dir)
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

        clean_null_base = strip_chunk_index(get_base_name(null_path))
        clean_alt_base = strip_chunk_index(get_base_name(alt_path))

        if args.output is None:
            if suffix_label:
                out_name = f"{clean_null_base}_{clean_alt_base}_{suffix_label}.fa"
            else:
                out_name = f"{clean_null_base}_{clean_alt_base}.fa"
            output_path = output_subdir / out_name
        elif not args.output.is_absolute():
            output_path = output_subdir / args.output.name
        else:
            output_path = args.output

        print(f"Concatenating 12 sequence files (6Mbp total)...", file=sys.stderr)
        with open(output_path, "w") as f:
            for header in base_headers:
                f.write(f">{header}\n")
                seq_concat = "".join(fa[header] for fa in fastas)
                f.write(f"{seq_concat}\n")

        print(f"Writing output FASTA: {output_path}", file=sys.stderr)
        print("Done.", file=sys.stderr)

    elif len(pos) == 12:
        # Sequential concatenation of 12 explicit files
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

        clean_null_base = strip_chunk_index(get_base_name(first_null)) if first_null else None
        clean_alt_base = strip_chunk_index(get_base_name(first_alt)) if first_alt else None

        if args.output is None:
            if clean_null_base and clean_alt_base:
                out_name = f"{clean_null_base}_{clean_alt_base}_{pattern_str}.fa"
            elif clean_null_base:
                out_name = f"{clean_null_base}_{pattern_str}.fa"
            else:
                out_name = f"{clean_alt_base}_{pattern_str}.fa"
            output_path = output_subdir / out_name
        elif not args.output.is_absolute():
            output_path = output_subdir / args.output.name
        else:
            output_path = args.output

        print(f"Concatenating 12 sequence files (6Mbp total)...", file=sys.stderr)
        with open(output_path, "w") as f:
            for header in base_headers:
                f.write(f">{header}\n")
                seq_concat = "".join(fa[header] for fa in fastas)
                f.write(f"{seq_concat}\n")

        print(f"Writing output FASTA: {output_path}", file=sys.stderr)
        print("Done.", file=sys.stderr)

    else:
        parser.print_help()
        sys.exit("\nError: Must specify (null, alt, [pattern_or_count]) or exactly 12 files to concatenate a 6Mbp sequence.")


if __name__ == "__main__":
    main()

