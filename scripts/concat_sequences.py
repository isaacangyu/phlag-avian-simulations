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


def get_connection_dir() -> Path:
    val = os.environ.get("CONNECTION_DIR")
    if not val:
        env_vars = load_env_vars()
        val = env_vars.get("CONNECTION_DIR")
    if val:
        return Path(val).expanduser().resolve()
    sys.exit("Error: CONNECTION_DIR environment variable is not set. Initialize a shared directory for simulations as input to phlag.")


def get_output_dir() -> Path:
    return get_connection_dir() / "simulations"


def get_phlag_dir() -> Path:
    val = os.environ.get("PHLAG_DIR")
    if not val:
        env_vars = load_env_vars()
        val = env_vars.get("PHLAG_DIR")
    if val:
        p = Path(val).expanduser().resolve()
        if p.exists():
            return p
    cand = (Path(__file__).resolve().parent.parent / "phlag").resolve()
    if cand.exists():
        return cand
    cand_home = (Path.home() / "phlag").resolve()
    if cand_home.exists():
        return cand_home
    sys.exit("Error: PHLAG_DIR environment variable is not set or phlag directory does not exist. Please clone the phlag repository.")


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


def get_experiment_folder_name(path_or_stem: str | Path | None) -> str:
    """
    Extracts experiment folder name from path or stem, e.g.:
      'alt-Nyctibiidae_10X_up-1.ts' -> 'Nyctibiidae_10X_up'
      'anomaly_case-recombination_increase-N36-500Kb-1.ts' -> 'N36_recombination_up'
      'null-1.ts'                  -> 'null'
    """
    if not path_or_stem:
        return "default"
    
    p = Path(path_or_stem)
    if get_base_name(p).isdigit():
        for part in reversed(p.parts[:-1]):
            if part.lower() not in ("ts", "simulated", "concat", "simulations", "store", ".", ".."):
                res = get_experiment_folder_name(part)
                if res and res != "null":
                    return res

    try:
        from simulate import generate_new_filename
        std_name = generate_new_filename(path_or_stem)
        clean = strip_chunk_index(get_base_name(std_name))
    except Exception:
        clean = strip_chunk_index(get_base_name(path_or_stem))

    if clean.startswith("alt-"):
        return clean[4:]
    elif clean.startswith("null-"):
        return "null"
    elif clean in ("null", "alt", "default"):
        return clean
    return clean


def resolve_fasta_file(f: Path, output_dir: Path) -> Path:
    bname = get_base_name(f)
    exp_name = get_experiment_folder_name(f)

    # Standard candidate paths
    fasta_candidates = [
        f if f.suffix in (".fa", ".fasta") else f.with_suffix(".fa"),
        output_dir / exp_name / "simulated" / f"{bname}.fa",
    ]
    
    chunk_idx = get_chunk_index(bname)
    if chunk_idx:
        fasta_candidates.insert(1, output_dir / exp_name / "simulated" / f"{chunk_idx}.fa")

    fasta_candidates.extend([
        output_dir / exp_name / "concat" / f"{bname}.fa",
        output_dir / exp_name / f"{bname}.fa",
        output_dir / "simulated" / f"{bname}.fa",
        f.parent / f"{bname}.fa",
    ])

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


def get_file_type(f: str | Path) -> str:
    s = str(f).lower()
    name = Path(f).name.lower()
    if name.startswith("a") or "anomaly" in s or "alt" in s or "anom" in s or "admixture" in s or "recombination" in s or "10x" in s:
        return "a"
    return "n"


def strip_chunk_index(stem: str) -> str:
    """Removes trailing chunk index or embedded chunk index from base names."""
    clean_stem = get_base_name(Path(stem))
    return re.sub(r"-(?:neoaves-)?\d+(?:-[^/]+)?$", "", re.sub(r"-\d+$", "", clean_stem))


def get_chunk_index(stem: str) -> str:
    """Extracts numeric chunk index from base name (e.g., '5' from null-5 or 'a2' or legacy default-500Kb-5)."""
    m = re.search(r"-(\d+)(?:-[^/]+)?$", stem)
    if not m:
        m = re.search(r"-(\d+)-", stem)
    if not m:
        m = re.search(r"(\d+)", stem)
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
    exp_name = get_experiment_folder_name(base_path)
    found_indices = set()

    search_dirs = [
        output_dir / exp_name / "ts",
        output_dir / exp_name / "simulated",
        output_dir / exp_name / "concat",
        output_dir / exp_name,
        output_dir / "ts",
        Path(base_path).parent if Path(base_path).parent.exists() else Path("."),
        Path("."),
    ]
    
    phlag_dir = get_phlag_dir()
    search_dirs.append(phlag_dir / "ts")

    for sdir in search_dirs:
        if not sdir.exists() or not sdir.is_dir():
            continue
        try:
            for p in sdir.iterdir():
                p_bname = get_base_name(p).lower()
                # If searching for null, match null- or default-
                matches_base = False
                if sdir.name in ("ts", "simulated", "concat") and sdir.parent.name.lower() == exp_name.lower():
                    matches_base = True
                elif sdir.name.lower() == exp_name.lower():
                    matches_base = True
                elif clean_base in ("null", "default") or clean_base.startswith("null") or clean_base.startswith("default"):
                    matches_base = p_bname.startswith("null") or p_bname.startswith("default")
                else:
                    p_clean = strip_chunk_index(p_bname)
                    matches_base = (clean_base in p_bname) or (clean_base in p_clean)

                if matches_base:
                    if p_bname.isdigit():
                        found_indices.add(int(p_bname))
                    else:
                        m = re.search(r"-(\d+)(?:-[^/]+)?$", p_bname)
                        if not m:
                            m = re.search(r"-(\d+)-", p_bname)
                        if m:
                            found_indices.add(int(m.group(1)))
        except Exception:
            pass

    return sorted(list(found_indices))


import random


def parse_proportion(prop_raw: str | float | None) -> float | None:
    if prop_raw is None or prop_raw == "":
        return None
    if isinstance(prop_raw, (int, float)):
        val = float(prop_raw)
    else:
        s = str(prop_raw).strip()
        if "/" in s:
            parts = s.split("/")
            val = float(parts[0]) / float(parts[1])
        else:
            val = float(s)

    if val > 1.0 and val <= 12:
        return val / 12.0
    return val


def parse_alt_pattern(pattern_raw: str | None) -> list[str]:
    """
    Parses alt pattern string (e.g. 'a2a1a3' or 'a1a2a3' or '2,1,3') into list of alt chunk index strings.
    If pattern_raw is None or empty, defaults to increasing chunk indices ['1', '2', ..., '12'].
    """
    if not pattern_raw:
        return [str(i) for i in range(1, 13)]

    pattern_str = str(pattern_raw).strip().lower()
    matches = re.findall(r"a?(\d+)", pattern_str)
    if not matches:
        return [str(i) for i in range(1, 13)]

    return matches


def format_prop_decimal(prop_val: float) -> str:
    val_str = f"{round(prop_val, 4):g}"
    if "." not in val_str:
        val_str = f"{val_str}.0"
    return val_str


def parse_pattern_to_12_chunks(
    pattern_raw: str | None, prop_raw: str | float | None = None
) -> tuple[list[tuple[str, str]], str, float]:
    """
    Parses pattern string and/or proportion into 12 chunk items.
    Returns (pattern_items, suffix_label, prop_val).
      - If pattern_raw contains 'n' and specifies 12 chunks (e.g. 'a1n2a3n4...'), it is used directly.
      - Otherwise, alt chunks are placed via random insertion rule (random start_idx) in 12 chunks.
    """
    prop_val = parse_proportion(prop_raw)
    pattern_str = str(pattern_raw).strip() if pattern_raw else ""

    # 1. Full 12-chunk pattern with nulls specified
    if "n" in pattern_str.lower():
        try:
            pattern_items = parse_indexed_pattern(pattern_str)
            if len(pattern_items) == 12:
                num_alt = sum(1 for t, _ in pattern_items if t == "a")
                if prop_val is None:
                    prop_val = num_alt / 12.0
                suffix_label = build_indexed_pattern([Path(f"{t}{idx}") for t, idx in pattern_items])
                return pattern_items, suffix_label, prop_val
        except SystemExit:
            pass

    # 2. Alt-only or omitted nulls -> random insertion rule
    alt_chunks = parse_alt_pattern(pattern_raw)

    if prop_val is not None:
        num_alt_chunks = int(round(prop_val * 12))
    elif pattern_raw:
        num_alt_chunks = len(alt_chunks)
    else:
        num_alt_chunks = 3  # Default 3 alt chunks out of 12 (0.25 proportion)

    num_alt_chunks = max(0, min(12, num_alt_chunks))
    if prop_val is None:
        prop_val = num_alt_chunks / 12.0

    # Fill remaining up to num_alt_chunks with unused increasing indices if pattern specified fewer
    if len(alt_chunks) < num_alt_chunks:
        used = set(alt_chunks)
        for i in range(1, 13):
            if len(alt_chunks) >= num_alt_chunks:
                break
            if str(i) not in used:
                alt_chunks.append(str(i))

    selected_alt_chunks = alt_chunks[:num_alt_chunks]

    # Random insertion rule: pick a random starting index for the contiguous block of alt chunks
    start_idx = random.randint(0, max(0, 12 - num_alt_chunks))

    pattern_items = []
    alt_ptr = 0
    for i in range(12):
        if start_idx <= i < start_idx + num_alt_chunks:
            pattern_items.append(("a", selected_alt_chunks[alt_ptr]))
            alt_ptr += 1
        else:
            pattern_items.append(("n", str(i + 1)))

    suffix_label = build_indexed_pattern([Path(f"{t}{idx}") for t, idx in pattern_items])
    return pattern_items, suffix_label, prop_val


def expand_pattern_to_target_files(
    null_path: Path, alt_path: Path, pattern_raw: str | None, output_dir: Path, proportion: str | float | None = None
) -> tuple[list[Path], str, float]:
    """
    Expands a pattern string and proportion into concrete target file paths for 12 chunks (6Mbp total).
    Returns (target_files, suffix_label, prop_val).
    """
    pattern_items, suffix_label, prop_val = parse_pattern_to_12_chunks(pattern_raw, prop_raw=proportion)
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
            target_files.append(construct_indexed_path(base_p, num_str))
            if indices and int(num_str) in indices:
                ptr = (indices.index(int(num_str)) + 1) % len(indices)
            else:
                ptr += 1

        if ftype == "n":
            null_ptr = ptr
        else:
            alt_ptr = ptr

    return target_files, suffix_label, prop_val


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Concatenate FASTA (or TS) sequences into a 6Mbp sequence (12 chunks) using pattern or random alternative insertion."
    )

    parser.add_argument(
        "positional_args",
        type=str,
        nargs="*",
        help="Sequence files to concatenate: (null, alt, [pattern_or_proportion]) or 12 individual files.",
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
        "-p", "--prop", "--proportion",
        dest="proportion",
        type=str,
        default=None,
        help="Proportion of alternative sequence (e.g. 0.25, 0.5, 3/12)."
    )
    parser.add_argument(
        "--pattern", type=str, default=None, help="Alternative pattern string specifying alt chunk order (e.g. 'a2a1a3' or 'a1a2a3')."
    )
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output FASTA path.")

    args = parser.parse_args(argv)

    output_dir = get_output_dir()

    pos = [Path(p) for p in args.positional_args]

    null_path = None
    alt_path = None
    pattern_raw = args.pattern
    proportion_raw = args.proportion
    is_pair_mode = False

    if args.null and args.alt:
        is_pair_mode = True
        null_path = args.null
        alt_path = args.alt
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
            p3_str = pos[2].name
            if parse_proportion(p3_str) is not None:
                proportion_raw = p3_str
            else:
                pattern_raw = p3_str
    elif len(pos) == 1 or (args.alt and not args.null):
        is_pair_mode = True
        alt_path = args.alt if args.alt else pos[0]
        if null_path is None:
            null_path = Path("null")

    if is_pair_mode:
        files_to_concat, suffix_label, prop_val = expand_pattern_to_target_files(
            null_path, alt_path, pattern_raw, output_dir, proportion=proportion_raw
        )
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

        exp_folder_name = get_experiment_folder_name(alt_path if alt_path else null_path)
        exp_dir = output_dir / exp_folder_name

        prop_dec_str = format_prop_decimal(prop_val)
        concat_dir = exp_dir / f"concat_{prop_dec_str}"
        concat_dir.mkdir(parents=True, exist_ok=True)

        file_types = [get_file_type(f) for f in files_to_concat]
        pattern_str = build_indexed_pattern(files_to_concat)
        pattern_txt_content = pattern_str if pattern_str else (suffix_label if suffix_label else "")

        if args.output is None:
            output_name = f"{pattern_txt_content}_{prop_dec_str}.fa" if pattern_txt_content else f"{exp_folder_name}_{prop_dec_str}.fa"
            output_path = concat_dir / output_name
        elif not args.output.is_absolute():
            output_name = args.output.name
            if not output_name.endswith(f"_{prop_dec_str}.fa") and not output_name.endswith(f"_{prop_dec_str}.fasta"):
                stem = get_base_name(output_name)
                ext = Path(output_name).suffix or ".fa"
                output_name = f"{stem}_{prop_dec_str}{ext}"
            output_path = concat_dir / output_name
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

        first_alt = next((rf for rf, ft in zip(resolved_files, file_types) if ft == "a"), None)
        first_null = next((rf for rf, ft in zip(resolved_files, file_types) if ft == "n"), None)

        exp_folder_name = get_experiment_folder_name(first_alt if first_alt else first_null)
        exp_dir = output_dir / exp_folder_name
        concat_dir = exp_dir / "concat"
        concat_dir.mkdir(parents=True, exist_ok=True)

        if args.output is None:
            output_name = f"{pattern_str}.fa" if pattern_str else f"{exp_folder_name}.fa"
            output_path = concat_dir / output_name
        elif not args.output.is_absolute():
            output_path = concat_dir / args.output.name
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

