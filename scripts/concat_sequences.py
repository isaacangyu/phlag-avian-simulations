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
      'alt-Nyctibiidae_10X_up-1.ts' -> 'Nyctibiidae'
      'N554_10X_up-1.ts' -> 'N554'
      'Strigiformes_N297_admixture_rate090-time7099554-1.ts' -> 'Strigiformes_N297_rate090-time7099554'
      'null-1.ts' -> 'null'
    """
    if not path_or_stem:
        return "default"
    
    p = Path(path_or_stem)
    reserved = {"ts", "simulated", "concat", "simulations", "store", ".", "..", "up", "down", "low", "high", "10x", "recombination", "admixture", "iang", "home", "drive2", "root", "none"}

    bname = get_base_name(p)
    if bname.isdigit() or bname.lower() in reserved:
        for part in reversed(p.parts[:-1]):
            part_clean = part.strip("-_")
            if part_clean.lower() not in reserved and not part_clean.isdigit():
                res = get_experiment_folder_name(part)
                if res and res != "default":
                    return res
        return "default"

    std_name = None
    if "simulate" in sys.modules and hasattr(sys.modules["simulate"], "generate_new_filename"):
        std_name = sys.modules["simulate"].generate_new_filename(path_or_stem)
    else:
        try:
            from simulate import generate_new_filename
            std_name = generate_new_filename(path_or_stem)
        except Exception:
            pass

    if std_name:
        clean = strip_chunk_index(get_base_name(std_name))
    else:
        clean = strip_chunk_index(get_base_name(path_or_stem))

    if clean.startswith("alt-"):
        clean = clean[4:]
    elif clean.startswith("null-") or clean == "null":
        return "null"

    clean = clean.replace("anomaly_case-", "").replace("alt-", "")
    
    meta = {}
    if "simulate" in sys.modules and hasattr(sys.modules["simulate"], "parse_filename"):
        meta = sys.modules["simulate"].parse_filename(p)
    else:
        try:
            from simulate import parse_filename
            meta = parse_filename(p)
        except Exception:
            pass

    branch = meta.get("branch")
    change = meta.get("change")
    rate = meta.get("rate")
    time = meta.get("time")

    if change == "admixture" and branch and rate and time:
        rate_str = rate if rate.startswith("rate") else f"rate{rate}"
        time_str = time if time.startswith("time") else f"time{time}"
        return f"{branch}_{rate_str}-{time_str}"
    elif branch and branch.lower() not in reserved and not branch.isdigit():
        return branch

    clean = clean.replace("_admixture", "").replace("admixture_", "").replace("admixture", "")
    clean = clean.replace("_10X_up", "").replace("_10X_down", "").replace("_10x_up", "").replace("_10x_down", "")
    clean = clean.replace("_recombination_up", "").replace("_recombination_down", "").replace("_recombination_increase", "").replace("_recombination_decrease", "")
    clean = clean.replace("_up", "").replace("_down", "")
    clean = clean.strip("-_")

    if clean.lower() in reserved:
        return "default"

    return clean if clean else "default"


def get_experiment_relative_path(path_or_stem: str | Path | None) -> Path:
    """
    Returns relative directory path under simulations/ with two levels of categories:
      - 10X/up/<exp_name> or 10X/down/<exp_name>
      - recombination/up/<exp_name> or recombination/down/<exp_name>
      - admixture/low/<exp_name> (time < 4) or admixture/high/<exp_name> (time >= 4)
      - null -> null
    """
    if not path_or_stem:
        return Path("default")

    p = Path(path_or_stem)
    exp_name = get_experiment_folder_name(p)
    if exp_name in ("null", "default"):
        return Path(exp_name)

    meta = {}
    if "simulate" in sys.modules and hasattr(sys.modules["simulate"], "parse_filename"):
        meta = sys.modules["simulate"].parse_filename(p)
    else:
        try:
            from simulate import parse_filename
            meta = parse_filename(p)
        except Exception:
            pass

    change = meta.get("change") or ""
    str_lower = str(p).lower()
    exp_lower = exp_name.lower()
    change_lower = change.lower()
    time_str = meta.get("time") or ""

    # 1. 10X
    if "10x" in str_lower or "10x" in change_lower:
        sub = "up" if ("up" in str_lower or "up" in change_lower) else "down"
        return Path("10X") / sub / exp_name

    # 2. Recombination
    if "recombination" in str_lower or "recombination" in change_lower:
        is_down = "suppression" in str_lower or "decrease" in str_lower or "down" in str_lower or "down" in change_lower
        sub = "down" if is_down else "up"
        return Path("recombination") / sub / exp_name

    # 3. Admixture
    if "admixture" in str_lower or change_lower == "admixture":
        div_time = None
        if time_str:
            m_num = re.search(r"\d+", time_str)
            if m_num:
                div_time = float(m_num.group())
        if div_time is None:
            m_time = re.search(r"time[-_]?(\d+)", str_lower)
            if m_time:
                div_time = float(m_time.group(1))

        if div_time is not None:
            if div_time > 1000:
                div_time = div_time / 1000000.0
            sub = "low" if div_time < 4.0 else "high"
        else:
            sub = "low"
        return Path("admixture") / sub / exp_name

    return Path(exp_name)


def resolve_fasta_file(f: Path, output_dir: Path) -> Path:
    bname = get_base_name(f)
    exp_name = get_experiment_folder_name(f)
    exp_rel = get_experiment_relative_path(f)

    # Standard candidate paths (checking 2-level path first, then legacy flat path)
    fasta_candidates = [
        f if f.suffix in (".fa", ".fasta") else f.with_suffix(".fa"),
        output_dir / exp_rel / "simulated" / f"{bname}.fa",
        output_dir / exp_name / "simulated" / f"{bname}.fa",
    ]
    
    chunk_idx = get_chunk_index(bname)
    if chunk_idx:
        fasta_candidates.insert(1, output_dir / exp_rel / "simulated" / f"{chunk_idx}.fa")
        fasta_candidates.insert(2, output_dir / exp_name / "simulated" / f"{chunk_idx}.fa")

    fasta_candidates.extend([
        output_dir / exp_rel / "concat" / f"{bname}.fa",
        output_dir / exp_name / "concat" / f"{bname}.fa",
        output_dir / exp_rel / f"{bname}.fa",
        output_dir / exp_name / f"{bname}.fa",
        output_dir / "simulated" / f"{bname}.fa",
        f.parent / f"{bname}.fa",
    ])

    for candidate in fasta_candidates:
        if candidate.exists() and candidate.is_file():
            return candidate

    sys.exit(f"Error: Sequence FASTA file does not exist for base name: '{bname}' ({f})")


def ensure_null_12mb_block(output_dir: Path) -> Path:
    """
    Creates/retrieves a 12Mbps (12 x 500Kb chunks) null concatenated block saved in output_dir/null/n1n2n3n4n5n6n7n8n9n10n11n12.fa.
    """
    null_dir = output_dir / "null" / "concat"
    null_dir.mkdir(parents=True, exist_ok=True)
    null_12mb_path = null_dir / "n1n2n3n4n5n6n7n8n9n10n11n12.fa"
    legacy_path = output_dir / "null" / "null_12mb.fa"

    if null_12mb_path.exists():
        return null_12mb_path
    if legacy_path.exists():
        return legacy_path

    # Concatenate null chunks 1..12
    null_chunk_files = [resolve_fasta_file(Path(f"null-{i}"), output_dir) for i in range(1, 13)]
    fastas = [read_fasta(f) for f in null_chunk_files]

    base_headers = list(fastas[0].keys())
    for idx, fa in enumerate(fastas[1:], start=2):
        missing = set(base_headers) - set(fa.keys())
        extra = set(fa.keys()) - set(base_headers)
        if missing or extra:
            sys.exit(f"Error: Header mismatch in null chunk #{idx} ({null_chunk_files[idx-1]}): missing={missing}, extra={extra}")

    print(f"Creating 12Mbps null block -> {null_12mb_path}...", file=sys.stderr)
    with open(null_12mb_path, "w") as f:
        for header in base_headers:
            f.write(f">{header}\n")
            seq_concat = "".join(fa[header] for fa in fastas)
            f.write(f"{seq_concat}\n")

    return null_12mb_path


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
    """Extracts numeric chunk index (1-12) from base name (e.g., '1' from '...-500Kb-1' or 'null-5' or '1.ts')."""
    s = Path(stem).stem if ("/" in stem or "\\" in stem or stem.endswith(".ts") or stem.endswith(".fa")) else stem
    # 1. Matches explicit trailing chunk index after dash (e.g. -1, -2, -500Kb-1)
    m = re.search(r"-(\d{1,2})$", s)
    if not m:
        # 2. Matches index before extension or after 500Kb
        m = re.search(r"(?:-500[Kk][Bb]|-null|-alt|-default)-(\d{1,2})$", s)
    if not m:
        # 3. Matches standalone single/double digit index (e.g., "1", "2", "3")
        m = re.fullmatch(r"(\d{1,2})", s)
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

    exp_rel = get_experiment_relative_path(base_path)
    search_dirs = [
        output_dir / exp_rel / "ts",
        output_dir / exp_rel / "simulated",
        output_dir / exp_rel / "concat",
        output_dir / exp_rel,
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


def build_interleaved_pattern_string(np_pattern: str, ap_pattern: str, intervals: list[tuple[float, float]] | None) -> str:
    np_matches = re.findall(r'n\d+', np_pattern.lower())
    null_chunks = np_matches if np_matches else [f"n{i}" for i in range(1, 13)]

    ap_matches = re.findall(r'a\d+', ap_pattern.lower())
    alt_chunks = ap_matches if ap_matches else [f"a{i}" for i in range(1, 4)]

    if not intervals:
        intervals = [(0.0, 25.0)]

    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    
    total_alt_pct = sum(e - s for s, e in sorted_intervals)

    M = len(null_chunks)
    null_chunk_starts = [i * (100.0 / M) for i in range(M)]

    K = len(alt_chunks)
    alt_chunk_starts = [j * (total_alt_pct / K) for j in range(K)] if total_alt_pct > 0 else []

    tokens = []
    last_null_idx = -1
    last_alt_idx = -1

    current_null_prog = 0.0
    current_alt_prog = 0.0

    curr_pos = 0.0
    for st, end in sorted_intervals:
        null_len = st - curr_pos
        if null_len > 0:
            null_end_prog = current_null_prog + null_len
            for i in range(M):
                if i > last_null_idx and null_chunk_starts[i] < null_end_prog - 1e-7:
                    tokens.append(null_chunks[i])
                    last_null_idx = i
            current_null_prog = null_end_prog

        alt_len = end - st
        if alt_len > 0 and K > 0:
            alt_end_prog = current_alt_prog + alt_len
            for j in range(K):
                if j > last_alt_idx and alt_chunk_starts[j] < alt_end_prog - 1e-7:
                    tokens.append(alt_chunks[j])
                    last_alt_idx = j
            current_alt_prog = alt_end_prog

        curr_pos = end

    remaining_null_len = 100.0 - curr_pos
    if remaining_null_len > 0:
        null_end_prog = current_null_prog + remaining_null_len
        for i in range(M):
            if i > last_null_idx and null_chunk_starts[i] < null_end_prog - 1e-7:
                tokens.append(null_chunks[i])
                last_null_idx = i

    return "".join(tokens)


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


def parse_interval_pattern(pattern_raw: str) -> list[tuple[int, int]] | None:
    """
    Parses start1_end1_start2_end2... pattern (e.g. '500_1000_1500_2000' or '500kb_1000kb')
    into list of (start_bp, end_bp) tuples.
    Returns None if pattern_raw is not an interval pattern.
    """
    if not pattern_raw:
        return None
    s = pattern_raw.strip().lower()
    # Check if string matches start_end format with underscores/dashes/commas
    parts = [p.strip() for p in re.split(r"[_,-]", s) if p.strip()]
    if not parts or len(parts) % 2 != 0:
        return None

    intervals = []
    for i in range(0, len(parts), 2):
        p_start, p_end = parts[i], parts[i + 1]
        try:
            mult_start = 1000 if "kb" in p_start else (1000000 if "mb" in p_start else 1)
            mult_end = 1000 if "kb" in p_end else (1000000 if "mb" in p_end else 1)
            val_start = int(float(re.sub(r"[^\d.]", "", p_start)) * mult_start)
            val_end = int(float(re.sub(r"[^\d.]", "", p_end)) * mult_end)
            if val_start < val_end:
                intervals.append((val_start, val_end))
            else:
                return None
        except ValueError:
            return None

    return intervals if intervals else None


def format_prop_decimal(prop_val: float) -> str:
    val_str = f"{round(prop_val, 4):g}"
    if "." not in val_str:
        val_str = f"{val_str}.0"
    return val_str


def parse_pattern_to_12_chunks(
    pattern_raw: str | None, prop_raw: str | float | None = None
) -> tuple[list[tuple[str, str]], str, float]:
    """
    Parses chunk-based pattern string and/or proportion into 12 chunk items.
    Returns (pattern_items, suffix_label, prop_val).
    """
    prop_val = parse_proportion(prop_raw)
    pattern_str = str(pattern_raw).strip() if pattern_raw else ""

    # Check if pattern string lacks chunk indices (e.g. 'a', 'a1a2', etc.)
    has_explicit_indices = bool(re.search(r"[na]\d+", pattern_str.lower()))

    # If pattern string has no chunk indices, use 12Mbps null concatenated block for insertion
    if pattern_str and not has_explicit_indices and "n" not in pattern_str.lower():
        # User specified alt pattern without chunk indices (e.g. 'a') or positional intervals
        alt_chunks = parse_alt_pattern(pattern_raw)
        if prop_val is not None:
            num_alt_chunks = int(round(prop_val * 12))
        else:
            num_alt_chunks = len(alt_chunks)
        num_alt_chunks = max(1, min(12, num_alt_chunks))
        if prop_val is None:
            prop_val = num_alt_chunks / 12.0

        pattern_items = [("n1n2n3n4n5n6n7n8n9n10n11n12", "")]
        for i in range(1, num_alt_chunks + 1):
            pattern_items.append(("a", str(i)))
        suffix_label = f"n1n2n3n4n5n6n7n8n9n10n11n12_{pattern_str}" if pattern_str else f"n1n2n3n4n5n6n7n8n9n10n11n12_a{num_alt_chunks}"
        return pattern_items, suffix_label, prop_val

    # 1. Pattern with explicit 'n' specified
    if "n" in pattern_str.lower():
        try:
            parsed = parse_indexed_pattern(pattern_str)
            if len(parsed) <= 12:
                pattern_items = list(parsed)
                # Fill remaining positions up to 12 with increasing null chunk indices ('n')
                if len(pattern_items) < 12:
                    curr_pos = len(pattern_items)
                    for idx in range(curr_pos + 1, 13):
                        pattern_items.append(("n", str(idx)))
                
                # Assign default indices to any unindexed 'n' or 'a' entries based on position/occurrence
                final_items = []
                for i, (t, idx) in enumerate(pattern_items, start=1):
                    if not idx:
                        idx = str(i)
                    final_items.append((t, idx))
                pattern_items = final_items

                num_alt = sum(1 for t, _ in pattern_items if t == "a")
                if prop_raw is not None and prop_val is not None:
                    requested_alt_chunks = int(round(prop_val * 12))
                    if requested_alt_chunks > num_alt:
                        sys.exit(f"Error: Specified proportion {prop_raw} ({requested_alt_chunks} alt chunk(s)) is greater than the number of alternative ('a') chunks in pattern '{pattern_raw}' ({num_alt} chunk(s)).")
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
        if pattern_raw and len(alt_chunks) < num_alt_chunks:
            sys.exit(f"Error: Specified proportion {prop_raw} ({num_alt_chunks} alt chunk(s)) is greater than the number of alternative ('a') chunks in pattern '{pattern_raw}' ({len(alt_chunks)} chunk(s)).")
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

    # Alt chunk insertion rule: start alt chunks from position 0
    start_idx = 0

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
    null_indices = find_existing_chunk_indices(null_path if null_path else Path("null"), output_dir)
    alt_indices = find_existing_chunk_indices(alt_path, output_dir) if alt_path else []

    null_ptr = 0
    alt_ptr = 0
    target_files = []

    for ftype, num_str in pattern_items:
        base_p = (null_path if null_path else Path("null")) if ftype == "n" else (alt_path if alt_path else Path("alt"))
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
        "--np", "--null-pattern", type=str, default="n1n2n3n4n5n6n7n8n9n10n11n12", help="Null pattern string specifying null chunk order."
    )
    parser.add_argument(
        "--ap", "--alt-pattern", type=str, default="a1a2a3", help="Alternative pattern string specifying alt chunk order."
    )
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output FASTA path.")

    args = parser.parse_args(argv)

    output_dir = get_output_dir()

    pos = [Path(p) for p in args.positional_args]

    is_pair_mode = False
    null_path = args.null
    alt_path = args.alt

    # Extract any positional arguments that are intervals
    intervals_raw = []
    file_pos = []
    for p in pos:
        if bool(re.match(r"^[\d\.,\-_kmb]+$", p.name.lower().strip())):
            intervals_raw.append(p.name)
        else:
            file_pos.append(p)
    
    args.interval = ",".join(intervals_raw) if intervals_raw else None

    if null_path or alt_path or args.interval or len(file_pos) <= 3:
        is_pair_mode = True
        if file_pos and not null_path:
            null_path = file_pos[0]
        if len(file_pos) > 1 and not alt_path:
            alt_path = file_pos[1]
        
        if null_path is None:
            null_path = Path("null")

    if is_pair_mode:
        intervals = parse_interval_pattern(args.interval) if args.interval else None
        
        # 1. Load the null base sequence
        if null_path and str(null_path).lower() != "null":
            base_null_fa = resolve_fasta_file(null_path, output_dir)
            null_fastas = read_fasta(base_null_fa)
        else:
            matches = re.findall(r'[an](\d+)', args.np.lower())
            null_chunk_indices = [int(m) for m in matches] if matches else [i for i in range(1, 13)]
            
            null_fastas = {}
            for idx in null_chunk_indices:
                cf = output_dir / "null" / "simulated" / f"{idx}.fa"
                if not cf.exists():
                    cf = ensure_null_12mb_block(output_dir)
                    null_fastas = read_fasta(cf)
                    break
                else:
                    rf = read_fasta(cf)
                    for header, seq in rf.items():
                        if header not in null_fastas:
                            null_fastas[header] = ""
                        null_fastas[header] += seq

        base_headers = list(null_fastas.keys())

        # 2. Load the alt chunks
        alt_fastas = {}
        if alt_path and str(alt_path).lower() != "null":
            matches = re.findall(r'[an](\d+)', args.ap.lower())
            alt_chunk_indices = [int(m) for m in matches] if matches else [1, 2, 3]
            
            for idx in alt_chunk_indices:
                cf = construct_indexed_path(alt_path, str(idx))
                rf = resolve_fasta_file(cf, output_dir)
                chunk_data = read_fasta(rf)
                for header in base_headers:
                    if header not in alt_fastas:
                        alt_fastas[header] = ""
                    alt_fastas[header] += chunk_data.get(header, "")

        exp_rel_path = get_experiment_relative_path(alt_path if alt_path and str(alt_path).lower() != "null" else null_path)
        exp_folder_name = exp_rel_path.name
        exp_dir = output_dir / exp_rel_path
        concat_dir = exp_dir / "concat"
        concat_dir.mkdir(parents=True, exist_ok=True)

        if intervals:
            total_null_bp = sum(len(seq) for seq in null_fastas.values()) // len(base_headers) if base_headers else 6000000
            
            parsed_intervals = []
            pct_intervals = []
            for start_val, end_val in intervals:
                if end_val <= 100:
                    st = int((start_val / 100.0) * total_null_bp)
                    length = int(((end_val - start_val) / 100.0) * total_null_bp)
                    pct_intervals.append((float(start_val), float(end_val)))
                else:
                    st = start_val
                    length = end_val - start_val
                    pct_st = (start_val / float(total_null_bp)) * 100.0
                    pct_end = (end_val / float(total_null_bp)) * 100.0
                    pct_intervals.append((pct_st, pct_end))
                parsed_intervals.append((st, length))

            total_alt_bp = sum(length for _, length in parsed_intervals)

            is_custom_np = args.np.lower() != "n1n2n3n4n5n6n7n8n9n10n11n12"
            is_custom_ap = args.ap.lower() != "a1a2a3"

            if args.output is None:
                interval_labels = []
                if intervals:
                    for st_v, end_v in intervals:
                        st_s = f"{int(st_v)}" if int(st_v) == st_v else f"{st_v}"
                        end_s = f"{int(end_v)}" if int(end_v) == end_v else f"{end_v}"
                        interval_labels.append(f"{st_s}-{end_s}")
                interval_str_label = ",".join(interval_labels) if interval_labels else ""

                if is_custom_np or is_custom_ap:
                    pattern_stem = build_interleaved_pattern_string(args.np, args.ap, pct_intervals if intervals else None)
                    if interval_str_label:
                        output_stem = f"{pattern_stem};{interval_str_label}"
                    else:
                        output_stem = pattern_stem
                else:
                    output_stem = interval_str_label if interval_str_label else exp_folder_name
                output_name = f"{output_stem}.fa"
            else:
                output_name = args.output.name
            output_path = concat_dir / output_name

            if output_path.exists() and output_path.stat().st_size > 0:
                print(f"[Concat Sequences] Concatenated FASTA already exists: {output_path}", file=sys.stderr)
                return

            print(f"Applying {len(intervals)} interval insertions (total length: {total_alt_bp} bp) to {args.np} block...", file=sys.stderr)
            with open(output_path, "w") as f:
                for header in base_headers:
                    f.write(f">{header}\n")
                    null_seq = null_fastas[header]
                    alt_seq = alt_fastas.get(header, null_seq)
                    
                    alt_sub = alt_seq[:total_alt_bp] if len(alt_seq) >= total_alt_bp else (alt_seq * (total_alt_bp // max(1, len(alt_seq)) + 1))[:total_alt_bp]
                    trimmed_null_length = max(0, len(null_seq) - total_alt_bp)

                    final_seq = []
                    null_idx = 0
                    alt_idx = 0
                    
                    sorted_intervals = sorted(parsed_intervals, key=lambda x: x[0])
                    curr_pos = 0
                    
                    for st, length in sorted_intervals:
                        null_len = st - curr_pos
                        if null_len > 0:
                            final_seq.append(null_seq[null_idx:null_idx + null_len])
                            null_idx += null_len
                            curr_pos += null_len
                        
                        alt_chunk = alt_sub[alt_idx:alt_idx + length]
                        final_seq.append(alt_chunk)
                        alt_idx += length
                        curr_pos += length
                    
                    remaining_null_len = trimmed_null_length - null_idx
                    if remaining_null_len > 0:
                        final_seq.append(null_seq[null_idx:null_idx + remaining_null_len])

                    f.write(f"{''.join(final_seq)}\n")

            print(f"Writing output FASTA: {output_path}", file=sys.stderr)
            print("Done.", file=sys.stderr)
            return
        else:
            is_pure_null = (not alt_path or str(alt_path).lower() == "null")
            np_str = "" if args.np == "n1n2n3n4n5n6n7n8n9n10n11n12" else args.np
            ap_str = "" if args.ap == "a1a2a3" else args.ap

            if args.output is None:
                parts = [np_str] if is_pure_null else [np_str, ap_str]
                parts = [p for p in parts if p]
                output_name = f"{'_'.join(parts)}.fa" if parts else f"{exp_folder_name}.fa"
            else:
                output_name = args.output.name
            output_path = concat_dir / output_name

            if output_path.exists() and output_path.stat().st_size > 0:
                print(f"[Concat Sequences] Concatenated FASTA already exists: {output_path}", file=sys.stderr)
                return
            
            with open(output_path, "w") as f:
                for header in base_headers:
                    f.write(f">{header}\n")
                    if is_pure_null:
                        f.write(f"{null_fastas[header]}\n")
                    else:
                        null_seq = null_fastas[header]
                        alt_seq = alt_fastas.get(header, "")
                        
                        total_alt_bp = len(alt_seq)
                        trimmed_null_length = max(0, len(null_seq) - total_alt_bp)
                        base_seq = list(null_seq[:trimmed_null_length])
                        
                        # Default insertion: prefix alt sequence
                        base_seq[0:0] = list(alt_seq)
                        
                        f.write(f"{''.join(base_seq[:len(null_seq)])}\n")
            print(f"Writing base FASTA: {output_path}", file=sys.stderr)
            print("Done.", file=sys.stderr)
            return

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

        exp_rel_path = get_experiment_relative_path(first_alt if first_alt else first_null)
        exp_folder_name = exp_rel_path.name
        exp_dir = output_dir / exp_rel_path
        concat_dir = exp_dir / "concat"
        concat_dir.mkdir(parents=True, exist_ok=True)

        if args.output is None:
            output_name = f"{pattern_str}.fa" if pattern_str else f"{exp_folder_name}.fa"
            output_path = concat_dir / output_name
        elif not args.output.is_absolute():
            output_path = concat_dir / args.output.name
        else:
            output_path = args.output

        if output_path.exists() and output_path.stat().st_size > 0:
            print(f"[Concat Sequences] Concatenated FASTA already exists: {output_path}", file=sys.stderr)
            return

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

