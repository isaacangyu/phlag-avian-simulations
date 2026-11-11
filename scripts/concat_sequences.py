from __future__ import annotations

import os
import re
import sys
import argparse
from typing import Union
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPTS_DIR.parent

_MISSING_ALT_CHUNK_WARNINGS: list[str] = []


def get_missing_alt_chunk_warnings() -> list[str]:
    """Warnings recorded by concatenate_pair() for alt chunks that were requested but unavailable."""
    return list(_MISSING_ALT_CHUNK_WARNINGS)


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
    for ext in (".ts", ".fa", ".fasta", ".gtrees"):
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
      - admixture/low/<exp_name> (time <= 4,000,000) or admixture/high/<exp_name> (time > 4,000,000)
      - null -> null
    Classification comes solely from parse_filename()'s `change`/`time` fields (never from a raw
    substring search over the path string), and there is no flat fallback: an experiment whose
    category can't be determined is an error, not a guess.
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

    change = (meta.get("change") or "").strip()

    if change in ("10X_up", "10X_down"):
        sub = "up" if change == "10X_up" else "down"
        return Path("10X") / sub / exp_name

    if change in ("recombination_up", "recombination_down"):
        sub = "up" if change == "recombination_up" else "down"
        return Path("recombination") / sub / exp_name

    if change == "admixture":
        time_str = meta.get("time") or ""
        m_num = re.search(r"\d+", time_str)
        if not m_num:
            sys.exit(f"Error: Could not determine admixture divergence time for '{path_or_stem}' - cannot classify as low/high.")
        div_time = float(m_num.group())
        sub = "high" if div_time > 4_000_000 else "low"
        return Path("admixture") / sub / exp_name

    sys.exit(f"Error: Could not determine experiment category (10X/recombination/admixture) for '{path_or_stem}' (change={change!r}).")


GTREES_EXT = ".gtrees"

# Pre-made gene tree files use their own change names, e.g. 'popsize_decrease_10x-N109-500Kb-1'
# or 'recombination_suppression-N115-500Kb-1'. canonical_gtrees_stem() rewrites them into the .ts
# naming parse_filename() already understands, so classification still goes through
# get_experiment_relative_path() and nothing else.
GTREES_POPSIZE_DIRECTIONS = {"increase": "up", "up": "up", "decrease": "down", "down": "down", "suppression": "down"}
_GTREES_ALT_RE = re.compile(r"^(popsize|recombination)_(increase|decrease|suppression|up|down)(?:_10x)?-([A-Za-z0-9_]+)(.*)$", re.IGNORECASE)
# Admixture gene tree names carry no rate; every admixture run used rate 0.90.
GTREES_ADMIXTURE_RATE = "rate090"
_GTREES_ADMIXTURE_RE = re.compile(r"^(?:simulations-)?admixture-([A-Za-z0-9_]+)-(\d+)(.*)$", re.IGNORECASE)
_GTREES_NULL_RE = re.compile(r"^(?:default|null)(?:-neoaves)?(?:-500[Kk][Bb])?(-\d+)?$", re.IGNORECASE)
# Pre-joined files like 'popsize_decrease_10x-concat-N109' are just chunks 1-3 back to back.
_GTREES_JOINED_RE = re.compile(r"^(popsize|recombination)_[A-Za-z0-9_]+-concat-[A-Za-z0-9_]+$", re.IGNORECASE)


def is_joined_gtrees_file(path_or_name: str | Path) -> bool:
    """True for a pre-joined gene tree file ('<change>-concat-<branch>.gtrees'), which has no chunk index."""
    return bool(_GTREES_JOINED_RE.match(get_base_name(path_or_name)))


def canonical_gtrees_stem(path_or_name: str | Path) -> str:
    """
    Maps a gene tree file name onto the equivalent .ts-style stem:
      'popsize_increase_10x-Accipitriformes-500Kb-1' -> 'anomaly_case-Accipitriformes_10X_up-500Kb-1'
      'popsize_decrease_10x-N109-500Kb-1'            -> 'anomaly_case-N109_10X_down-500Kb-1'
      'recombination_increase_10x-N113-500Kb-1'      -> 'anomaly_case-recombination_increase_N113-500Kb-1'
      'recombination_suppression-N115-500Kb-1'       -> 'anomaly_case-recombination_suppression_N115-500Kb-1'
      'simulations-admixture-N90_Columbinapicui-2099554-500Kb-1'
                                                     -> 'N90_Columbinapicui_admixture_rate090-time2099554-500Kb-1'
      'default-500Kb-3' / 'null-3'                   -> 'null-3'
    Any other name is returned unchanged (minus its extension).
    """
    stem = get_base_name(path_or_name)
    m = _GTREES_ALT_RE.match(stem)
    if m:
        kind, direction, branch, tail = m.group(1).lower(), m.group(2).lower(), m.group(3), m.group(4)
        if kind == "popsize":
            return f"anomaly_case-{branch}_10X_{GTREES_POPSIZE_DIRECTIONS[direction]}{tail}"
        return f"anomaly_case-recombination_{direction}_{branch}{tail}"
    m = _GTREES_ADMIXTURE_RE.match(stem)
    if m:
        branch, time, tail = m.groups()
        return f"{branch}_admixture_{GTREES_ADMIXTURE_RATE}-time{time}{tail}"
    m = _GTREES_NULL_RE.match(stem)
    if m:
        return f"null{m.group(1) or ''}"
    return stem


def classify_gtrees_file(path: str | Path) -> tuple[Path, str]:
    """Returns (experiment relative path, chunk index) for a gene tree file, e.g.
    'popsize_decrease_10x-N109-500Kb-1.gtrees' -> (Path('10X/down/N109'), '1')."""
    p = Path(path)
    stem = canonical_gtrees_stem(p)
    return get_experiment_relative_path(p.parent / stem), get_chunk_index(stem)


def resolve_fasta_file(f: Path, output_dir: Path) -> Path:
    """Returns the single deterministic simulated FASTA path for a base name/chunk, erroring if missing."""
    bname = get_base_name(f)
    exp_rel = get_experiment_relative_path(f)
    chunk_idx = get_chunk_index(bname)
    name = f"{chunk_idx}.fa" if chunk_idx else f"{bname}.fa"

    candidate = output_dir / exp_rel / "simulated" / name
    if candidate.exists() and candidate.is_file():
        return candidate

    sys.exit(f"Error: Sequence FASTA file does not exist: {candidate}")


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
    Finds existing numeric chunk indices for a base name under its single deterministic
    experiment directory (ts/, simulated/, concat/) - no legacy-layout, CWD, or external
    directory scanning.
    """
    exp_rel = get_experiment_relative_path(base_path)
    found_indices = set()

    for subdir in ("ts", "simulated", "concat"):
        sdir = output_dir / exp_rel / subdir
        if not sdir.is_dir():
            continue
        for p in sdir.iterdir():
            stem = get_base_name(p)
            if stem.isdigit():
                found_indices.add(int(stem))

    return sorted(found_indices)


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


def interval_label(interval_str: str) -> str:
    """
    Computes the auto-generated concat filename stem for a given interval string, using the
    exact same parsing/formatting `concatenate_pair` uses when np/ap are left at their defaults.
    Single source of truth so callers (e.g. a structure audit) can never drift from what
    `concatenate_pair` actually names its output file.
    """
    intervals = parse_interval_pattern(interval_str)
    if not intervals:
        sys.exit(f"Error: Invalid interval string: {interval_str}")
    labels = []
    for st_v, end_v in intervals:
        st_s = f"{int(st_v)}" if int(st_v) == st_v else f"{st_v}"
        end_s = f"{int(end_v)}" if int(end_v) == end_v else f"{end_v}"
        labels.append(f"{st_s}-{end_s}")
    return ",".join(labels)


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


def concat_output_stem(
    alt_path: Path | None,
    exp_folder_name: str,
    np_pattern: str,
    ap_pattern: str,
    interval: str | None,
    pct_intervals: list[tuple[float, float]] | None,
) -> str:
    """
    The concat file stem (no extension) shared by concatenate_pair() (.fa) and
    concatenate_gtrees_pair() (.gtrees), so both always name the same interval the same way.
    pct_intervals is only used for custom np/ap patterns.
    """
    if interval:
        is_custom_np = np_pattern.lower() != "n1n2n3n4n5n6n7n8n9n10n11n12"
        is_custom_ap = ap_pattern.lower() != "a1a2a3"
        interval_str_label = interval_label(interval)
        if is_custom_np or is_custom_ap:
            pattern_stem = build_interleaved_pattern_string(np_pattern, ap_pattern, pct_intervals)
            return f"{pattern_stem};{interval_str_label}" if interval_str_label else pattern_stem
        return interval_str_label if interval_str_label else exp_folder_name

    is_pure_null = (not alt_path or str(alt_path).lower() == "null")
    np_str = "" if np_pattern == "n1n2n3n4n5n6n7n8n9n10n11n12" else np_pattern
    ap_str = "" if ap_pattern == "a1a2a3" else ap_pattern
    parts = [np_str] if is_pure_null else [np_str, ap_str]
    parts = [p for p in parts if p]
    return "_".join(parts) if parts else exp_folder_name


DEFAULT_NP_PATTERN = "n1n2n3n4n5n6n7n8n9n10n11n12"
NULL_ORDER_EXT = ".null_order"


def shuffled_null_order(exp_rel_path: Path, interval: str | None, ap_pattern: str) -> list[int]:
    """
    A random order of null chunks 1..12 for one concat file. Each file gets its own shuffle,
    seeded from its experiment path, interval and alt pattern, so a rebuild gives the same order.
    """
    import hashlib
    key = f"{exp_rel_path.as_posix()}|{interval or ''}|{ap_pattern.lower()}"
    seed = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big")
    order = list(range(1, 13))
    random.Random(seed).shuffle(order)
    return order


_NULL_CHUNK_CACHE: dict[Path, dict[str, str]] = {}


def _read_null_chunk(path: Path) -> dict[str, str]:
    """Reads a null chunk FASTA once per process (the same 12 chunks are reused by every concat)."""
    if path not in _NULL_CHUNK_CACHE:
        _NULL_CHUNK_CACHE[path] = read_fasta(path)
    return _NULL_CHUNK_CACHE[path]


def concatenate_pair(
    null_path: Path | None,
    alt_path: Path | None,
    output_dir: Path,
    np_pattern: str = "n1n2n3n4n5n6n7n8n9n10n11n12",
    ap_pattern: str = "a1a2a3",
    proportion: str | float | None = None,
    interval: str | None = None,
) -> Path:
    """
    Concatenates a null sequence with an alt sequence into a 6Mbp FASTA, written to
    output_dir/<exp_rel>/concat/<auto-generated-name>.fa. Called directly (no argv/CLI
    round-trip) by both this module's own `main()` and by `simulate.py`.

    With the default np pattern, the 12 null chunks are joined in a random order picked for this
    file (see shuffled_null_order); the order is saved next to it as <name>.null_order. A custom
    np pattern is used as given. Null-only concat files are not written: the null chunks are only
    joined in memory while building an alt concat.
    """
    if not alt_path or str(alt_path).lower() == "null":
        sys.exit("Error: Null-only concat files are no longer written. Give an alt experiment to concatenate against null.")

    intervals = parse_interval_pattern(interval) if interval else None

    exp_rel_path = get_experiment_relative_path(alt_path if alt_path and str(alt_path).lower() != "null" else null_path)
    exp_folder_name = exp_rel_path.name
    exp_dir = output_dir / exp_rel_path
    concat_dir = exp_dir / "concat"

    def _compute_output_path(pct_intervals):
        stem = concat_output_stem(alt_path, exp_folder_name, np_pattern, ap_pattern, interval, pct_intervals)
        return concat_dir / f"{stem}.fa"

    # The output filename only depends on sequence content when an interval is given in raw
    # base-pair form (>100) - every other case (including everything simulate.py's default sweep
    # produces, which is always a percentage interval) can be named without touching any FASTA
    # data, so an already-existing target can be detected before paying for a multi-hundred-MB read.
    needs_bp_length_for_naming = bool(intervals) and any(end_val > 100 for _, end_val in intervals)
    output_path = None
    if not needs_bp_length_for_naming:
        output_path = _compute_output_path([(float(s), float(e)) for s, e in intervals] if intervals else None)
        if output_path.exists() and output_path.stat().st_size > 0:
            print(f"[Concat Sequences] Concatenated FASTA already exists: {output_path}", file=sys.stderr)
            return output_path

    concat_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load the null base sequence
    null_chunk_indices = None
    if null_path and str(null_path).lower() != "null":
        base_null_fa = resolve_fasta_file(null_path, output_dir)
        null_fastas = read_fasta(base_null_fa)
    else:
        if np_pattern.lower() == DEFAULT_NP_PATTERN:
            null_chunk_indices = shuffled_null_order(exp_rel_path, interval, ap_pattern)
        else:
            matches = re.findall(r'[an](\d+)', np_pattern.lower())
            null_chunk_indices = [int(m) for m in matches] if matches else [i for i in range(1, 13)]
        print(f"[Concat Sequences] Null chunk order: {' '.join(map(str, null_chunk_indices))}", file=sys.stderr)

        chunk_paths = [output_dir / "null" / "simulated" / f"{idx}.fa" for idx in null_chunk_indices]
        for cf in chunk_paths:
            if not cf.exists():
                sys.exit(f"Error: Null chunk FASTA does not exist: {cf}")
        chunks = [_read_null_chunk(cf) for cf in chunk_paths]
        null_fastas = {header: "".join(c[header] for c in chunks) for header in chunks[0]}

    base_headers = list(null_fastas.keys())

    if needs_bp_length_for_naming:
        total_null_bp_for_naming = sum(len(seq) for seq in null_fastas.values()) // len(base_headers) if base_headers else 6000000
        pct_intervals = [
            (float(s), float(e)) if e <= 100 else (s / float(total_null_bp_for_naming) * 100.0, e / float(total_null_bp_for_naming) * 100.0)
            for s, e in intervals
        ]
        output_path = _compute_output_path(pct_intervals)
        if output_path.exists() and output_path.stat().st_size > 0:
            print(f"[Concat Sequences] Concatenated FASTA already exists: {output_path}", file=sys.stderr)
            return output_path

    # 2. Load the alt chunks. Not every experiment has data for every requested chunk index
    # (e.g. only 1 tree sequence was ever simulated for it) - that's expected, not fatal: load
    # whichever requested chunks exist, log the ones that don't, and only error if none do.
    alt_fastas = {}
    if alt_path and str(alt_path).lower() != "null":
        matches = re.findall(r'[an](\d+)', ap_pattern.lower())
        alt_chunk_indices = [int(m) for m in matches] if matches else [1, 2, 3]

        loaded_count = 0
        for idx in alt_chunk_indices:
            cf = construct_indexed_path(alt_path, str(idx))
            try:
                rf = resolve_fasta_file(cf, output_dir)
            except SystemExit as e:
                warning = f"Alt chunk {idx} unavailable for '{alt_path}': {e}"
                print(f"[Concat Sequences] Warning: {warning}", file=sys.stderr)
                _MISSING_ALT_CHUNK_WARNINGS.append(warning)
                continue
            chunk_data = read_fasta(rf)
            loaded_count += 1
            for header in base_headers:
                alt_fastas.setdefault(header, "")
                alt_fastas[header] += chunk_data.get(header, "")

        if loaded_count == 0:
            sys.exit(f"Error: None of the requested alt chunks {alt_chunk_indices} exist for '{alt_path}'.")

    if intervals:
        total_null_bp = sum(len(seq) for seq in null_fastas.values()) // len(base_headers) if base_headers else 6000000

        parsed_intervals = []
        for start_val, end_val in intervals:
            if end_val <= 100:
                st = int((start_val / 100.0) * total_null_bp)
                length = int(((end_val - start_val) / 100.0) * total_null_bp)
            else:
                st = start_val
                length = end_val - start_val
            parsed_intervals.append((st, length))

        total_alt_bp = sum(length for _, length in parsed_intervals)

        print(f"Applying {len(intervals)} interval insertions (total length: {total_alt_bp} bp) to {np_pattern} block...", file=sys.stderr)
        tmp_path = output_path.with_name(output_path.name + ".tmp")
        with open(tmp_path, "w") as f:
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

        _finish_concat_write(tmp_path, output_path, null_chunk_indices)
        print(f"Writing output FASTA: {output_path}", file=sys.stderr)
        print("Done.", file=sys.stderr)
        return output_path

    tmp_path = output_path.with_name(output_path.name + ".tmp")
    with open(tmp_path, "w") as f:
        for header in base_headers:
            f.write(f">{header}\n")
            null_seq = null_fastas[header]
            alt_seq = alt_fastas.get(header, "")
            # Default insertion: prefix alt sequence, trim the null tail to keep the length.
            f.write(f"{(alt_seq + null_seq[:max(0, len(null_seq) - len(alt_seq))])[:len(null_seq)]}\n")
    _finish_concat_write(tmp_path, output_path, null_chunk_indices)
    print(f"Writing base FASTA: {output_path}", file=sys.stderr)
    print("Done.", file=sys.stderr)
    return output_path


def _finish_concat_write(tmp_path: Path, output_path: Path, null_chunk_indices: list[int] | None) -> None:
    """Moves a finished concat into place and saves the null chunk order used next to it."""
    if null_chunk_indices is not None:
        output_path.with_suffix(NULL_ORDER_EXT).write_text("".join(f"n{i}" for i in null_chunk_indices) + "\n")
    tmp_path.replace(output_path)


GTREES_CHUNK_BP = 500_000


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


def prune_outgroups(newick: str) -> str:
    """
    Removes OUTGROUP_TAXA leaves from one Newick tree so every gene tree has the same Neoaves
    taxon set as the Neoaves-filtered FASTAs. Trees without outgroups are returned unchanged.
    """
    if not any(taxon in newick for taxon in OUTGROUP_TAXA):
        return newick

    import warnings
    import treeswift

    tree = treeswift.read_tree_newick(newick)
    outgroups = {nd.label for nd in tree.traverse_leaves()} & OUTGROUP_TAXA
    with warnings.catch_warnings():
        # Removing the outgroups leaves unifurcations whose support labels treeswift warns about.
        warnings.simplefilter("ignore")
        pruned = tree.extract_tree_without(outgroups)
        pruned.suppress_unifurcations()
    pruned.root.edge_length = None
    pruned.root.label = None
    return pruned.newick().removeprefix("[&R] ")


def read_gtrees(filepath: Path) -> list[str]:
    """Reads a gene tree file: one Newick tree per non-empty line, in genome order."""
    trees = []
    with open(filepath, "r") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            if not line.endswith(";"):
                sys.exit(f"Error: Line {line_no} of {filepath} is not a complete Newick tree (no trailing ';').")
            trees.append(line)
    if not trees:
        sys.exit(f"Error: Gene tree file is empty: {filepath}")
    return trees


def _read_gtrees_chunks(chunk_dir: Path, indices: list[int], label: str) -> list[str]:
    """Reads <chunk_dir>/<idx>.gtrees for every idx in order (outgroups pruned), erroring if any is missing."""
    missing = [chunk_dir / f"{idx}{GTREES_EXT}" for idx in indices if not (chunk_dir / f"{idx}{GTREES_EXT}").is_file()]
    if missing:
        sys.exit(f"Error: {label} gene tree chunk file(s) missing: {', '.join(str(p) for p in missing)}")
    trees = []
    for idx in indices:
        trees.extend(_read_pruned_gtrees(chunk_dir / f"{idx}{GTREES_EXT}"))
    return trees


_GTREES_CHUNK_CACHE: dict[Path, list[str]] = {}


def _read_pruned_gtrees(path: Path) -> list[str]:
    """One gene tree chunk with outgroups pruned, read once per process."""
    if path not in _GTREES_CHUNK_CACHE:
        _GTREES_CHUNK_CACHE[path] = [prune_outgroups(t) for t in read_gtrees(path)]
    return _GTREES_CHUNK_CACHE[path]


def concatenate_gtrees_pair(
    null_path: Path | None,
    alt_path: Path | None,
    output_dir: Path,
    np_pattern: str = "n1n2n3n4n5n6n7n8n9n10n11n12",
    ap_pattern: str = "a1a2a3",
    proportion: str | float | None = None,
    interval: str | None = None,
) -> Path:
    """
    Gene tree counterpart of concatenate_pair(). Builds output_dir/<exp_rel>/concat/<stem>.gtrees
    from null/gtrees/<n>.gtrees (np order) with <exp_rel>/gtrees/<a>.gtrees (ap order) inserted,
    using the same insert-and-shift rule and the same output stem as the MSA concat, so tree i of
    concat/<stem>.gtrees lines up with window i of concat/<stem>.fa:
      out = null[0:st] + alt[0:len] + null[st:N - len]   (N = total null trees)
    Percent intervals use int(start/100*N) and int((end-start)/100*N); base-pair intervals (>100)
    assume each chunk is GTREES_CHUNK_BP long. With the default np pattern, the null chunks use
    the same random order as the matching concat/<stem>.fa (checked against its <stem>.null_order
    file). Null-only concats are not written. `proportion` is accepted for parity and, like in
    concatenate_pair, does not change the output. Skips if the target already exists.
    """
    if null_path and str(null_path).lower() != "null":
        sys.exit(f"Error: Gene tree concat only supports the chunked null base ('null'), got '{null_path}'.")
    if not alt_path or str(alt_path).lower() == "null":
        sys.exit("Error: Null-only concat files are no longer written. Give an alt experiment to concatenate against null.")

    intervals = parse_interval_pattern(interval) if interval else None

    exp_source = Path(alt_path)
    # Classify an experiment dir by its path under output_dir only, so ancestor dir names above
    # output_dir can never leak into parse_filename()'s branch guess.
    if exp_source.is_absolute() and exp_source.is_relative_to(output_dir):
        exp_source = exp_source.relative_to(output_dir)
    if exp_source.name.lower().endswith(GTREES_EXT):
        exp_source = exp_source.parent / canonical_gtrees_stem(exp_source)
    exp_rel_path = get_experiment_relative_path(exp_source)
    concat_dir = output_dir / exp_rel_path / "concat"

    if np_pattern.lower() == DEFAULT_NP_PATTERN:
        null_chunk_indices = shuffled_null_order(exp_rel_path, interval, ap_pattern)
    else:
        null_matches = re.findall(r'[an](\d+)', np_pattern.lower())
        null_chunk_indices = [int(m) for m in null_matches] if null_matches else list(range(1, 13))
    total_bp = GTREES_CHUNK_BP * len(null_chunk_indices)

    pct_intervals = None
    if intervals:
        pct_intervals = [
            (float(s), float(e)) if e <= 100 else (s / float(total_bp) * 100.0, e / float(total_bp) * 100.0)
            for s, e in intervals
        ]
    stem = concat_output_stem(alt_path, exp_rel_path.name, np_pattern, ap_pattern, interval, pct_intervals)
    output_path = concat_dir / f"{stem}{GTREES_EXT}"
    if output_path.exists() and output_path.stat().st_size > 0:
        print(f"[Concat Gene Trees] Concatenated gene trees already exist: {output_path}", file=sys.stderr)
        return output_path

    # The matching .fa must have used the same null order, or the trees won't line up with it.
    order_path = output_path.with_suffix(NULL_ORDER_EXT)
    order_str = "".join(f"n{i}" for i in null_chunk_indices)
    if order_path.exists() and order_path.read_text().strip() != order_str:
        sys.exit(f"Error: Null order {order_str} does not match {order_path} ({order_path.read_text().strip()}).")
    print(f"[Concat Gene Trees] Null chunk order: {' '.join(map(str, null_chunk_indices))}", file=sys.stderr)

    null_trees = _read_gtrees_chunks(output_dir / "null" / "gtrees", null_chunk_indices, "Null")
    num_trees = len(null_trees)

    alt_matches = re.findall(r'[an](\d+)', ap_pattern.lower())
    alt_chunk_indices = [int(m) for m in alt_matches] if alt_matches else [1, 2, 3]
    alt_trees = _read_gtrees_chunks(output_dir / exp_rel_path / "gtrees", alt_chunk_indices, f"Alt ({exp_rel_path})")

    if intervals:
        parsed_intervals = []
        for start_val, end_val in intervals:
            if end_val <= 100:
                st = int((start_val / 100.0) * num_trees)
                length = int(((end_val - start_val) / 100.0) * num_trees)
            else:
                st = int(start_val / float(total_bp) * num_trees)
                length = int((end_val - start_val) / float(total_bp) * num_trees)
            parsed_intervals.append((st, length))

        total_alt = sum(length for _, length in parsed_intervals)
        if len(alt_trees) < total_alt:
            sys.exit(f"Error: Interval '{interval}' needs {total_alt} alt gene trees but only {len(alt_trees)} exist for '{exp_rel_path}' (chunks {ap_pattern}).")
        alt_sub = alt_trees[:total_alt]
        trimmed_null_length = max(0, num_trees - total_alt)

        final_trees = []
        null_idx = 0
        alt_idx = 0
        curr_pos = 0
        for st, length in sorted(parsed_intervals, key=lambda x: x[0]):
            null_len = st - curr_pos
            if null_len > 0:
                final_trees.extend(null_trees[null_idx:null_idx + null_len])
                null_idx += null_len
                curr_pos += null_len
            final_trees.extend(alt_sub[alt_idx:alt_idx + length])
            alt_idx += length
            curr_pos += length

        remaining_null_len = trimmed_null_length - null_idx
        if remaining_null_len > 0:
            final_trees.extend(null_trees[null_idx:null_idx + remaining_null_len])
    else:
        # Default insertion (no interval): prefix the alt trees, like concatenate_pair.
        final_trees = (alt_trees + null_trees[:max(0, num_trees - len(alt_trees))])[:num_trees]

    concat_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(output_path.name + ".tmp")
    with open(tmp_path, "w") as f:
        for tree in final_trees:
            f.write(f"{tree}\n")
    if not order_path.exists():
        order_path.write_text(order_str + "\n")
    tmp_path.replace(output_path)
    print(f"[Concat Gene Trees] Wrote {len(final_trees)} gene trees: {output_path}", file=sys.stderr)
    return output_path


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
        return concatenate_pair(
            null_path, alt_path, output_dir,
            np_pattern=args.np, ap_pattern=args.ap,
            proportion=args.proportion, interval=args.interval,
        )

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

        output_name = f"{pattern_str}.fa" if pattern_str else f"{exp_folder_name}.fa"
        output_path = concat_dir / output_name

        if output_path.exists() and output_path.stat().st_size > 0:
            print(f"[Concat Sequences] Concatenated FASTA already exists: {output_path}", file=sys.stderr)
            return output_path

        print(f"Concatenating 12 sequence files (6Mbp total)...", file=sys.stderr)
        with open(output_path, "w") as f:
            for header in base_headers:
                f.write(f">{header}\n")
                seq_concat = "".join(fa[header] for fa in fastas)
                f.write(f"{seq_concat}\n")

        print(f"Writing output FASTA: {output_path}", file=sys.stderr)

        print("Done.", file=sys.stderr)

        return output_path

    else:
        parser.print_help()
        sys.exit("\nError: Must specify (null, alt, [pattern_or_count]) or exactly 12 files to concatenate a 6Mbp sequence.")


if __name__ == "__main__":
    main()

