#!/usr/bin/env python3
"""
Simulate command pipeline. Four modes, selected by exactly one positional target:
  (no target)      Ensures every already-organized node experiment directory (10X/recombination/
                    admixture; null is skipped) is complete: renames any non-numerically-named
                    .ts files, simulates/filters every chunk, ensures mapping, and runs the
                    default six-proportion concat sweep. Never touches loose .ts files outside a
                    node directory - use -r for those.
  -r/--run-all      Organizes every loose tree sequence file, simulates/maps/concats each
                    experiment with the default six-proportion sweep. No other arguments.
  null              Simulates/filters/concats the null experiment only.
  <relative_dir>    Simulates/maps/concats an already-existing experiment directory (relative
                    to the simulations dir) against null, applying --np/--ap/-p/interval if given.
"""

from __future__ import annotations

import re
import sys
import argparse
from pathlib import Path

# Ensure parent directory of scripts is in sys.path
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

REPO_ROOT = SCRIPTS_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from simulate_sequences import simulate_ts_to_fasta, get_output_dir
from concat_sequences import (
    get_experiment_relative_path,
    get_chunk_index,
    concatenate_pair,
    get_missing_alt_chunk_warnings,
    parse_proportion,
    concatenate_gtrees_pair,
    OUTGROUP_TAXA,
)
from estimate_gene_trees import ensure_estimated_gene_trees, DEFAULT_SEED


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


def write_fasta(sequences: dict, filepath: Path):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n{seq}\n")


def get_admixture_source_clade(tree_path: Path, target_node: str, source_node: str) -> str:
    """
    Determines which quadripartition clade (P1, P2, P3, Po) the source node belongs to
    relative to the target node's branch in the species tree.
    Returns clade label in lowercase: '1', '2', '3', or 'o'.
    """
    import treeswift as ts
    tree = ts.read_tree_newick(tree_path)
    label_to_node = tree.label_to_node(selection="all")

    if target_node not in label_to_node:
        return "o"
    if source_node not in label_to_node:
        return "o"

    nd_target = label_to_node[target_node]
    nd_source = label_to_node[source_node]

    # Source leaf set (or leaf itself)
    source_leaves = set(nd.get_label() for nd in nd_source.traverse_leaves()) if not nd_source.is_leaf() else {nd_source.get_label()}

    children = nd_target.child_nodes()
    if len(children) >= 2:
        p1_leaves = set(nd.get_label() for nd in children[0].traverse_leaves())
        p2_leaves = set()
        for child in children[1:]:
            p2_leaves.update(nd.get_label() for nd in child.traverse_leaves())

        if source_leaves.issubset(p1_leaves):
            return "1"
        if source_leaves.issubset(p2_leaves):
            return "2"

    parent_node = nd_target.get_parent()
    if parent_node:
        p3_leaves = set()
        for child in parent_node.child_nodes():
            if child != nd_target:
                p3_leaves.update(nd.get_label() for nd in child.traverse_leaves())
        if source_leaves.issubset(p3_leaves):
            return "3"

    return "o"


def generate_mapping(tree_path: Path, node_label: str, output_path: Path):
    import treeswift as ts
    tree = ts.read_tree_newick(tree_path)
    label_to_node = tree.label_to_node(selection="all")
    
    target_node = node_label.split("_")[0] if "_" in node_label else node_label
    if target_node not in label_to_node and node_label in label_to_node:
        target_node = node_label
        
    if target_node not in label_to_node:
        raise ValueError(f"Error: Internal node '{target_node}' not found in species tree '{tree_path}'.")
        
    nd_target = label_to_node[target_node]
    
    # Quadripartition around the branch ABOVE target internal node (edge between parent_node and nd_target):
    # P1 & P2 = Leaves under children of target internal node
    children = nd_target.child_nodes()
    if len(children) < 2:
        raise ValueError(f"Error: Internal node '{target_node}' has fewer than 2 children in species tree.")
        
    p1_leaves = set(nd.get_label() for nd in tree.extract_subtree(children[0]).traverse_leaves())
    p2_leaves = set()
    for child in children[1:]:
        p2_leaves.update(nd.get_label() for nd in tree.extract_subtree(child).traverse_leaves())
    
    # P3 = Leaves under sister clade(s) of target internal node (attached to parent_node)
    parent_node = nd_target.get_parent()
    p3_leaves = set()
    if parent_node:
        for child in parent_node.child_nodes():
            if child != nd_target:
                p3_leaves.update(nd.get_label() for nd in tree.extract_subtree(child).traverse_leaves())
                
    # Po = All remaining leaves in species tree
    all_leaves = set(nd.get_label() for nd in tree.traverse_leaves())
    po_leaves = all_leaves - p1_leaves - p2_leaves - p3_leaves
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for leaf in sorted(p1_leaves):
            f.write(f"{leaf}\tP1\n")
        for leaf in sorted(p2_leaves):
            f.write(f"{leaf}\tP2\n")
        for leaf in sorted(p3_leaves):
            f.write(f"{leaf}\tP3\n")
        for leaf in sorted(po_leaves):
            f.write(f"{leaf}\tPo\n")
            
    print(f"Successfully generated mapping at '{output_path}':")
    print(f"  P1: {len(p1_leaves)} leaves")
    print(f"  P2: {len(p2_leaves)} leaves")
    print(f"  P3: {len(p3_leaves)} leaves")
    print(f"  Po: {len(po_leaves)} leaves")


def get_base_name(path: str | Path | None) -> str:
    """Extracts base name from a file path or string, stripping .ts, .fa, or .fasta extensions if present."""
    if path is None:
        return ""
    name = Path(path).name
    for ext in (".ts", ".fa", ".fasta"):
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return Path(path).stem


DEFAULT_CENTERED_PROPORTIONS = [0.02, 0.05, 0.1, 0.15, 0.20, 0.25]


def proportion_to_centered_interval(proportion: float) -> str:
    """Converts a proportion of total sequence length into a percentage interval string centered at the midpoint (50%)."""
    half_width = (proportion * 100) / 2
    start = 50 - half_width
    end = 50 + half_width
    return f"{start:g}-{end:g}"


KNOWN_CHANGES_MAP = {
    "recombination_increase": "recombination_up",
    "recombination_up": "recombination_up",
    "recombination_decrease": "recombination_down",
    "recombination_suppression": "recombination_down",
    "recombination_down": "recombination_down",
    "10x_up": "10X_up",
    "10X_up": "10X_up",
    "10x_down": "10X_down",
    "10X_down": "10X_down",
    "admixture": "admixture",
}


def parse_filename(path_or_str: str | Path) -> dict:
    """
    Parses simulation filename into structured metadata:
      - prefix: 'null' or 'alt'
      - is_neoaves: bool
      - branch: str | None
      - change: str | None ('recombination_up', '10X_up', 'admixture', etc.)
      - rate: str | None
      - time: str | None
      - index: str | None
      - ext: str
      - is_concat: bool
    """
    filename = Path(path_or_str).name
    ext = ""
    for e in (".ts", ".fa", ".fasta", ".tsv"):
        if filename.lower().endswith(e):
            ext = e
            filename = filename[:-len(e)]
            break
    if not ext:
        ext = Path(path_or_str).suffix

    is_neoaves = "neoaves" in filename.lower()
    is_concat = filename.startswith("concat_") or bool(re.search(r"[an]\d+[an]\d+", filename))

    prefix = "alt" if ("anomaly" in filename.lower() or filename.startswith("alt")) else "null"
    index = ""
    m_idx = re.search(r"-(\d+)$", filename)
    if m_idx:
        index = m_idx.group(1)
        filename = filename[:m_idx.start()]

    branch = None
    change = None
    rate = None
    time = None

    # Handle admixture: rate and time (with or without explicit 'admixture' keyword)
    m_adm = re.search(r"(?:admixture[-_])?([A-Za-z0-9_]+?)[-_](?:admixture[-_])?(rate\d+|0\d+)[-_](time\d+|\d+)", filename, re.IGNORECASE)
    if m_adm and m_adm.group(1).lower() != "null":
        branch = m_adm.group(1)
        change = "admixture"
        rate = m_adm.group(2)
        time = m_adm.group(3)

    if not change:
        # Check 10X_up / 10X_down
        m_10x = re.search(r"([A-Za-z0-9_]+?)[-_](10X_up|10X_down|10x_up|10x_down)", filename, re.IGNORECASE)
        if m_10x:
            branch = m_10x.group(1)
            change = m_10x.group(2).replace("10x_", "10X_")
        else:
            m_10x_rev = re.search(r"(10X_up|10X_down|10x_up|10x_down)[-_]([A-Za-z0-9_]+)", filename, re.IGNORECASE)
            if m_10x_rev:
                change = m_10x_rev.group(1).replace("10x_", "10X_")
                branch = m_10x_rev.group(2).split("-500K")[0].split("-500k")[0]

    if not change:
        # Check recombination_increase / recombination_suppression / recombination_decrease / recombination_up / recombination_down
        m_rec = re.search(r"recombination_(increase|decrease|suppression|up|down)[-_]([A-Za-z0-9_]+)", filename, re.IGNORECASE)
        if m_rec:
            c_raw = f"recombination_{m_rec.group(1).lower()}"
            change = KNOWN_CHANGES_MAP.get(c_raw, c_raw)
            branch = m_rec.group(2).split("-500K")[0].split("-500k")[0]
        else:
            m_rec_rev = re.search(r"([A-Za-z0-9_]+)[-_]recombination_(increase|decrease|suppression|up|down)", filename, re.IGNORECASE)
            if m_rec_rev:
                branch = m_rec_rev.group(1)
                c_raw = f"recombination_{m_rec_rev.group(2).lower()}"
                change = KNOWN_CHANGES_MAP.get(c_raw, c_raw)

    # Infer change and branch from parent directory structure if not in filename
    full_path_str = str(path_or_str).lower()
    if not change and not filename.startswith("null") and filename.lower() != "null":
        if "10x/up" in full_path_str or "10x_up" in full_path_str:
            change = "10X_up"
        elif "10x/down" in full_path_str or "10x_down" in full_path_str:
            change = "10X_down"
        elif "recombination/up" in full_path_str or "recombination_up" in full_path_str:
            change = "recombination_up"
        elif "recombination/down" in full_path_str or "recombination_down" in full_path_str:
            change = "recombination_down"
        elif "admixture" in full_path_str:
            change = "admixture"

    # Infer change and branch from parent directory path if missing or numeric
    full_path_str = str(path_or_str).lower()
    if not change and not filename.startswith("null") and filename.lower() != "null":
        if "10x/up" in full_path_str or "10x_up" in full_path_str:
            change = "10X_up"
        elif "10x/down" in full_path_str or "10x_down" in full_path_str:
            change = "10X_down"
        elif "recombination/up" in full_path_str or "recombination_up" in full_path_str:
            change = "recombination_up"
        elif "recombination/down" in full_path_str or "recombination_down" in full_path_str:
            change = "recombination_down"
        elif "admixture" in full_path_str:
            change = "admixture"

    if not branch or branch.isdigit():
        try:
            p_obj = Path(path_or_str)
            reserved_folder_names = {"ts", "simulated", "concat", "simulations", "store", "10x", "recombination", "admixture", "up", "down", "low", "high", "iang", "home", "drive2", "root", "none", ".", ".."}
            for part in reversed(p_obj.parts[:-1]):
                if not part or part in ("/", "\\"):
                    continue
                part_clean = part.strip("-_")
                if part_clean.lower() not in reserved_folder_names and not part_clean.isdigit():
                    m_part_adm = re.search(r"(?:admixture[-_])?([A-Za-z0-9_]+?)[-_](?:admixture[-_])?(rate\d+|0\d+)[-_](time\d+|\d+)", part, re.IGNORECASE)
                    if m_part_adm:
                        branch = m_part_adm.group(1)
                        if not change:
                            change = "admixture"
                        if not rate:
                            rate = m_part_adm.group(2)
                        if not time:
                            time = m_part_adm.group(3)
                        break

                    m_part_10x = re.search(r"([A-Za-z0-9_]+?)[-_](10X_up|10X_down|10x_up|10x_down)", part, re.IGNORECASE)
                    if m_part_10x:
                        branch = m_part_10x.group(1)
                        if not change:
                            change = m_part_10x.group(2).replace("10x_", "10X_")
                        break

                    m_part_rec = re.search(r"([A-Za-z0-9_]+)[-_]recombination_(increase|decrease|suppression|up|down)", part, re.IGNORECASE)
                    if m_part_rec:
                        branch = m_part_rec.group(1)
                        if not change:
                            c_raw = f"recombination_{m_part_rec.group(2).lower()}"
                            change = KNOWN_CHANGES_MAP.get(c_raw, c_raw)
                        break

                    clean_part = part.split("-500K")[0].split("-500k")[0].split("_rate")[0].strip("-_")
                    if clean_part and not clean_part.isdigit() and clean_part.lower() not in reserved_folder_names:
                        branch = clean_part
                        break
        except Exception:
            pass

    # If branch is still None/numeric and filename is a valid branch name
    reserved = {"ts", "simulated", "concat", "simulations", "store", "10x", "recombination", "admixture", "up", "down", "low", "high", "iang", "home", "drive2", "root", "none", "null", "default", "."}
    if (not branch or branch.isdigit()) and filename.lower() not in reserved and not filename.isdigit():
        branch = filename

    if branch and (branch.isdigit() or branch.lower() in reserved):
        branch = None

    # Clean up branch name
    if branch:
        branch = re.sub(r"-500[Kk][Bb].*$", "", branch)
        branch = re.sub(r"-\d+$", "", branch)
        branch = branch.strip("-_")

    if (change or branch) and not filename.startswith("null") and filename.lower() != "null":
        prefix = "alt"

    # Format rate and time
    if rate and not rate.startswith("rate"):
        rate = f"rate{rate}"
    if time and not time.startswith("time"):
        time = f"time{time}"

    return {
        "raw_stem": filename,
        "prefix": prefix,
        "is_neoaves": is_neoaves,
        "branch": branch,
        "change": change,
        "rate": rate,
        "time": time,
        "index": index,
        "ext": ext,
        "is_concat": is_concat,
    }


def extract_branch_name(path_or_str: str | Path) -> str | None:
    """
    Extracts the branch name from a simulation filename.
    Returns branch name string (e.g. 'Nyctibiidae', 'N249', 'Strigiformes_N297') or None.
    """
    meta = parse_filename(path_or_str)
    return meta["branch"]


def generate_new_filename(path_or_str: str | Path) -> str:
    """
    Converts filename to new standard format:
      - 10X / recombination: branch-index.ext (e.g. N554-1.ts)
      - admixture: branch_rate-time-index.ext (e.g. Strigiformes_N297_rate090-time7099554-1.ts)
      - null: null-index.ext
    """
    meta = parse_filename(path_or_str)
    ext = meta["ext"]
    raw = meta["raw_stem"]

    if meta["is_concat"]:
        new_name = raw
        new_name = new_name.replace("default-500Kb", "null")
        new_name = new_name.replace("default-neoaves-500Kb", "null")
        new_name = new_name.replace("default-neoaves", "null")
        new_name = new_name.replace("default", "null")
        new_name = new_name.replace("anomaly_case-", "")
        new_name = new_name.replace("alt-", "")
        new_name = new_name.replace("-500Kb", "")
        new_name = new_name.replace("-500kb", "")
        new_name = new_name.replace("recombination_increase", "")
        new_name = new_name.replace("recombination_decrease", "")
        new_name = new_name.replace("recombination_suppression", "")
        new_name = new_name.replace("recombination_up", "")
        new_name = new_name.replace("recombination_down", "")
        new_name = new_name.replace("10X_up", "")
        new_name = new_name.replace("10X_down", "")
        new_name = new_name.replace("10x_up", "")
        new_name = new_name.replace("10x_down", "")
        new_name = new_name.replace("admixture_", "")
        new_name = new_name.replace("_admixture", "")
        new_name = new_name.replace("rate090_time7099554", "rate090-time7099554")
        return f"{new_name}{ext}"

    prefix = meta["prefix"]
    idx_str = f"-{meta['index']}" if meta["index"] else ""

    if meta["prefix"] == "null" or raw.startswith("null") or raw.lower() == "null":
        return f"null{idx_str}{ext}"

    branch = meta["branch"]
    change = meta["change"]

    if change == "admixture" and meta["rate"] and meta["time"]:
        rate_str = meta["rate"] if meta["rate"].startswith("rate") else f"rate{meta['rate']}"
        time_str = meta["time"] if meta["time"].startswith("time") else f"time{meta['time']}"
        return f"{branch}_{rate_str}-{time_str}{idx_str}{ext}"
    elif branch:
        return f"{branch}{idx_str}{ext}"
    else:
        return f"null{idx_str}{ext}"


def rename_all_files(target_dir: str | Path, dry_run: bool = False) -> list[tuple[Path, Path]]:
    """Batch renames files in target_dir according to standard scheme."""
    target_path = Path(target_dir)
    if not target_path.exists():
        return []

    renamed = []
    files = [p for p in target_path.rglob("*") if p.is_file() and not p.name.startswith(".")]

    for old_path in sorted(files):
        if old_path.name == "ape.fa":
            continue
        new_filename = generate_new_filename(old_path)
        new_path = old_path.parent / new_filename
        
        if old_path != new_path:
            if not dry_run:
                if not new_path.exists():
                    old_path.rename(new_path)
            renamed.append((old_path, new_path))

    return renamed


def extract_branch_name_from_filename(path: str | Path | None) -> str | None:
    """Extracts branch name from filename. Retains full compound branch names (e.g. N491_N482 or Strigiformes_N297)."""
    if not path:
        return None
    extracted = extract_branch_name(path)

    if not extracted:
        clean_name = get_base_name(path)
        m = re.search(r"[-_]([A-Za-z0-9_]+)_", clean_name)
        if m:
            extracted = m.group(1)

    return extracted


def ensure_mapping(node_label: str, exp_rel: Path, output_dir: Path) -> Path:
    """
    Checks if the mapping file exists at its single deterministic location
    (output_dir/exp_rel/<tree>_<label>_mapping.tsv), otherwise generates it. Tries the Neoaves
    species tree first, then the full Avian species tree; errors if the internal node is in
    neither.
    """
    import treeswift as ts

    target_node = node_label.split("_")[0] if "_" in node_label else node_label

    tree_candidates = [
        (REPO_ROOT / "speciestree" / "main-neoaves-numgen.nwk", "neoaves"),
        (REPO_ROOT / "speciestree" / "main-avian-numgen.nwk", "avian"),
    ]

    tree_path = None
    tree_stem = None
    for candidate_path, candidate_stem in tree_candidates:
        if not candidate_path.exists():
            continue
        tree = ts.read_tree_newick(candidate_path)
        if target_node in tree.label_to_node(selection="all"):
            tree_path = candidate_path
            tree_stem = candidate_stem
            break

    if tree_path is None:
        sys.exit(
            f"Error: Internal node '{target_node}' not found in either species tree "
            f"(main-neoaves-numgen.nwk or main-avian-numgen.nwk)."
        )

    formatted_label = node_label
    if "_" in node_label:
        parts = node_label.split("_")
        source_node = parts[1]
        clade = get_admixture_source_clade(tree_path, target_node, source_node)
        formatted_label = f"{target_node}:{clade}"

    target_name = f"{tree_stem}_{formatted_label}_mapping.tsv"
    mapping_dir = output_dir / exp_rel
    mapping_path = mapping_dir / target_name

    if mapping_path.exists():
        print(f"[Mapping] Mapping file already exists: {mapping_path}", file=sys.stderr)
        return mapping_path

    mapping_dir.mkdir(parents=True, exist_ok=True)
    print(f"[Mapping] Generating mapping for node='{node_label}' (tree={tree_stem}) -> {mapping_path}...", file=sys.stderr)
    generate_mapping(tree_path, node_label, mapping_path)
    return mapping_path


def ensure_sequence_simulated(ts_path: Path, output_dir: Path, mutation_rate: float = 1.9086e-08) -> tuple[Path, bool]:
    """
    Ensures the tree sequence at ts_path (expected at output_dir/<exp_rel>/ts/<name>.ts) is
    simulated to output_dir/<exp_rel>/simulated/<name>.fa, simulating it if missing. No
    candidate-path search and no move/rename side effect - errors if ts_path itself is missing.
    Returns (target_fa, freshly_simulated) - freshly_simulated is False when the FASTA already
    existed, so callers can skip re-verifying work (e.g. Neoaves filtering) already done for it.
    """
    exp_rel = get_experiment_relative_path(ts_path)
    name = get_base_name(ts_path)
    sim_dir = output_dir / exp_rel / "simulated"
    target_fa = sim_dir / f"{name}.fa"

    if target_fa.exists():
        print(f"[Simulate] Simulated FASTA already exists: {target_fa}", file=sys.stderr)
        return target_fa, False

    if not ts_path.exists():
        sys.exit(f"Error: Tree sequence file does not exist: {ts_path}")

    sim_dir.mkdir(parents=True, exist_ok=True)
    print(f"[Simulate] Simulating {ts_path} -> {target_fa}...", file=sys.stderr)
    simulate_ts_to_fasta(ts_path, target_fa, mutation_rate=mutation_rate)
    return target_fa, True


def ensure_neoaves_filtered(fasta_path: Path, freshly_simulated: bool = True) -> Path:
    """
    Filters a freshly-simulated FASTA for Neoaves (removes outgroup taxa) in-place if needed.
    Note: Alt sequences do not require Neoaves filtering (already simulated on Neoaves tree).
    A FASTA that already existed on disk (freshly_simulated=False) was already filtered the first
    time it was simulated, since simulation and filtering are always paired - re-verifying its
    content here would mean re-reading the whole (often 100+MB) file for nothing on every run.
    """
    if not freshly_simulated:
        return fasta_path

    base_name = get_base_name(fasta_path)

    # Alt sequences do not need Neoaves filtering (silent check)
    if base_name.startswith("alt") or "anomaly" in base_name.lower():
        return fasta_path

    fasta_data = read_fasta(fasta_path)
    has_outgroup = any(h in OUTGROUP_TAXA for h in fasta_data)

    if not has_outgroup:
        print(f"[3/4 Check Neoaves Filter] FASTA is already Neoaves-filtered: {fasta_path}", file=sys.stderr)
        return fasta_path

    print(f"[3/4 Check Neoaves Filter] Filtering Neoaves species from {fasta_path}...", file=sys.stderr)
    filtered_data = {h: seq for h, seq in fasta_data.items() if h not in OUTGROUP_TAXA}
    write_fasta(filtered_data, fasta_path)
    print(f"[3/4 Check Neoaves Filter] Filter complete: retained {len(filtered_data)} / {len(fasta_data)} species.", file=sys.stderr)
    return fasta_path


def _rename_ts_dir(ts_dir: Path) -> list[Path]:
    """
    Renames any non-numerically-named .ts file in ts_dir to its canonical chunk-index name
    (matching what organize_loose_ts_files() would have produced), leaving already-numeric names
    untouched. Errors if a chunk index can't be determined, or the canonical name is already
    taken by a different file. Returns the final list of .ts files in ts_dir.
    """
    if not ts_dir.is_dir():
        return []

    ts_files = []
    for ts_file in sorted(ts_dir.iterdir()):
        if not ts_file.is_file():
            continue
        stem = get_base_name(ts_file)
        if stem.isdigit():
            ts_files.append(ts_file)
            continue

        chunk_idx = get_chunk_index(stem) or get_chunk_index(Path(generate_new_filename(ts_file)).stem)
        if not chunk_idx:
            sys.exit(f"Error: Could not determine chunk index for '{ts_file}'.")

        canonical = ts_dir / f"{chunk_idx}.ts"
        if canonical.exists() and canonical != ts_file:
            sys.exit(f"Error: Cannot rename '{ts_file}' -> '{canonical}': target already exists.")

        print(f"[Simulate] Renaming {ts_file.name} -> {canonical.name}", file=sys.stderr)
        ts_file.rename(canonical)
        ts_files.append(canonical)

    return ts_files


_GTREES_WARNINGS: list[str] = []


def ensure_gtrees_concat(null_path: Path, alt_path: Path | None, output_dir: Path, **concat_kwargs) -> Path | None:
    """
    Builds the pre-made gene tree concat (concat/<name>.gtrees) matching a sequence concat, when
    null/gtrees/ and (for an alt experiment) <exp>/gtrees/ exist. Missing chunks or too few alt
    trees are recorded as run-summary warnings instead of stopping the pipeline.
    """
    if not (output_dir / "null" / "gtrees").is_dir():
        return None
    if alt_path is not None and not (Path(alt_path) / "gtrees").is_dir():
        return None
    try:
        return concatenate_gtrees_pair(null_path, alt_path, output_dir, **concat_kwargs)
    except SystemExit as e:
        # The same missing-chunk error repeats for every interval of an experiment - report it once.
        exp_label = Path(alt_path).relative_to(output_dir) if alt_path is not None else "null"
        warning = f"Gene tree concat skipped for {exp_label}: {e}"
        if warning not in _GTREES_WARNINGS:
            print(f"[Simulate] Warning: {warning}", file=sys.stderr)
            _GTREES_WARNINGS.append(warning)
        return None


def _concat(null_path: Path, alt_path: Path | None, output_dir: Path, **concat_kwargs) -> Path:
    """concatenate_pair() plus its matching pre-made gene tree concat. Returns the FASTA path."""
    fa_path = concatenate_pair(null_path, alt_path, output_dir, **concat_kwargs)
    ensure_gtrees_concat(null_path, alt_path, output_dir, **concat_kwargs)
    return fa_path


def _gtrees_only_sweep(output_dir: Path, exp_rels: set[Path], proportions: list[float]) -> None:
    """Default six-proportion .gtrees concat sweep for experiments that only received new .gtrees
    files this run (no new .ts files, so the sequence pipeline didn't visit them)."""
    for exp_rel in sorted(exp_rels):
        alt_path = None if exp_rel == Path("null") else output_dir / exp_rel
        print(f"\n[Simulate -r] Running default gene tree concat sweep for {exp_rel}...", file=sys.stderr)
        for p in proportions:
            ensure_gtrees_concat(Path("null"), alt_path, output_dir, interval=proportion_to_centered_interval(p))


def _print_run_summary() -> None:
    """Prints every alt-chunk-unavailable warning recorded during this run, or a success message
    if there were none (i.e. every requested chunk was found and simulated)."""
    warnings = get_missing_alt_chunk_warnings() + _GTREES_WARNINGS
    if warnings:
        print(f"\n[Simulate] {len(warnings)} alt chunk warning(s) during this run:")
        for w in warnings:
            print(f"  {w}")
    else:
        print(f"\n[Simulate] Success: all requested chunks were fully simulated.")


def _finish_experiment(output_dir: Path, exp_dir: Path, ts_files: list[Path], mutation_rate: float, proportions: list[float]) -> list[Path]:
    """
    Shared by run_all() and run_ensure_all(): derives exp_dir's branch/node label strictly from
    its directory name, ensures its mapping, simulates+filters every given ts file, then runs the
    default six-proportion concat sweep. Returns the concat files created.
    """
    rel_dir = exp_dir.relative_to(output_dir)
    node_label = extract_branch_name_from_filename(rel_dir) or extract_branch_name(rel_dir)
    if not node_label:
        sys.exit(f"Error: Could not derive a branch/internal-node label from directory name '{rel_dir}'.")
    ensure_mapping(node_label, rel_dir, output_dir)

    for ts_path in ts_files:
        sim_fa, freshly_simulated = ensure_sequence_simulated(ts_path, output_dir, mutation_rate=mutation_rate)
        ensure_neoaves_filtered(sim_fa, freshly_simulated)

    return [
        _concat(Path("null"), exp_dir, output_dir, interval=proportion_to_centered_interval(p))
        for p in proportions
    ]


def run_all(output_dir: Path, mutation_rate: float, proportions: list[float] = DEFAULT_CENTERED_PROPORTIONS) -> list[Path]:
    """
    -r/--run-all: organizes every loose .ts file into its experiment directory, simulates only
    those files to FASTA (always Neoaves-filtered), ensures each experiment's mapping, then runs
    the default six-proportion concat sweep for each experiment. Single fixed behavior - the
    caller (main()) rejects any pattern/interval override before this is called.
    """
    from organize_simulations import organize_loose_ts_files, organize_loose_gtrees_files
    all_ts_files = organize_loose_ts_files(output_dir)
    gtrees_exp_rels = {p.parent.parent.relative_to(output_dir) for p in organize_loose_gtrees_files(output_dir)}

    if not all_ts_files:
        print(f"No loose tree sequence (.ts) files found under '{output_dir}'.", file=sys.stderr)
        _gtrees_only_sweep(output_dir, gtrees_exp_rels, proportions)
        _print_run_summary()
        return []

    print(f"[Simulate -r] Found {len(all_ts_files)} loose tree sequence file(s) under '{output_dir}'. Processing...", file=sys.stderr)

    exp_to_ts = {}
    for ts_path in all_ts_files:
        exp_rel = get_experiment_relative_path(ts_path)
        exp_to_ts.setdefault(exp_rel, []).append(ts_path)

    concat_files_created = []
    for exp_rel, ts_files in exp_to_ts.items():
        print(f"\n==========================================", file=sys.stderr)
        print(f"[Simulate -r] Processing Experiment: {exp_rel}", file=sys.stderr)
        print(f"==========================================", file=sys.stderr)

        exp_dir = output_dir / exp_rel

        if exp_rel == Path("null"):
            for ts_path in ts_files:
                sim_fa, freshly_simulated = ensure_sequence_simulated(ts_path, output_dir, mutation_rate=mutation_rate)
                ensure_neoaves_filtered(sim_fa, freshly_simulated)
            print(f"\n[Simulate -r] Running default concat sweep for null...", file=sys.stderr)
            for proportion in proportions:
                interval = proportion_to_centered_interval(proportion)
                concat_files_created.append(_concat(Path("null"), None, output_dir, interval=interval))
        else:
            print(f"\n[Simulate -r] Running default concat sweep for {exp_rel}...", file=sys.stderr)
            concat_files_created.extend(_finish_experiment(output_dir, exp_dir, ts_files, mutation_rate, proportions))

    _gtrees_only_sweep(output_dir, gtrees_exp_rels - set(exp_to_ts), proportions)

    print(f"\n[Simulate -r] All simulations completed successfully.", file=sys.stderr)
    print(f"\n[Simulate -r] Concat files created ({len(concat_files_created)}):")
    for cf in concat_files_created:
        print(f"  {cf}")
    _print_run_summary()
    return concat_files_created


def run_ensure_all(output_dir: Path, mutation_rate: float, proportions: list[float] = DEFAULT_CENTERED_PROPORTIONS) -> list[Path]:
    """
    No-target mode: for every already-organized node experiment directory
    (10X/{up,down}/*, recombination/{up,down}/*, admixture/{low,high}/*) - null is skipped, and
    loose .ts files outside a node directory are never touched, that's -r's job - renames any
    non-numerically-named .ts file in its ts/ folder, simulates/filters every chunk, ensures its
    mapping, and runs the default six-proportion concat sweep.
    """
    exp_dirs = []
    for category in ("10X", "recombination"):
        for direction in ("up", "down"):
            base = output_dir / category / direction
            if base.is_dir():
                exp_dirs.extend(sorted(p for p in base.iterdir() if p.is_dir()))
    for level in ("low", "high"):
        base = output_dir / "admixture" / level
        if base.is_dir():
            exp_dirs.extend(sorted(p for p in base.iterdir() if p.is_dir()))

    concat_files_created = []
    for exp_dir in exp_dirs:
        rel_dir = exp_dir.relative_to(output_dir)
        print(f"\n[Simulate] Ensuring experiment complete: {rel_dir}", file=sys.stderr)
        ts_files = _rename_ts_dir(exp_dir / "ts")
        concat_files_created.extend(_finish_experiment(output_dir, exp_dir, ts_files, mutation_rate, proportions))

    print(f"\n[Simulate] {len(exp_dirs)} experiment director{'y' if len(exp_dirs) == 1 else 'ies'} ensured complete.")
    print(f"[Simulate] Concat files created ({len(concat_files_created)}):")
    for cf in concat_files_created:
        print(f"  {cf}")
    _print_run_summary()
    return concat_files_created


def run_null(output_dir: Path, np_pattern: str, ap_pattern: str, proportion, interval, mutation_rate: float) -> list[Path]:
    """The literal 'null' target: simulates/filters the needed null chunks, then concatenates -
    an explicit interval if given, else the default six-proportion sweep."""
    matches = re.findall(r'[an](\d+)', np_pattern.lower())
    null_chunk_indices = [int(m) for m in matches] if matches else list(range(1, 13))
    for idx in null_chunk_indices:
        ts_path = output_dir / "null" / "ts" / f"{idx}.ts"
        fa, freshly_simulated = ensure_sequence_simulated(ts_path, output_dir, mutation_rate=mutation_rate)
        ensure_neoaves_filtered(fa, freshly_simulated)

    if interval:
        return [_concat(Path("null"), None, output_dir, np_pattern=np_pattern, ap_pattern=ap_pattern, proportion=proportion, interval=interval)]

    return [
        _concat(Path("null"), None, output_dir, np_pattern=np_pattern, ap_pattern=ap_pattern, interval=proportion_to_centered_interval(p))
        for p in DEFAULT_CENTERED_PROPORTIONS
    ]


def run_relative_dir(output_dir: Path, rel_dir: Path, np_pattern: str, ap_pattern: str, proportion, interval, mutation_rate: float) -> list[Path]:
    """
    The <relative_dir> target: must be a relative path to an already-existing experiment
    directory under output_dir - this mode never creates a new experiment from scratch, that's
    what -r is for. Simulates/maps the alt chunks in that directory, ensures the shared null
    chunks needed, then concatenates - an explicit interval if given, else the default
    six-proportion sweep. Branch/mapping is derived strictly from the directory name.
    """
    if rel_dir.is_absolute():
        sys.exit(f"Error: Target must be a relative path under the simulations directory, got absolute path: {rel_dir}")

    exp_dir = output_dir / rel_dir
    if not exp_dir.is_dir():
        sys.exit(f"Error: Experiment directory does not exist: {exp_dir}. Run -r first to organize new tree sequences.")

    computed_exp_rel = get_experiment_relative_path(exp_dir)
    if computed_exp_rel != rel_dir:
        sys.exit(
            f"Error: '{rel_dir}' does not match its deterministic experiment path "
            f"'{computed_exp_rel}'. The directory may be misplaced or misnamed."
        )

    node_label = extract_branch_name_from_filename(rel_dir) or extract_branch_name(rel_dir)
    if not node_label:
        sys.exit(f"Error: Could not derive a branch/internal-node label from directory name '{rel_dir}'.")

    ensure_mapping(node_label, rel_dir, output_dir)

    alt_matches = re.findall(r'[an](\d+)', ap_pattern.lower())
    alt_chunk_indices = [int(m) for m in alt_matches] if alt_matches else [1, 2, 3]
    for idx in alt_chunk_indices:
        ts_path = exp_dir / "ts" / f"{idx}.ts"
        fa, freshly_simulated = ensure_sequence_simulated(ts_path, output_dir, mutation_rate=mutation_rate)
        ensure_neoaves_filtered(fa, freshly_simulated)

    null_matches = re.findall(r'[an](\d+)', np_pattern.lower())
    null_chunk_indices = [int(m) for m in null_matches] if null_matches else list(range(1, 13))
    for idx in null_chunk_indices:
        ts_path = output_dir / "null" / "ts" / f"{idx}.ts"
        fa, freshly_simulated = ensure_sequence_simulated(ts_path, output_dir, mutation_rate=mutation_rate)
        ensure_neoaves_filtered(fa, freshly_simulated)

    if interval:
        return [_concat(Path("null"), exp_dir, output_dir, np_pattern=np_pattern, ap_pattern=ap_pattern, proportion=proportion, interval=interval)]

    return [
        _concat(Path("null"), exp_dir, output_dir, np_pattern=np_pattern, ap_pattern=ap_pattern, interval=proportion_to_centered_interval(p))
        for p in DEFAULT_CENTERED_PROPORTIONS
    ]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Run the sequence simulation pipeline. Four modes, selected by exactly one "
            "positional target: (no target) ensures every already-organized node experiment "
            "directory (10X/recombination/admixture; null is skipped) is complete - renames "
            "non-numeric .ts files, simulates/filters every chunk, ensures mapping, and runs "
            "the default six-proportion concat sweep; -r/--run-all organizes every loose tree "
            "sequence file into its experiment directory first, then does the same for every "
            "experiment including null; 'null' simulates/filters/concats the null experiment "
            "only; <relative_dir> simulates/maps/concats an already-existing experiment "
            "directory (relative to the simulations dir) against null, honoring --np/--ap/-p/"
            "interval overrides."
        )
    )
    parser.add_argument(
        "target",
        type=str,
        nargs="?",
        default=None,
        help="'null', or a relative path (under the simulations dir) to an existing experiment directory. Omit entirely to ensure every node experiment directory is complete (see -r to organize new loose tree sequences first).",
    )
    parser.add_argument(
        "interval",
        type=str,
        nargs="?",
        default=None,
        help="Optional explicit insertion interval (e.g. '40-60'). If omitted, the default six-proportion centered sweep is used. Only valid with target 'null' or <relative_dir>.",
    )
    parser.add_argument(
        "-p", "--prop", "--proportion",
        dest="proportion",
        type=str,
        default=None,
        help="Proportion of alternative sequence (e.g. 0.1, 0.25, 3/12). On its own it replaces the default six-proportion sweep with that proportion's centered interval (0.1 -> 45-55). Works with every mode."
    )
    parser.add_argument(
        "--np", "--null-pattern", type=str, default="n1n2n3n4n5n6n7n8n9n10n11n12", help="Null pattern string specifying null chunk order."
    )
    parser.add_argument(
        "--ap", "--alt-pattern", type=str, default="a1a2a3", help="Alternative pattern string specifying alt chunk order."
    )
    parser.add_argument(
        "-mrate", "--mutation-rate",
        type=float,
        default=1.9086e-08,
        help="Mutation rate for sequence simulation."
    )
    parser.add_argument(
        "-r", "--recursive", "--run-all",
        action="store_true",
        dest="run_all",
        help="Organize every loose tree sequence file and simulate/map/concat each experiment (including null) with the default centered-proportion sweep. Accepts no target, interval, or --np/--ap/-p overrides."
    )
    parser.add_argument(
        "-g", "--gene-trees",
        action="store_true",
        help="Last step: estimate gene trees (IQ-TREE 2, GTR+G, --fast) from 500bp windows every 8.5kb of the simulated null and alt chunks, then insert the alt trees for each concat interval of this run: <exp_dir>/genetrees/<interval>.gtrees. Works with every mode."
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help="Seed for the random window start offsets used by --gene-trees."
    )
    parser.add_argument(
        "--gt-threads",
        type=int,
        default=None,
        help="Max parallel IQ-TREE jobs for --gene-trees (default: all CPUs)."
    )
    args = parser.parse_args(argv)

    output_dir = get_output_dir()

    # -p on its own means the centered interval for that proportion (e.g. 0.1 -> 45-55) instead
    # of the default six-proportion sweep.
    proportions = [parse_proportion(args.proportion)] if args.proportion else DEFAULT_CENTERED_PROPORTIONS
    if args.proportion and not args.interval:
        args.interval = proportion_to_centered_interval(proportions[0])

    has_overrides = bool(
        (args.interval and not args.proportion)
        or args.np != "n1n2n3n4n5n6n7n8n9n10n11n12"
        or args.ap != "a1a2a3"
    )

    if args.run_all:
        if args.target is not None or has_overrides:
            sys.exit("Error: -r/--run-all does not accept a target, interval, or --np/--ap overrides - it organizes and simulates every experiment with the default sweep, or only -p's centered interval if given.")
        concat_files = run_all(output_dir, args.mutation_rate, proportions)
    elif args.target is None:
        if has_overrides:
            sys.exit("Error: No target given. Pass 'null' or a relative experiment directory to use --np/--ap/interval; omit them to ensure every node experiment directory is complete (-p is allowed).")
        concat_files = run_ensure_all(output_dir, args.mutation_rate, proportions)
    else:
        if args.target == "null":
            concat_files = run_null(output_dir, args.np, args.ap, args.proportion, args.interval, args.mutation_rate)
        else:
            concat_files = run_relative_dir(output_dir, Path(args.target), args.np, args.ap, args.proportion, args.interval, args.mutation_rate)

        print(f"\n[Simulate] Concat files created ({len(concat_files)}):")
        for cf in concat_files:
            print(f"  {cf}")
        _print_run_summary()

    if args.gene_trees:
        # Gene trees always use the default null/alt chunk order; only each concat file's
        # interval (the part of its name after any ';') carries over.
        exp_to_intervals = {}
        for cf in concat_files:
            exp_to_intervals.setdefault(cf.parent.parent, []).append(cf.stem.split(";")[-1])
        print(f"\n[Simulate] Estimating gene trees for {len(exp_to_intervals)} experiment(s)...", file=sys.stderr)
        for exp_dir, intervals in exp_to_intervals.items():
            is_null = exp_dir == output_dir / "null"
            ensure_estimated_gene_trees(output_dir, None if is_null else exp_dir, None if is_null else intervals, args.gt_threads, args.seed)


if __name__ == "__main__":
    main()

