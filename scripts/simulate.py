#!/usr/bin/env python3
"""
Simulate command pipeline:
1. Checks if mapping is generated (generates via internal generate_mapping if missing).
2. Checks if tree sequences (.ts) are simulated to FASTA (simulates via simulate_sequences if missing).
3. Checks if FASTA sequences are filtered for Neoaves (filters outgroup taxa via internal filter_neoaves if missing).
4. Runs sequence concatenation via concat_sequences.
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
    construct_indexed_path,
    find_existing_chunk_indices,
    get_phlag_dir,
    get_experiment_folder_name,
    get_experiment_relative_path,
    get_chunk_index,
    main as concat_main,
)
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


def normalize_and_rename_file(file_path: str | Path | None) -> Path | None:
    """
    Checks if a file path or legacy filename should be renamed to standard format.
    Renames the file on disk if it exists and needs renaming, preserving existing files.
    """
    if not file_path:
        return None
    p = Path(file_path)
    if not p.exists():
        return p
    new_name = generate_new_filename(p)
    if new_name == p.name:
        return p

    new_path = p.parent / new_name
    if p != new_path:
        if new_path.exists():
            print(f"[Auto-Rename] Target standard file already exists: {new_path}, keeping {p.name}", file=sys.stderr)
            return new_path
        else:
            print(f"[Auto-Rename] Standardizing filename on disk: {p.name} -> {new_name}", file=sys.stderr)
            try:
                p.rename(new_path)
            except Exception as e:
                print(f"[Auto-Rename Warning] Could not rename {p} -> {new_path}: {e}", file=sys.stderr)
                return p
            return new_path
    return p


def ensure_mapping(tree_path: Path, node_label: str, mapping_path: Path = None, exp_name: str = None, output_dir: Path = None) -> Path:
    """Checks if mapping file exists, otherwise generates it using species tree & internal node."""
    tree_stem = "neoaves"
    formatted_label = node_label
    if "_" in node_label:
        parts = node_label.split("_")
        target_node = parts[0]
        source_node = parts[1]
        clade = get_admixture_source_clade(tree_path, target_node, source_node)
        formatted_label = f"{target_node}:{clade}"

    target_name = f"{tree_stem}_{formatted_label}_mapping.tsv"
    
    if output_dir is None:
        output_dir = get_output_dir()
    phlag_dir = get_phlag_dir()
    exp_rel = get_experiment_relative_path(exp_name or node_label) if (exp_name or node_label) else None
    desired_dir = output_dir / exp_rel if exp_rel else (output_dir / exp_name if exp_name else output_dir)

    if mapping_path is None:
        candidates = []
        if exp_rel:
            candidates.append(output_dir / exp_rel / target_name)
            candidates.append(phlag_dir / "simulations" / exp_rel / target_name)
        if exp_name:
            candidates.append(output_dir / exp_name / target_name)
            candidates.append(phlag_dir / "simulations" / exp_name / target_name)
        candidates.append(output_dir / target_name)

        for cand in candidates:
            if cand.exists():
                mapping_path = cand
                break
                
        if mapping_path is None:
            desired_dir.mkdir(parents=True, exist_ok=True)
            mapping_path = desired_dir / target_name

    # Relocate mapping file to standard directory structure if it currently exists in a non-relativized path
    if mapping_path.exists():
        if desired_dir and mapping_path.parent != desired_dir:
            desired_dir.mkdir(parents=True, exist_ok=True)
            desired_path = desired_dir / target_name
            if not desired_path.exists():
                print(f"[Simulate] Moving mapping file to directory format: {mapping_path} -> {desired_path}", file=sys.stderr)
                mapping_path.rename(desired_path)
                mapping_path = desired_path
            else:
                mapping_path = desired_path
        print(f"[1/4 Check Mapping] Mapping file already exists: {mapping_path}", file=sys.stderr)
    else:
        if not tree_path.exists():
            sys.exit(f"Error: Species tree file does not exist: {tree_path}")
        print(f"[1/4 Check Mapping] Mapping file missing. Generating mapping for tree='{tree_path}' node='{node_label}' -> {mapping_path}...", file=sys.stderr)
        generate_mapping(tree_path, node_label, mapping_path)
    return mapping_path


def ensure_sequence_simulated(input_file: Path, output_dir: Path, mutation_rate: float = 1.9086e-08, skip_neoaves: bool = False, design: str = None) -> Path:
    """
    Checks if a sequence file is simulated as FASTA in output_dir/<exp_rel>/simulated/<base_name>.fa.
    If only .ts exists (or path to .ts), simulates it to FASTA.
    """
    base_name = get_base_name(input_file)
    exp_name = get_experiment_folder_name(input_file)
    exp_rel = get_experiment_relative_path(input_file)
    chunk_idx = get_chunk_index(base_name)
    simulated_name = f"{chunk_idx}.fa" if chunk_idx else f"{base_name}.fa"

    # Throw error if alternative sequence chunk index > 3
    s_lower = str(input_file).lower()
    is_alt_by_name = input_file.name.lower().startswith("a") or any(kw in s_lower for kw in ("anomaly", "alt", "anom", "admixture", "recombination", "10x"))
    is_alt = is_alt_by_name or (exp_name != "null" and input_file.name.startswith("alt"))
    if is_alt and chunk_idx and chunk_idx.isdigit() and int(chunk_idx) > 3:
        sys.exit(f"Error: Cannot simulate alternative chunk '{input_file.name}': chunk index {chunk_idx} exceeds maximum allowed index of 3.")

    sim_dir = output_dir / exp_rel / "simulated"
    sim_dir.mkdir(parents=True, exist_ok=True)

    # If already a FASTA file that exists
    if input_file.exists() and input_file.is_file() and (input_file.suffix == ".fa" or str(input_file).endswith(".fa")):
        print(f"[2/4 Check Simulation] FASTA sequence file exists: {input_file}", file=sys.stderr)
        return input_file

    target_fa = sim_dir / simulated_name
    if target_fa.exists():
        print(f"[2/4 Check Simulation] Simulated FASTA already exists: {target_fa}", file=sys.stderr)
        return target_fa

    # Look for tree sequence (.ts) file
    ts_path = input_file
    if not (ts_path.exists() and ts_path.is_file() and (ts_path.suffix == ".ts" or str(ts_path).endswith(".ts"))):
        phlag_dir = get_phlag_dir()
        ts_candidates = []
        if chunk_idx:
            ts_candidates.extend([
                output_dir / exp_rel / "ts" / f"{chunk_idx}.ts",
                output_dir / exp_rel / f"{chunk_idx}.ts",
                output_dir / exp_rel / "ts" / f"default-500Kb-{chunk_idx}.ts",
                output_dir / exp_rel / "ts" / f"default-{chunk_idx}.ts",
                output_dir / exp_rel / "ts" / f"null-{chunk_idx}.ts",
                output_dir / exp_rel / "ts" / f"favian-500Kb-{chunk_idx}.ts",
                output_dir / exp_name / "ts" / f"{chunk_idx}.ts",
                output_dir / exp_name / f"{chunk_idx}.ts",
                output_dir / exp_name / "ts" / f"default-500Kb-{chunk_idx}.ts",
                output_dir / exp_name / "ts" / f"default-{chunk_idx}.ts",
                output_dir / exp_name / "ts" / f"null-{chunk_idx}.ts",
                output_dir / exp_name / "ts" / f"favian-500Kb-{chunk_idx}.ts",
                output_dir / "null" / "ts" / f"default-500Kb-{chunk_idx}.ts",
                output_dir / "null" / "ts" / f"default-{chunk_idx}.ts",
                output_dir / "null" / "ts" / f"null-{chunk_idx}.ts",
                output_dir / "null" / "ts" / f"{chunk_idx}.ts",
            ])
        ts_candidates.extend([
            output_dir / exp_rel / "ts" / f"{base_name}.ts",
            output_dir / exp_rel / f"{base_name}.ts",
            output_dir / exp_name / "ts" / f"{base_name}.ts",
            output_dir / exp_name / f"{base_name}.ts",
            output_dir / "ts" / f"{base_name}.ts",
            input_file.with_suffix(".ts"),
            input_file.parent / f"{base_name}.ts",
        ])

        # Also check parents up to output_dir for loose .ts files
        p_curr = input_file.parent if input_file.is_file() else input_file
        while p_curr and p_curr != output_dir and p_curr != p_curr.parent:
            if chunk_idx:
                ts_candidates.append(p_curr / f"{chunk_idx}.ts")
            ts_candidates.append(p_curr / f"{base_name}.ts")
            p_curr = p_curr.parent
        
        if design:
            ts_candidates.insert(0, phlag_dir / design / "ts" / f"{base_name}.ts")
            if chunk_idx:
                ts_candidates.insert(0, phlag_dir / design / "ts" / f"default-500Kb-{chunk_idx}.ts")
                ts_candidates.insert(0, phlag_dir / design / "ts" / f"{chunk_idx}.ts")
        
        ts_candidates.append(phlag_dir / "ts" / f"{base_name}.ts")
        
        for cand in ts_candidates:
            if cand.exists():
                ts_path = cand
                break

    if ts_path.exists() and (ts_path.suffix == ".ts" or str(ts_path).endswith(".ts")):
        desired_ts_dir = output_dir / exp_rel / "ts" if exp_rel else (output_dir / exp_name / "ts" if exp_name else ts_path.parent)
        desired_name = ts_path.name
        if exp_name == "null" or ts_path.name.startswith("default"):
            desired_name = f"null-{chunk_idx}.ts" if chunk_idx else "null.ts"
        elif chunk_idx and not ts_path.name.startswith(base_name):
            desired_name = f"{chunk_idx}.ts"
            
        desired_path = desired_ts_dir / desired_name
        if ts_path != desired_path:
            desired_ts_dir.mkdir(parents=True, exist_ok=True)
            if desired_path.exists():
                print(f"[Simulate] Overwriting existing tree sequence: {ts_path} -> {desired_path}", file=sys.stderr)
                desired_path.unlink()
            else:
                print(f"[Simulate] Moving tree sequence {ts_path} -> {desired_path}", file=sys.stderr)
            ts_path.rename(desired_path)
            ts_path = desired_path

        print(f"[2/4 Check Simulation] Simulating sequence for {ts_path} -> {target_fa}...", file=sys.stderr)
        simulate_ts_to_fasta(ts_path, target_fa, mutation_rate=mutation_rate)
        return target_fa
    elif (input_file.parent / simulated_name).exists():
        return input_file.parent / simulated_name
    elif (input_file.parent / f"{base_name}.fa").exists():
        return input_file.parent / f"{base_name}.fa"
    else:
        sys.exit(f"Error: Neither FASTA nor tree sequence file exists for base name: '{base_name}' ({input_file})")


def ensure_neoaves_filtered(fasta_path: Path) -> Path:
    """
    Checks if a FASTA sequence file is filtered for Neoaves (outgroup taxa removed).
    Filters it in-place (writing to the same FASTA file) if not already filtered.
    Note: Alt sequences do not require Neoaves filtering (already simulated on Neoaves tree).
    """
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


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run sequence simulation pipeline: check mapping, simulate tree sequences, filter Neoaves, and concat sequences."
    )

    parser.add_argument(
        "positional_args",
        type=str,
        nargs="*",
        help="Sequence base names or files to concatenate, or (null_base, alt_base, pattern).",
    )

    parser.add_argument("-n", "--null", type=Path, help="Base name or path for null sequence (without .ts or .fa).")
    parser.add_argument(
        "-a",
        "--alt",
        "--alternative",
        "--anomalous",
        type=Path,
        dest="alt",
        help="Base name or path for alternative/anomalous sequence (without .ts or .fa).",
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
    parser.add_argument(
        "-t", "--tree", "--input-tree",
        dest="input_tree",
        type=Path,
        default=None,
        help="Species tree path for mapping check."
    )
    parser.add_argument(
        "-b", "--branch", "--node", "--internal-node",
        dest="internal_node",
        type=str,
        default=None,
        help="Internal node label for mapping check (extracted from filename if omitted)."
    )
    parser.add_argument(
        "--extract-branch", "--branch-only",
        action="store_true",
        dest="extract_branch_only",
        help="Output extracted branch name for given sequence file(s) and exit."
    )
    parser.add_argument(
        "-m", "--mapping",
        type=Path,
        default=None,
        help="Output or existing mapping TSV path."
    )
    parser.add_argument(
        "-mrate", "--mutation-rate",
        type=float,
        default=1.9086e-08,
        help="Mutation rate for sequence simulation."
    )
    parser.add_argument(
        "--skip-neoaves",
        action="store_true",
        help="Skip Neoaves outgroup filtering."
    )
    parser.add_argument(
        "--design",
        type=str,
        default=None,
        help="High-level design name (model parameters) under PHLAG_DIR to look for tree sequences."
    )
    parser.add_argument(
        "-r", "--recursive", "--run-all",
        action="store_true",
        dest="run_all",
        help="Simulate all tree sequence (.ts) files recursively in store/simulations."
    )
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output FASTA path.")
    args = parser.parse_args(argv)

    output_dir = get_output_dir()

    if args.run_all:
        all_ts_files = []
        try:
            from organize_simulations import organize_loose_ts_files
            all_ts_files = organize_loose_ts_files(output_dir)
        except ImportError as e:
            print(f"[Warning] Could not import organize_simulations: {e}", file=sys.stderr)

        if not all_ts_files:
            print(f"No loose tree sequence (.ts) files found under '{output_dir}'.", file=sys.stderr)
            sys.exit(0)

        print(f"[Simulate -r] Found {len(all_ts_files)} loose tree sequence file(s) under '{output_dir}'. Processing...", file=sys.stderr)
        
        exp_to_ts = {}
        for ts_path in all_ts_files:
            exp_name = get_experiment_folder_name(ts_path)
            exp_to_ts.setdefault(exp_name, []).append(ts_path)

        for exp_name, ts_files in exp_to_ts.items():
            print(f"\n==========================================", file=sys.stderr)
            print(f"[Simulate -r] Processing Experiment: {exp_name}", file=sys.stderr)
            print(f"==========================================", file=sys.stderr)
            
            node_label = args.internal_node
            if not node_label and exp_name != "null":
                for ts_path in ts_files:
                    extracted = extract_branch_name_from_filename(ts_path) or extract_branch_name(ts_path)
                    if extracted:
                        node_label = extracted
                        break

            tree_path = args.input_tree or (REPO_ROOT / "speciestree" / "main-neoaves-numgen.nwk")
            if node_label and exp_name != "null":
                ensure_mapping(tree_path, node_label, args.mapping, exp_name=str(ts_files[0]), output_dir=output_dir)

            for ts_path in ts_files:
                sim_fa = ensure_sequence_simulated(ts_path, output_dir, mutation_rate=args.mutation_rate, skip_neoaves=args.skip_neoaves, design=args.design)
                if not args.skip_neoaves:
                    ensure_neoaves_filtered(sim_fa)

            print(f"\n[Simulate -r] Running default concatenation for {exp_name}...", file=sys.stderr)
            concat_argv = ["-a", str(ts_files[0])] if exp_name != "null" else []
            concat_argv.extend(["40-65"])
            if args.design:
                concat_argv.extend(["--design", args.design])
            concat_main(concat_argv)

        print(f"\n[Simulate -r] All simulations completed successfully.", file=sys.stderr)
        sys.exit(0)

    if args.extract_branch_only:
        targets = []
        if args.null:
            targets.append(args.null)
        if args.alt:
            targets.append(args.alt)
        targets.extend([Path(p) for p in args.positional_args if p])

        if not targets:
            sys.exit("Error: Must specify sequence file(s) or -n/-a to extract branch from.")

        for target in targets:
            bname = extract_branch_name_from_filename(target) or extract_branch_name(target)
            print(f"{target}: {bname}")
        sys.exit(0)

    if args.null:
        args.null = normalize_and_rename_file(args.null)
    if args.alt:
        args.alt = normalize_and_rename_file(args.alt)
    pos = [normalize_and_rename_file(Path(p)) for p in args.positional_args if p]

    # Determine input files for simulation / filtering / concat
    coord_parts = []
    file_pos = []
    for p in pos:
        name_clean = p.name.lower().strip()
        is_interval = bool(re.match(r"^[\d\.,\-_kmb]+$", name_clean))
        if is_interval and not p.exists():
            coord_parts.append(p.name)
        else:
            file_pos.append(p)

    null_path = args.null
    if not null_path and len(file_pos) > 1:
        null_path = file_pos.pop(0)

    alt_path = args.alt
    if not alt_path and file_pos:
        alt_path = file_pos.pop(0)

    if not null_path:
        null_path = Path("null")

    if not coord_parts:
        coord_parts = ["40-65"]
        print(f"No interval provided. Defaulting to interval: {coord_parts[0]}", file=sys.stderr)

    interval_str = ",".join(coord_parts)

    # is_pattern_mode is True if we are doing any kind of concatenation
    # that involves `--np`, `--ap`, intervals, or proportions.
    is_pattern_mode = True if (alt_path or coord_parts or args.np != "n1n2n3n4n5n6n7n8n9n10n11n12" or args.ap != "a1a2a3" or args.proportion) else False

    target_files = []
    if is_pattern_mode:
        if alt_path:
            # Parse needed alt chunks from --ap
            matches = re.findall(r'[an](\d+)', args.ap.lower())
            alt_chunk_indices = [int(m) for m in matches] if matches else [1, 2, 3]
            exp_rel = get_experiment_relative_path(alt_path)
            for idx in alt_chunk_indices:
                cf = output_dir / exp_rel / "ts" / f"{idx}.ts"
                target_files.append(cf)
        else:
            print("Note: No alternative sequence explicitly provided. Simulating null sequence blocks.", file=sys.stderr)
        
        # Parse needed null chunks from --np
        matches = re.findall(r'[an](\d+)', args.np.lower())
        null_chunk_indices = [int(m) for m in matches] if matches else [i for i in range(1, 13)]
        for idx in null_chunk_indices:
            cf = output_dir / "null" / "ts" / f"{idx}.ts"
            target_files.append(cf)
    else:
        target_files = [null_path] + [alt_path] if alt_path else [null_path]

    for tf in target_files:
        if tf and tf.name != "null":
            ensure_sequence_simulated(tf, output_dir, mutation_rate=args.mutation_rate, skip_neoaves=args.skip_neoaves, design=args.design)

    # 1. Check Mapping
    tree_path = args.input_tree
    if tree_path is None:
        tree_path = REPO_ROOT / "speciestree" / "main-neoaves-numgen.nwk"

    # Extract branch name from filename/path if available
    node_label = args.internal_node
    if not node_label:
        candidates_to_check = pos + ([alt_path] if alt_path else []) + target_files
        for cand_file in candidates_to_check:
            if cand_file and str(cand_file).lower() != "null":
                extracted = extract_branch_name_from_filename(cand_file) or extract_branch_name(cand_file)
                if extracted and not extracted.isdigit() and extracted not in ("null", "default", "ts", "simulated", "concat", "store"):
                    node_label = extracted
                    break

    exp_names = set()
    for tf in target_files:
        exp_names.add(get_experiment_folder_name(tf))
    if is_pattern_mode:
        if alt_path:
            exp_names.add(get_experiment_folder_name(alt_path))
        if null_path:
            exp_names.add(get_experiment_folder_name(null_path))

    if "null" in exp_names and len(exp_names) > 1:
        exp_names.remove("null")

    if not exp_names:
        exp_names.add("default")

    is_null_only = exp_names == {"null"}

    if not is_null_only:
        if not node_label or node_label.isdigit() or node_label in ("null", "default", "ts", "simulated", "concat", "store"):
            sys.exit(f"Error: Could not extract a valid internal node label from '{alt_path or (pos[0] if pos else 'input')}'. Please specify --node / --branch explicitly.")

    if not is_null_only:
        exp_targets = set()
        if alt_path:
            exp_targets.add(str(alt_path))
        if pos:
            for p in pos:
                if str(p).lower() != "null":
                    exp_targets.add(str(p))
        if not exp_targets:
            for tf in target_files:
                if str(tf).lower() != "null":
                    exp_targets.add(str(tf))
        
        for target in sorted(exp_targets):
            ensure_mapping(tree_path, node_label, args.mapping, exp_name=target, output_dir=output_dir)

    # 2. Check Sequence Simulation & 3. Check Neoaves Filter
    processed_files = []
    for tf in target_files:
        sim_fa = ensure_sequence_simulated(tf, output_dir, mutation_rate=args.mutation_rate, skip_neoaves=args.skip_neoaves, design=args.design)
        if not args.skip_neoaves:
            filtered_fa = ensure_neoaves_filtered(sim_fa)
            processed_files.append(filtered_fa)
        else:
            processed_files.append(sim_fa)

    # 4. Run Concat Sequences
    print(f"\n[4/4 Concat Sequences] Concatenating sequence files...", file=sys.stderr)
    
    if is_pattern_mode:
        concat_argv = ["-n", str(null_path)]
        if alt_path:
            concat_argv.extend(["-a", str(alt_path)])
        concat_argv.extend(["--np", args.np, "--ap", args.ap])

        if coord_parts:
            concat_argv.extend(coord_parts)
        if args.proportion:
            concat_argv.extend(["-p", str(args.proportion)])
    else:
        concat_argv = [str(f) for f in processed_files]

    if args.output:
        concat_argv.extend(["-o", str(args.output)])
    if args.design:
        concat_argv.extend(["--design", args.design])

    concat_main(concat_argv)


if __name__ == "__main__":
    main()

