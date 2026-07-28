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
import treeswift as ts

# Ensure parent directory of scripts is in sys.path
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

REPO_ROOT = SCRIPTS_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from simulate_sequences import simulate_ts_to_fasta, get_output_dir
from concat_sequences import (
    parse_indexed_pattern,
    construct_indexed_path,
    expand_pattern_to_target_files,
    get_file_type,
    get_phlag_dir,
    get_experiment_folder_name,
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


def generate_mapping(tree_path: Path, node_label: str, output_path: Path):
    tree = ts.read_tree_newick(tree_path)
    label_to_node = tree.label_to_node(selection="all")
    
    if node_label not in label_to_node:
        sys.exit(f"Error: Internal node '{node_label}' not found in species tree '{tree_path}'.")
        
    nd_target = label_to_node[node_label]
    
    # P2 = Leaves under target internal node
    p2_leaves = set(nd.get_label() for nd in tree.extract_subtree(nd_target).traverse_leaves())
    
    # P1 = Leaves under sister clade(s) of target internal node
    parent_node = nd_target.get_parent()
    p1_leaves = set()
    if parent_node:
        for child in parent_node.child_nodes():
            if child != nd_target:
                p1_leaves.update(nd.get_label() for nd in tree.extract_subtree(child).traverse_leaves())
                
    # P3 = Leaves under sister clade(s) of parent node
    p3_leaves = set()
    grandparent_node = parent_node.get_parent() if parent_node else None
    if grandparent_node:
        for child in grandparent_node.child_nodes():
            if child != parent_node:
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
    "recombination_decrease": "recombination_down",
    "recombination_suppression": "recombination_down",
    "10x_up": "10X_up",
    "10x_down": "10X_down",
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
    is_concat = filename.startswith("concat_") or ("_a" in filename and ("n" in filename or "a" in filename) and ("_n" in filename or "_a" in filename or filename.endswith("_na")))

    prefix = "alt" if ("anomaly" in filename.lower() or filename.startswith("alt")) else "null"
    index = ""
    m_idx = re.search(r"-(\d+)$", filename)
    if m_idx:
        index = m_idx.group(1)

    branch = None
    change = None
    rate = None
    time = None

    # Handle admixture: rate and time
    m_adm = re.search(r"admixture[-_]([A-Za-z0-9_]+?)[-_](rate\d+|0\d+)[-_](time\d+|\d+)", filename)
    if not m_adm:
        m_adm = re.search(r"([A-Za-z0-9_]+?)[-_]admixture[-_](rate\d+|0\d+)[-_](time\d+|\d+)", filename)
        if m_adm:
            branch = m_adm.group(1)
            change = "admixture"
            rate = m_adm.group(2)
            time = m_adm.group(3)
    else:
        change = "admixture"
        branch = m_adm.group(1)
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
            change = KNOWN_CHANGES_MAP.get(c_raw, "recombination_up")
            branch = m_rec.group(2).split("-500K")[0].split("-500k")[0]
        else:
            m_rec_rev = re.search(r"([A-Za-z0-9_]+)[-_]recombination_(increase|decrease|suppression|up|down)", filename, re.IGNORECASE)
            if m_rec_rev:
                branch = m_rec_rev.group(1)
                c_raw = f"recombination_{m_rec_rev.group(2).lower()}"
                change = KNOWN_CHANGES_MAP.get(c_raw, "recombination_up")

    # Clean up branch name
    if branch:
        branch = re.sub(r"-500[Kk][Bb].*$", "", branch)
        branch = re.sub(r"-\d+$", "", branch)
        branch = branch.strip("-_")

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
        null/alt-(neoaves)-branch_change_(rate-time)-index.ext
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
        new_name = new_name.replace("anomaly_case-", "alt-")
        new_name = new_name.replace("-500Kb", "")
        new_name = new_name.replace("-500kb", "")
        new_name = new_name.replace("recombination_increase", "recombination_up")
        new_name = new_name.replace("recombination_decrease", "recombination_down")
        new_name = new_name.replace("recombination_suppression", "recombination_down")
        new_name = new_name.replace("rate090_time7099554", "rate090-time7099554")
        return f"{new_name}{ext}"

    prefix = meta["prefix"]
    idx_str = f"-{meta['index']}" if meta["index"] else ""

    if prefix == "null" or not meta["change"]:
        return f"null{idx_str}{ext}"

    neoaves_str = "-neoaves" if meta["is_neoaves"] else ""
    branch = meta["branch"]
    change = meta["change"]

    if change == "admixture" and meta["rate"] and meta["time"]:
        change_spec = f"admixture_{meta['rate']}-{meta['time']}"
    else:
        change_spec = change

    return f"alt{neoaves_str}-{branch}_{change_spec}{idx_str}{ext}"


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
                if new_path.exists() and new_path != old_path:
                    old_path.unlink()
                else:
                    old_path.rename(new_path)
            renamed.append((old_path, new_path))

    return renamed


def extract_branch_name_from_filename(path: str | Path | None) -> str | None:
    """Extracts branch name from filename. For admixture / compound branches (e.g., Strigiformes_N297), extracts the left branch (e.g., Strigiformes)."""
    if not path:
        return None
    extracted = extract_branch_name(path)

    if not extracted:
        clean_name = get_base_name(path)
        m = re.search(r"[-_]([A-Za-z0-9]+)_", clean_name)
        if m:
            extracted = m.group(1)

    if extracted and "_" in extracted:
        left_branch = extracted.split("_")[0]
        if left_branch:
            return left_branch
    return extracted


def normalize_and_rename_file(file_path: str | Path | None) -> Path | None:
    """
    Checks if a file path or legacy filename should be renamed to standard format:
        null/alt-(neoaves)-branch_change_(rate-time)-index.ext
    Renames the file on disk if it exists and needs renaming, and returns the updated Path.
    """
    if not file_path:
        return None
    p = Path(file_path)
    new_name = generate_new_filename(p)
    if new_name == p.name:
        return p

    new_path = p.parent / new_name
    if p.exists() and p != new_path:
        if new_path.exists():
            print(f"[Auto-Rename] Target standard file already exists: {new_path}, removing legacy {p.name}", file=sys.stderr)
            try:
                p.unlink()
            except Exception:
                pass
            return new_path
        else:
            print(f"[Auto-Rename] Standardizing filename on disk: {p.name} -> {new_name}", file=sys.stderr)
            try:
                p.rename(new_path)
            except Exception as e:
                print(f"[Auto-Rename Warning] Could not rename {p} -> {new_path}: {e}", file=sys.stderr)
                return p
            return new_path
    return new_path


def ensure_mapping(tree_path: Path, node_label: str, mapping_path: Path = None) -> Path:
    """Checks if mapping file exists, otherwise generates it using species tree & internal node."""
    target_name = f"{tree_path.stem}_{node_label}_mapping.tsv"
    
    if mapping_path is None:
        phlag_dir = get_phlag_dir()
        phlag_mapping_path = phlag_dir / "mapping" / target_name
        local_mapping_path = REPO_ROOT / "mapping" / target_name
        
        if phlag_mapping_path.exists():
            mapping_path = phlag_mapping_path
        elif local_mapping_path.exists():
            mapping_path = local_mapping_path
        else:
            mapping_dir = phlag_dir / "mapping"
            mapping_dir.mkdir(parents=True, exist_ok=True)
            mapping_path = mapping_dir / target_name

    if mapping_path.exists():
        print(f"[1/4 Check Mapping] Mapping file already exists: {mapping_path}", file=sys.stderr)
    else:
        if not tree_path.exists():
            sys.exit(f"Error: Species tree file does not exist: {tree_path}")
        print(f"[1/4 Check Mapping] Mapping file missing. Generating mapping for tree='{tree_path}' node='{node_label}' -> {mapping_path}...", file=sys.stderr)
        generate_mapping(tree_path, node_label, mapping_path)
    return mapping_path


def ensure_sequence_simulated(input_file: Path, output_dir: Path, mutation_rate: float = 1.9086e-08, skip_neoaves: bool = False, design: str = None) -> Path:
    """
    Checks if a sequence file is simulated as FASTA in output_dir/<exp_name>/simulated/<base_name>.fa.
    If only .ts exists (or path to .ts), simulates it to FASTA.
    """
    base_name = get_base_name(input_file)
    exp_name = get_experiment_folder_name(input_file)
    chunk_idx = get_chunk_index(base_name)
    simulated_name = f"{chunk_idx}.fa" if chunk_idx else f"{base_name}.fa"

    sim_dir = output_dir / exp_name / "simulated"
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
        ts_candidates = [
            output_dir / exp_name / "ts" / f"{base_name}.ts",
            output_dir / "ts" / f"{base_name}.ts",
            REPO_ROOT / "store" / "simulations" / exp_name / "ts" / f"{base_name}.ts",
            REPO_ROOT / "store" / "simulations" / "ts" / f"{base_name}.ts",
            REPO_ROOT / "store" / "simulations" / exp_name / f"{base_name}.ts",
            input_file.with_suffix(".ts"),
            input_file.parent / f"{base_name}.ts",
        ]
        
        if design:
            ts_candidates.insert(0, phlag_dir / design / "ts" / f"{base_name}.ts")
        
        ts_candidates.append(phlag_dir / "ts" / f"{base_name}.ts")
        
        for cand in ts_candidates:
            if cand.exists():
                ts_path = cand
                break

    if ts_path.exists() and (ts_path.suffix == ".ts" or str(ts_path).endswith(".ts")):
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
        "-p", "--pattern", type=str, help="Pattern string (e.g., 'a1n5a2a3n8n1' or 'n1a3a2n2')."
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
    parser.add_argument("-o", "--output", type=Path, required=False, help="Output FASTA path.")
    args = parser.parse_args(argv)

    output_dir = get_output_dir()

    if output_dir.exists():
        ts_dirs = [p / "ts" for p in output_dir.iterdir() if p.is_dir()] + [output_dir / "ts"]
        for ts_dir in ts_dirs:
            if ts_dir.exists():
                rename_all_files(ts_dir)

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

    if args.null and not args.alt and pos:
        args.alt = pos[0]
        pos = pos[1:]
    elif args.alt and not args.null and pos:
        args.null = pos[0]
        pos = pos[1:]

    # Determine input files for simulation / filtering / concat
    is_pattern_mode = False
    if args.null and args.alt:
        is_pattern_mode = True
        null_path = args.null
        alt_path = args.alt
        pattern_raw = args.pattern
    elif len(pos) == 3 and not pos[2].exists():
        try:
            parse_indexed_pattern(pos[2].name)
            is_pattern_mode = True
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
            pattern_raw = pos[2].name
        except SystemExit:
            is_pattern_mode = False
    elif len(pos) == 1 or (len(pos) == 2 and not pos[1].exists()):
        # Single sequence file or single alt + pattern
        single_file = pos[0]
        p_type = get_file_type(single_file)
        if p_type == "a" or ("admixture" in single_file.name or "recombination" in single_file.name or "10X" in single_file.name or "10x" in single_file.name or single_file.name.startswith("N")):
            is_pattern_mode = True
            alt_path = single_file
            if len(pos) == 2 and not pos[1].exists():
                pattern_raw = pos[1].name
            else:
                pattern_raw = None
        elif len(pos) == 2 and not pos[1].exists():
            try:
                parse_indexed_pattern(pos[1].name)
                pattern_raw = pos[1].name
                alt_path = single_file
                is_pattern_mode = True
            except SystemExit:
                is_pattern_mode = False

        if is_pattern_mode:
            null_path = args.null if args.null else Path("null")

            target_files, _ = expand_pattern_to_target_files(null_path, alt_path, pattern_raw, output_dir, design=args.design)
        else:
            target_files = [single_file]

    target_files = [] if not ('target_files' in locals()) else target_files

    if is_pattern_mode and not target_files:
        target_files, _ = expand_pattern_to_target_files(null_path, alt_path, pattern_raw, output_dir, design=args.design)
    elif not is_pattern_mode and not target_files:
        if len(pos) >= 2:
            target_files = pos
        elif len(pos) == 1:
            target_files = [pos[0]]
        else:
            parser.print_help()
            sys.exit("\nError: Must specify sequence file(s), or (alt, pattern), or (null, alt, pattern).")

    # 1. Check Mapping
    tree_path = args.input_tree
    if tree_path is None:
        tree_path = REPO_ROOT / "speciestree" / "main-neoaves-numgen.nwk"

    # Extract branch name from filename if available
    node_label = args.internal_node
    if not node_label:
        candidates_to_check = target_files + pos
        if is_pattern_mode:
            candidates_to_check.extend([null_path, alt_path])

        for cand_file in candidates_to_check:
            extracted = extract_branch_name_from_filename(cand_file)
            if extracted:
                node_label = extracted
                break

    if not node_label:
        sys.exit("Error: Internal node label could not be determined. Please specify --node.")

    ensure_mapping(tree_path, node_label, args.mapping)

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
    print(f"\n[4/4 Concat Sequences] Concatenating {len(processed_files)} sequence files...", file=sys.stderr)
    
    # Build args for concat_main
    if is_pattern_mode:
        concat_argv = ["-n", str(null_path), "-a", str(alt_path)]
        if pattern_raw:
            concat_argv.extend(["-p", pattern_raw])
    else:
        concat_argv = [str(f) for f in processed_files]

    if args.output:
        concat_argv.extend(["-o", str(args.output)])
    if args.design:
        concat_argv.extend(["--design", args.design])

    concat_main(concat_argv)


if __name__ == "__main__":
    main()
