#!/usr/bin/env python3
"""
Organizes loose .ts files in the simulations directory into their respective experiment folders.
Creates ts/, simulated/, and concat/ subdirectories, and renames .ts files to <chunk_idx>.ts inside the ts/ folder.
Loose pre-made gene tree files (.gtrees) are moved the same way to <exp>/gtrees/<chunk_idx>.gtrees.
"""

import sys
import filecmp
import argparse
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from concat_sequences import get_experiment_folder_name, get_experiment_relative_path, get_chunk_index, classify_gtrees_file, is_joined_gtrees_file, GTREES_EXT
from simulate import generate_new_filename

def organize_loose_ts_files(sim_dir: Path) -> list[Path]:
    if not sim_dir.exists():
        return []

    # Find all .ts files under sim_dir that are not inside a 'ts' directory
    loose_files = [p for p in sim_dir.rglob("*.ts") if p.is_file() and p.parent.name != "ts"]
    if not loose_files:
        return []
        
    organized_paths = []
    for p in loose_files:
        new_name = generate_new_filename(p)
        exp_rel = get_experiment_relative_path(p)
        chunk_idx = get_chunk_index(p.name) or get_chunk_index(Path(new_name).stem)
        
        exp_dir = sim_dir / exp_rel
        exp_dir.mkdir(parents=True, exist_ok=True)
        (exp_dir / "ts").mkdir(parents=True, exist_ok=True)
        (exp_dir / "simulated").mkdir(parents=True, exist_ok=True)
        (exp_dir / "concat").mkdir(parents=True, exist_ok=True)
        
        final_name = f"{chunk_idx}.ts" if chunk_idx else new_name
        if not final_name.endswith(".ts"):
            final_name = f"{final_name}.ts"
        final_path = exp_dir / "ts" / final_name
        
        if p != final_path:
            if final_path.exists():
                print(f"[Organize] Overwriting {final_path} with loose file {p.relative_to(sim_dir)} -> {exp_rel}/ts/{final_name}", file=sys.stderr)
                final_path.unlink()
            else:
                print(f"[Organize] Moving {p.relative_to(sim_dir)} -> {exp_rel}/ts/{final_name}", file=sys.stderr)
            p.rename(final_path)
        organized_paths.append(final_path)

    return organized_paths

def organize_loose_gtrees_files(sim_dir: Path) -> list[Path]:
    """
    Moves every loose .gtrees file under sim_dir (top level included; not already inside a
    gtrees/ dir, and not a concat/ or genetrees/ output) to <exp>/gtrees/<chunk_idx>.gtrees. Refuses to
    overwrite a different existing file; an identical one is left in place and reported.
    Returns only the newly moved paths, so experiments whose gene trees already exist aren't
    re-concatenated.
    """
    if not sim_dir.exists():
        return []

    loose_files = sorted(
        p for p in sim_dir.rglob(f"*{GTREES_EXT}")
        if p.is_file() and p.parent.name not in ("gtrees", "concat", "genetrees")
    )
    joined = [p for p in loose_files if is_joined_gtrees_file(p)]
    if joined:
        print(f"[Organize] Skipping {len(joined)} pre-joined gene tree file(s) (*-concat-*{GTREES_EXT}), left in place", file=sys.stderr)
        loose_files = [p for p in loose_files if p not in joined]
    if not loose_files:
        return []

    # Plan every move first, so a conflict stops the run before anything is moved.
    moves = []
    targets = {}
    conflicts = []
    for p in loose_files:
        exp_rel, chunk_idx = classify_gtrees_file(p)
        if not chunk_idx:
            sys.exit(f"Error: Could not determine the chunk index of gene tree file '{p}'.")
        final_path = sim_dir / exp_rel / "gtrees" / f"{chunk_idx}{GTREES_EXT}"
        if final_path in targets:
            conflicts.append(f"{p} and {targets[final_path]} both map to {final_path}")
        elif final_path.exists() and not filecmp.cmp(p, final_path, shallow=False):
            conflicts.append(f"{p} -> {final_path} (a different file already exists there)")
        targets[final_path] = p
        moves.append((p, final_path, exp_rel))
    if conflicts:
        sys.exit("Error: Refusing to overwrite gene tree files:\n  " + "\n  ".join(conflicts))

    organized_paths = []
    for p, final_path, exp_rel in moves:
        if final_path.exists():
            print(f"[Organize] Identical file already at {exp_rel}/gtrees/{final_path.name}, leaving loose copy {p.relative_to(sim_dir)} in place", file=sys.stderr)
            continue
        final_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"[Organize] Moving {p.relative_to(sim_dir)} -> {exp_rel}/gtrees/{final_path.name}", file=sys.stderr)
        p.rename(final_path)
        organized_paths.append(final_path)

    return organized_paths

def main():
    parser = argparse.ArgumentParser(description="Organize new .ts files in simulations dir.")
    parser.add_argument("--dir", type=Path, help="Simulations dir")
    args = parser.parse_args()
    
    sim_dir = args.dir
    if not sim_dir.exists():
        sys.exit(f"Error: {sim_dir} does not exist.")
        
    organize_loose_ts_files(sim_dir)

if __name__ == "__main__":
    main()
