#!/usr/bin/env python3
"""
Organizes loose .ts files in the simulations directory into their respective experiment folders.
Creates ts/, simulated/, and concat/ subdirectories, and renames .ts files to <chunk_idx>.ts inside the ts/ folder.
"""

import sys
import argparse
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from concat_sequences import get_experiment_folder_name, get_experiment_relative_path, get_chunk_index
from simulate import generate_new_filename

def migrate_existing_experiments(sim_dir: Path):
    """Migrates experiment directories into the 2-level directory structure (10X/up, 10X/down, etc.)."""
    if not sim_dir.exists():
        return

    reserved_names = {"10X", "recombination", "admixture", "null", "default", "ts", "simulated", "concat", "store", ".", "..", "ape", "up", "down", "low", "high", "iang", "home", "drive2", "root", "none"}

    exp_folders = []
    for p in list(sim_dir.rglob("*")):
        if p.is_dir() and p.name not in reserved_names:
            if any((p / sub).exists() for sub in ("ts", "simulated", "concat")) or any(p.glob("*mapping*.tsv")):
                exp_folders.append(p)

    for item in exp_folders:
        if not item.exists():
            continue
        exp_rel = get_experiment_relative_path(item)
        if not exp_rel or exp_rel in (Path("default"), Path("."), Path("null")) or str(exp_rel).startswith("/"):
            continue
        target_dir = sim_dir / exp_rel
        if item.resolve() == target_dir.resolve():
            continue

        target_dir.parent.mkdir(parents=True, exist_ok=True)
        print(f"[Migrate] Moving experiment folder {item.relative_to(sim_dir)} -> {exp_rel}", file=sys.stderr)
        if target_dir.exists():
            for child in list(item.iterdir()):
                dest = target_dir / child.name
                if not dest.exists():
                    child.rename(dest)
            try:
                item.rmdir()
            except Exception:
                pass
        else:
            try:
                item.rename(target_dir)
            except Exception:
                target_dir.mkdir(parents=True, exist_ok=True)
                for child in list(item.iterdir()):
                    dest = target_dir / child.name
                    if not dest.exists():
                        child.rename(dest)
                try:
                    item.rmdir()
                except Exception:
                    pass

    for legacy_name in ("10X_up", "10X_down", "recombination_up", "recombination_down"):
        for p in list(sim_dir.rglob(legacy_name)):
            if p.is_dir():
                try:
                    p.rmdir()
                except Exception:
                    pass


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
