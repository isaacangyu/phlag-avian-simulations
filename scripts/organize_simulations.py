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

from concat_sequences import get_experiment_folder_name, get_chunk_index
from simulate import generate_new_filename

def organize_loose_ts_files(sim_dir: Path):
    if not sim_dir.exists():
        return
        
    loose_files = [p for p in sim_dir.glob("*.ts") if p.is_file()]
    if not loose_files:
        return
        
    for p in loose_files:
        new_name = generate_new_filename(p)
        exp_name = get_experiment_folder_name(new_name)
        chunk_idx = get_chunk_index(Path(new_name).stem)
        
        exp_dir = sim_dir / exp_name
        exp_dir.mkdir(parents=True, exist_ok=True)
        (exp_dir / "ts").mkdir(parents=True, exist_ok=True)
        (exp_dir / "simulated").mkdir(parents=True, exist_ok=True)
        (exp_dir / "concat").mkdir(parents=True, exist_ok=True)
        
        final_name = f"{chunk_idx}.ts" if chunk_idx else new_name
        final_path = exp_dir / "ts" / final_name
        
        print(f"[Organize] Moving {p.name} -> {exp_name}/ts/{final_name}", file=sys.stderr)
        p.rename(final_path)

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
