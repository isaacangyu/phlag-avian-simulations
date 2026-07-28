#!/usr/bin/env python3
"""
Utility module forwarding to simulate.py for branch extraction and filename parsing.
"""

from __future__ import annotations
import sys
from pathlib import Path

# Ensure parent directory is in sys.path
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from simulate import (
    KNOWN_CHANGES_MAP,
    parse_filename,
    extract_branch_name,
    generate_new_filename,
    rename_all_files,
)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Extract branch name from simulation filename or batch rename files to new format."
    )
    parser.add_argument("input_paths", nargs="*", help="File paths or filenames to extract branch from or rename.")
    parser.add_argument("-b", "--branch-only", action="store_true", help="Output only extracted branch name.")
    parser.add_argument("-r", "--rename", action="store_true", help="Batch rename files in given directories.")
    parser.add_argument("--dry-run", action="store_true", help="Show proposed renames without changing files.")

    args = parser.parse_args()

    if args.rename:
        from concat_sequences import get_large_files_dir
        dirs = [Path(p) for p in args.input_paths] if args.input_paths else [get_large_files_dir()]
        for d in dirs:
            print(f"Renaming files in {d}...", file=sys.stderr)
            renamed = rename_all_files(d, dry_run=args.dry_run)
            print(f"Completed {len(renamed)} renames in {d}.", file=sys.stderr)
    elif args.input_paths:
        for p in args.input_paths:
            if args.branch_only:
                bname = extract_branch_name(p)
                print(f"{p}: {bname}")
            else:
                meta = parse_filename(p)
                print(f"File: {p}")
                print(f"  Branch : {meta['branch']}")
                print(f"  Prefix : {meta['prefix']}")
                print(f"  Change : {meta['change']}")

if __name__ == "__main__":
    main()
