#!/usr/bin/env python3
"""
CLI script to extract experiment folder names from input strings or CSV files.
"""

from __future__ import annotations

import sys
import csv
import argparse
import re
from pathlib import Path

# Ensure scripts directory is in sys.path
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from concat_sequences import get_experiment_folder_name


def process_string(input_str: str) -> str:
    return get_experiment_folder_name(input_str)


def process_csv(
    csv_path: Path,
    col_name: str | None = None,
    output_path: Path | None = None,
    in_place: bool = False,
):
    if not csv_path.exists():
        sys.exit(f"Error: CSV file '{csv_path}' does not exist.")

    if in_place:
        output_path = csv_path
    elif output_path is None:
        output_path = csv_path.parent / f"{csv_path.stem}_experiment_folders.csv"

    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            sys.exit(f"Error: CSV file '{csv_path}' is empty.")

        target_col_idx = 0
        if col_name:
            if col_name in header:
                target_col_idx = header.index(col_name)
            else:
                try:
                    target_col_idx = int(col_name)
                except ValueError:
                    sys.exit(f"Error: Column '{col_name}' not found in CSV header: {header}")
        else:
            # Auto-detect plausible column name
            lowered = [c.lower().strip() for c in header]
            possible_names = ["filename", "file", "name", "input", "path", "sequence", "string"]
            found = False
            for p in possible_names:
                if p in lowered:
                    target_col_idx = lowered.index(p)
                    found = True
                    break
            if not found:
                target_col_idx = 0

        # Collect raw rows and calculate max underscore parts
        raw_rows = list(reader)
        max_parts = 0
        parsed_data = []

        for row in raw_rows:
            if not row:
                continue
            input_val = row[target_col_idx] if target_col_idx < len(row) else ""
            exp_folder = process_string(input_val) if input_val else ""
            parts = re.split(r'[_.-]', input_val) if input_val else []
            max_parts = max(max_parts, len(parts))
            parsed_data.append((row, exp_folder, parts))

        # Determine output header
        out_header = list(header)
        exp_col_name = "experiment_folder"
        if exp_col_name not in out_header:
            out_header.append(exp_col_name)
        
        part_col_start = len(out_header)
        for i in range(1, max_parts + 1):
            part_col_name = f"part_{i}"
            if part_col_name not in out_header:
                out_header.append(part_col_name)

        out_rows = []
        for row, exp_folder, parts in parsed_data:
            new_row = list(row)
            while len(new_row) < len(out_header):
                new_row.append("")
            
            exp_col_idx = out_header.index(exp_col_name)
            new_row[exp_col_idx] = exp_folder

            for i, p_val in enumerate(parts, start=1):
                part_col_name = f"part_{i}"
                p_col_idx = out_header.index(part_col_name)
                new_row[p_col_idx] = p_val

            out_rows.append(new_row)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(out_header)
        writer.writerows(out_rows)
    print(f"Wrote results to '{output_path}'.", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract experiment folder names from string inputs or CSV files."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        type=str,
        help="Input string(s) or filename(s) to extract experiment folder names from.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Path to CSV file containing input strings.",
    )
    parser.add_argument(
        "--column", "-c",
        type=str,
        default=None,
        help="Column name or index in CSV to extract from (defaults to auto-detect or first column).",
    )
    parser.add_argument(
        "-i", "--in-place", "--inplace",
        action="store_true",
        dest="in_place",
        help="Overwrite input CSV file in place.",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output CSV path when using --csv (defaults to <same_dir>/<stem>_experiment_folders.csv).",
    )

    args = parser.parse_args(argv)

    if args.csv:
        process_csv(args.csv, col_name=args.column, output_path=args.output, in_place=args.in_place)
    elif args.inputs:
        for inp in args.inputs:
            exp_folder = process_string(inp)
            print(f"{inp}\t{exp_folder}")
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
