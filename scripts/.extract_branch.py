#!/usr/bin/env python3
"""
Hidden script to extract branch names, parse simulation filenames,
and rename files into standard format:
    null/alt-(neoaves)-branch_change_(rate-time)-index.ext

Rules:
  - sizes are always 500kb (omitted in standard format)
  - null FASTA files are always 'null-index.fa'
  - prefixes: 'default'/'default-neoaves' -> 'null', 'anomaly_case' -> 'alt'
  - changes: {admixture, recombination_increase, 10X_up, recombination_suppression, 10X_down, recombination_decrease}
  - 'increase'/'decrease' renamed to 'up'/'down' (e.g. recombination_increase -> recombination_up)
  - admixture has rate and time in filename
"""

from __future__ import annotations

import re
import sys
import argparse
from pathlib import Path

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
    Note: For null FASTA files (.fa/.fasta), always outputs 'null-neoaves-index.fa' (no 'null-index.fa').
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
    """
    Batch renames files in target_dir according to the new naming scheme.
    """
    target_path = Path(target_dir)
    if not target_path.exists():
        print(f"Directory does not exist: {target_path}", file=sys.stderr)
        return []

    renamed = []
    files = [p for p in target_path.rglob("*") if p.is_file() and not p.name.startswith(".")]

    for old_path in sorted(files):
        if old_path.name == "ape.fa":
            continue
        new_filename = generate_new_filename(old_path)
        new_path = old_path.parent / new_filename
        
        if old_path != new_path:
            print(f"Renaming: {old_path.name} -> {new_filename}", file=sys.stderr)
            if not dry_run:
                if new_path.exists() and new_path != old_path:
                    old_path.unlink()
                else:
                    old_path.rename(new_path)
            renamed.append((old_path, new_path))

    return renamed


def main():
    parser = argparse.ArgumentParser(
        description="Extract branch name from simulation filename or batch rename files to new format."
    )
    parser.add_argument("input_paths", nargs="*", help="File paths or filenames to extract branch from or rename.")
    parser.add_argument("-b", "--branch-only", action="store_true", help="Output only extracted branch name.")
    parser.add_argument("-r", "--rename", action="store_true", help="Batch rename files in given directories.")
    parser.add_argument("--dry-run", action="store_true", help="Show proposed renames without changing files.")

    args = parser.parse_args()

    if args.rename:
        dirs = [Path(p) for p in args.input_paths] if args.input_paths else [Path("store")]
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
                new_f = generate_new_filename(p)
                print(f"File: {p}")
                print(f"  Branch : {meta['branch']}")
                print(f"  Change : {meta['change']}")
                print(f"  Rate   : {meta['rate']}")
                print(f"  Time   : {meta['time']}")
                print(f"  NewName: {new_f}\n")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
