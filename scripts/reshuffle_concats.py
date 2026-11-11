#!/usr/bin/env python3
"""
One-off: rebuild every alt concat (FASTA and gene trees) with the null chunks in a random order.

For each <exp>/concat/<interval>.fa and .gtrees (null/ is skipped), the old file is renamed to
<interval>_inc.<ext> (null chunks were in increasing order) and a new file is built with
concatenate_pair() / concatenate_gtrees_pair(), which shuffle the null chunks per file (the .fa
and .gtrees of one interval get the same order). Files whose name is not a plain interval (e.g.
a custom pattern) are left alone. Safe to re-run: finished files are skipped.

Usage: python scripts/reshuffle_concats.py [-j JOBS] [--ext .fa|.gtrees] [--dry-run] [exp_dir ...]
"""

from __future__ import annotations

import sys
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from concat_sequences import get_output_dir, concatenate_pair, concatenate_gtrees_pair, interval_label, parse_interval_pattern, GTREES_EXT

INC_SUFFIX = "_inc"


BUILDERS = {".fa": concatenate_pair, GTREES_EXT: concatenate_gtrees_pair}


def interval_stems(concat_dir: Path, ext: str) -> list[str]:
    """Interval stems of <ext> files in concat_dir (old or new)."""
    stems = set()
    for f in concat_dir.glob(f"*{ext}"):
        stem = f.stem.removesuffix(INC_SUFFIX)
        # Only plain interval names (e.g. '40-60', '70-80,85-100') map back to a concat call.
        if parse_interval_pattern(stem) and interval_label(stem) == stem:
            stems.add(stem)
    return sorted(stems)


def rebuild_experiment(output_dir: Path, exp_dir: Path, exts: list[str], dry_run: bool) -> list[str]:
    log = []
    concat_dir = exp_dir / "concat"
    for ext in exts:
        for stem in interval_stems(concat_dir, ext):
            new_f = concat_dir / f"{stem}{ext}"
            inc_f = concat_dir / f"{stem}{INC_SUFFIX}{ext}"
            if not inc_f.exists():
                log.append(f"rename {new_f} -> {inc_f.name}")
                if not dry_run:
                    new_f.rename(inc_f)
            if not new_f.exists():
                log.append(f"build  {new_f}")
                if not dry_run:
                    BUILDERS[ext](Path("null"), exp_dir, output_dir, interval=stem)
    return log


def find_experiments(output_dir: Path) -> list[Path]:
    return sorted(
        d.parent for d in output_dir.glob("*/*/*/concat")
        if d.is_dir()
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Rename old concat files to *_inc.<ext> and rebuild them with shuffled null chunks.")
    parser.add_argument("exp_dirs", nargs="*", help="Experiment dirs relative to the simulations dir (default: all).")
    parser.add_argument("-j", "--jobs", type=int, default=8, help="Experiments built in parallel.")
    parser.add_argument("--ext", choices=list(BUILDERS), action="append", help="Only this file type (repeatable; default: all).")
    parser.add_argument("--dry-run", action="store_true", help="Only print what would happen.")
    args = parser.parse_args(argv)

    output_dir = get_output_dir()
    exps = [output_dir / e for e in args.exp_dirs] if args.exp_dirs else find_experiments(output_dir)
    print(f"[Reshuffle] {len(exps)} experiment(s)", file=sys.stderr)

    failures = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(rebuild_experiment, output_dir, e, args.ext or list(BUILDERS), args.dry_run): e for e in exps}
        for i, fut in enumerate(as_completed(futures), start=1):
            exp = futures[fut].relative_to(output_dir)
            try:
                for line in fut.result():
                    print(f"  {line}")
                print(f"[Reshuffle] ({i}/{len(exps)}) done: {exp}", file=sys.stderr)
            except BaseException as e:  # sys.exit() from the concat code raises SystemExit
                failures.append(f"{exp}: {e}")
                print(f"[Reshuffle] ({i}/{len(exps)}) FAILED: {exp}: {e}", file=sys.stderr)

    if failures:
        print(f"\n[Reshuffle] {len(failures)} failure(s):")
        for f in failures:
            print(f"  {f}")
        sys.exit(1)
    print("\n[Reshuffle] All done.")


if __name__ == "__main__":
    main()
