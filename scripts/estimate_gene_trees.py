#!/usr/bin/env python3
"""
Estimates gene trees straight from the per-chunk simulated FASTAs (not from a concat MSA).

  1. Null trees (computed once, shared by every experiment): null chunks 1..12 are treated as one
     sequence. From a seeded random start offset in [0, WINDOW_BP + SKIP_BP), a WINDOW_BP window
     is taken every WINDOW_BP + SKIP_BP bp until the sequence ends, and IQ-TREE 2 infers one tree
     per window. Output: null/genetrees/null.gtrees (one Newick tree per line, genome order).
  2. Alt trees: same sampling over an experiment's alt chunks 1..3, with its own start offset.
     Output: <exp_dir>/genetrees/alt.gtrees.
  3. Insertion (same insert-and-shift as the MSA concat): the first alt trees are inserted at a
     centered interval (e.g. p=0.1 -> 45-55%), later null trees shift right and the null tail is
     trimmed so the tree count stays the same. Output: <exp_dir>/genetrees/<interval_label>.gtrees, named like
     the matching concat/<interval_label>.fa, plus a <interval_label>.sources.tsv sidecar that
     says where each line came from.

Per-window trees are kept in genetrees/<name>.work/ until every window is done, so an
interrupted run picks up where it left off. Existing non-empty outputs are never rebuilt.
"""

from __future__ import annotations

import os
import sys
import random
import shutil
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from concat_sequences import get_output_dir, parse_interval_pattern, parse_proportion, interval_label

WINDOW_BP = 500
SKIP_BP = 8000
STRIDE_BP = WINDOW_BP + SKIP_BP
DEFAULT_SEED = 1
NULL_CHUNKS = list(range(1, 13))
ALT_CHUNKS = [1, 2, 3]
IQTREE_BIN = "iqtree2"


def centered_interval(proportion: float) -> str:
    """Same as simulate.proportion_to_centered_interval (kept here so this module doesn't pull in tskit)."""
    half_width = (proportion * 100) / 2
    return f"{50 - half_width:g}-{50 + half_width:g}"


def draw_offsets(seed: int) -> tuple[int, int]:
    """(null_offset, alt_offset): two separate draws from one seeded RNG, so both are reproducible."""
    rng = random.Random(seed)
    return rng.randrange(STRIDE_BP), rng.randrange(STRIDE_BP)


def window_starts(total_bp: int, offset: int) -> list[int]:
    """Start positions of every full WINDOW_BP window, one every STRIDE_BP bp from offset."""
    return list(range(offset, total_bp - WINDOW_BP + 1, STRIDE_BP))


def _read_chunk(path: Path) -> dict[str, str]:
    """Reads a one-line-per-sequence chunk FASTA."""
    seqs = {}
    header = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                seqs[header] = ""
            elif header is not None:
                seqs[header] += line
    return seqs


def _chunk_length(path: Path) -> int:
    """Length of the first sequence in a chunk FASTA, read without loading the whole file."""
    with open(path) as f:
        f.readline()
        return len(f.readline().strip())


def _chunk_paths(exp_dir: Path, indices: list[int], required: bool) -> list[Path]:
    paths = []
    for idx in indices:
        cf = exp_dir / "simulated" / f"{idx}.fa"
        if cf.is_file():
            paths.append(cf)
        elif required:
            sys.exit(f"Error: Chunk FASTA does not exist: {cf}")
        else:
            print(f"[Gene Trees] Warning: alt chunk {idx} missing, skipping: {cf}", file=sys.stderr)
    if not paths:
        sys.exit(f"Error: No chunk FASTAs {indices} found under {exp_dir / 'simulated'}")
    return paths


def _write_windows(chunk_paths: list[Path], starts: list[int], out_files: list[Path]) -> None:
    """
    Cuts the given windows (global coordinates over the chunks joined in order) out of the chunk
    FASTAs and writes each one as its own FASTA. Chunks are read one at a time; windows that
    cross a chunk boundary are stitched from both sides.
    """
    headers: list[str] | None = None
    parts: list[dict[str, list[str]]] = [{} for _ in starts]
    chunk_start = 0
    for cf in chunk_paths:
        print(f"[Gene Trees] Reading {cf}...", file=sys.stderr)
        seqs = _read_chunk(cf)
        if headers is None:
            headers = list(seqs)
        elif set(seqs) != set(headers):
            sys.exit(f"Error: {cf} has different taxa than {chunk_paths[0]}")
        chunk_end = chunk_start + len(seqs[headers[0]])
        for w, st in enumerate(starts):
            lo, hi = max(st, chunk_start), min(st + WINDOW_BP, chunk_end)
            if lo < hi:
                for h in headers:
                    parts[w].setdefault(h, []).append(seqs[h][lo - chunk_start:hi - chunk_start])
        chunk_start = chunk_end
        del seqs

    for w, out_fa in enumerate(out_files):
        with open(out_fa, "w") as f:
            for h in headers:
                f.write(f">{h}\n{''.join(parts[w][h])}\n")


def _run_iqtree(locus_fa: Path) -> Path:
    """Runs IQ-TREE on one window, keeps only its .treefile, and returns that path."""
    prefix = locus_fa.with_suffix("")
    result = subprocess.run(
        [IQTREE_BIN, "-s", str(locus_fa), "-m", "GTR+G", "--fast", "-T", "1",
         "--prefix", str(prefix), "-seed", "1", "-quiet", "-redo"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    treefile = prefix.with_suffix(".treefile")
    if result.returncode != 0 or not treefile.exists():
        raise RuntimeError(f"IQ-TREE failed on {locus_fa} (exit {result.returncode}): {result.stderr.strip()}")

    for ext in (".iqtree", ".log", ".ckp.gz", ".bionj", ".mldist", ".uniqueseq.phy"):
        prefix.with_suffix(ext).unlink(missing_ok=True)
    locus_fa.unlink()
    return treefile


def _infer_window_trees(chunk_paths: list[Path], offset: int, out_path: Path, threads: int | None) -> Path:
    """
    Samples windows over chunk_paths (joined in order) from `offset` and writes one IQ-TREE gene
    tree per window to out_path, unless it already exists. Also writes <stem>.windows.tsv
    (tree line -> window start/end in the joined sequence).
    """
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"[Gene Trees] Gene trees already exist: {out_path}", file=sys.stderr)
        return out_path

    if shutil.which(IQTREE_BIN) is None:
        sys.exit(f"Error: '{IQTREE_BIN}' not found on PATH.")

    total_bp = sum(_chunk_length(cf) for cf in chunk_paths)
    starts = window_starts(total_bp, offset)
    if not starts:
        sys.exit(f"Error: {total_bp} bp is too short for one {WINDOW_BP} bp window at offset {offset}.")
    n = len(starts)

    work_dir = out_path.parent / f"{out_path.stem}.work"
    work_dir.mkdir(parents=True, exist_ok=True)
    # A resumed run must use the same windows, or old and new trees would be mixed.
    stamp = f"offset={offset} window={WINDOW_BP} skip={SKIP_BP} total_bp={total_bp}\n"
    stamp_path = work_dir / "windows.txt"
    if stamp_path.exists() and stamp_path.read_text() != stamp:
        sys.exit(f"Error: {work_dir} was started with different windows ({stamp_path.read_text().strip()}); "
                 f"now {stamp.strip()}. Remove {work_dir} to start over.")
    stamp_path.write_text(stamp)

    width = len(str(n))
    treefiles = [work_dir / f"{i:0{width}d}.treefile" for i in range(1, n + 1)]
    todo = [i for i, tf in enumerate(treefiles) if not tf.exists()]

    if todo:
        locus_files = [treefiles[i].with_suffix(".fa") for i in todo]
        _write_windows(chunk_paths, [starts[i] for i in todo], locus_files)

        workers = threads or os.cpu_count() or 1
        print(
            f"[Gene Trees] Inferring {len(todo)}/{n} gene trees ({WINDOW_BP} bp windows every "
            f"{STRIDE_BP} bp from offset {offset}, {total_bp} bp total, {workers} workers) -> {out_path}",
            file=sys.stderr,
        )
        done = n - len(todo)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_run_iqtree, fa) for fa in locus_files]
            for fut in as_completed(futures):
                fut.result()
                done += 1
                if done % 100 == 0 or done == n:
                    print(f"[Gene Trees] {done}/{n} done", file=sys.stderr)

    windows_path = out_path.with_suffix(".windows.tsv")
    with open(windows_path, "w") as f:
        f.write("line\tstart\tend\n")
        for i, st in enumerate(starts, 1):
            f.write(f"{i}\t{st}\t{st + WINDOW_BP}\n")

    tmp_path = out_path.with_suffix(".gtrees.tmp")
    with open(tmp_path, "w") as out:
        for tf in treefiles:
            out.write(tf.read_text().strip() + "\n")
    tmp_path.rename(out_path)
    shutil.rmtree(work_dir)
    print(f"[Gene Trees] Wrote {n} gene trees to {out_path}", file=sys.stderr)
    return out_path


def ensure_null_gene_trees(output_dir: Path, threads: int | None = None, seed: int = DEFAULT_SEED) -> Path:
    """null/genetrees/null.gtrees from null chunks 1..12."""
    null_dir = output_dir / "null"
    chunks = _chunk_paths(null_dir, NULL_CHUNKS, required=True)
    return _infer_window_trees(chunks, draw_offsets(seed)[0], null_dir / "genetrees" / "null.gtrees", threads)


def ensure_alt_gene_trees(exp_dir: Path, threads: int | None = None, seed: int = DEFAULT_SEED) -> Path:
    """<exp_dir>/genetrees/alt.gtrees from the experiment's alt chunks 1..3 (missing ones skipped)."""
    chunks = _chunk_paths(exp_dir, ALT_CHUNKS, required=False)
    return _infer_window_trees(chunks, draw_offsets(seed)[1], exp_dir / "genetrees" / "alt.gtrees", threads)


def to_interval(proportion_or_interval: str | float) -> str:
    """0.1 / '0.1' -> '45-55'; an interval string like '40-60' or '10-20_70-80' is returned as is."""
    s = str(proportion_or_interval).strip()
    if parse_interval_pattern(s) is not None:
        return s
    try:
        prop = parse_proportion(s)
    except (ValueError, ZeroDivisionError):
        sys.exit(f"Error: '{s}' is neither a proportion nor an interval like '40-60'.")
    if prop is None or not 0 < prop <= 1:
        sys.exit(f"Error: Proportion must be in (0, 1], got '{s}'.")
    return centered_interval(prop)


def insert_alt_trees(null_gtrees: Path, alt_gtrees: Path, interval: str, out_dir: Path) -> Path:
    """
    Writes out_dir/<interval_label>.gtrees exactly like the MSA concat, in tree units: alt trees
    (first ones, in order) are inserted at `interval` (percent of the null tree list), null trees
    shift right past each insert, and the null tail is trimmed so the count stays the same.
    Also writes <interval_label>.sources.tsv (line, source, 1-based index into null/alt.gtrees).
    """
    label = interval_label(interval)
    out_path = out_dir / f"{label}.gtrees"
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"[Gene Trees] Gene trees already exist: {out_path}", file=sys.stderr)
        return out_path

    intervals = parse_interval_pattern(interval)
    if any(end > 100 for _, end in intervals):
        sys.exit(f"Error: Gene tree insertion takes percent intervals (0-100), got '{interval}'.")

    null_trees = [t for t in null_gtrees.read_text().splitlines() if t.strip()]
    alt_trees = [t for t in alt_gtrees.read_text().splitlines() if t.strip()]
    n = len(null_trees)

    # Same arithmetic as concatenate_pair(), in tree units instead of bp.
    spans = sorted((int(st / 100.0 * n), int((end - st) / 100.0 * n)) for st, end in intervals)
    need = sum(length for _, length in spans)
    if need > len(alt_trees):
        sys.exit(f"Error: Interval {interval} needs {need} alt trees but {alt_gtrees} has only {len(alt_trees)}.")

    # Insert-and-shift, same loop as concatenate_pair(): null up to each interval start, then the
    # alt block, then null continues from where it stopped; the null tail is trimmed so the total
    # stays n.
    trimmed_null_length = max(0, n - need)
    sources: list[tuple[str, int]] = []
    null_idx = alt_idx = curr_pos = 0
    for st, length in spans:
        null_len = st - curr_pos
        if null_len > 0:
            sources += [("null", i) for i in range(null_idx, null_idx + null_len)]
            null_idx += null_len
            curr_pos += null_len
        sources += [("alt", i) for i in range(alt_idx, alt_idx + length)]
        alt_idx += length
        curr_pos += length
    remaining_null_len = trimmed_null_length - null_idx
    if remaining_null_len > 0:
        sources += [("null", i) for i in range(null_idx, null_idx + remaining_null_len)]

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(".gtrees.tmp")
    with open(tmp_path, "w") as f:
        for src, i in sources:
            f.write((null_trees if src == "null" else alt_trees)[i] + "\n")
    with open(out_path.with_suffix(".sources.tsv"), "w") as f:
        f.write("line\tsource\tindex\n")
        for line, (src, i) in enumerate(sources, 1):
            f.write(f"{line}\t{src}\t{i + 1}\n")
    tmp_path.rename(out_path)

    alt_lines = ", ".join(f"{st + 1}-{st + length}" for st, length in spans if length)
    print(f"[Gene Trees] Wrote {n} gene trees to {out_path} ({need} alt trees at lines {alt_lines or 'none'})", file=sys.stderr)
    return out_path


def ensure_estimated_gene_trees(
    output_dir: Path,
    exp_dir: Path | None,
    proportions_or_interval: str | float | list[str | float] | None = None,
    threads: int | None = None,
    seed: int = DEFAULT_SEED,
) -> list[Path]:
    """
    Ensures the null trees, and (when exp_dir is given) the experiment's alt trees plus one
    inserted tree file per proportion/interval. exp_dir may be absolute or relative to
    output_dir. Returns the paths of every tree file ensured, null first.
    """
    paths = [ensure_null_gene_trees(output_dir, threads, seed)]
    if exp_dir is None:
        return paths

    exp_dir = exp_dir if exp_dir.is_absolute() else output_dir / exp_dir
    if exp_dir.resolve() == (output_dir / "null").resolve():
        return paths
    if not exp_dir.is_dir():
        sys.exit(f"Error: Experiment directory does not exist: {exp_dir}")

    alt_gtrees = ensure_alt_gene_trees(exp_dir, threads, seed)
    paths.append(alt_gtrees)

    if proportions_or_interval is None:
        return paths
    items = proportions_or_interval if isinstance(proportions_or_interval, (list, tuple)) else [proportions_or_interval]
    for item in items:
        paths.append(insert_alt_trees(paths[0], alt_gtrees, to_interval(item), exp_dir / "genetrees"))
    return paths


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=f"Estimate gene trees ({WINDOW_BP} bp windows, {SKIP_BP} bp apart) from simulated chunk FASTAs with IQ-TREE 2."
    )
    parser.add_argument("target", help="'null', or an experiment dir relative to the simulations dir (e.g. 10X/up/N554).")
    parser.add_argument("interval", nargs="*", help="Percent interval(s) to insert alt trees into, e.g. 40-60.")
    parser.add_argument("-p", "--proportion", action="append", default=[], help="Centered alt proportion, e.g. 0.1 -> 45-55. Repeatable.")
    parser.add_argument("-t", "--threads", type=int, default=None, help="Max parallel IQ-TREE jobs (default: all CPUs).")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Seed for the random window start offsets.")
    parser.add_argument("-o", "--output-dir", type=Path, default=None, help="Simulations dir (default: $CONNECTION_DIR/simulations).")
    args = parser.parse_intermixed_args(argv)

    output_dir = args.output_dir or get_output_dir()
    exp_dir = None if args.target.strip("/").lower() == "null" else Path(args.target)
    items = args.proportion + args.interval
    if exp_dir is None and items:
        sys.exit("Error: Proportions/intervals need an experiment target, not 'null'.")

    for p in ensure_estimated_gene_trees(output_dir, exp_dir, items or None, args.threads, args.seed):
        print(p)


if __name__ == "__main__":
    main()
