# Avian gene tree simulations used to benchmark Phlag

This repository simulates avian sequence data to benchmark [Phlag](https://github.com/isaacangyu/phlag)'s
phylogenetic anomaly detection. 

## Environment setup

```bash
micromamba create -n phlag-avian-simulations python=3.9
micromamba activate phlag-avian-simulations
pip install .
```

## Getting started

Drag your `.ts` tree sequence files into `simulations/` (anywhere under it), then run:

```bash
simulate -r
simulate
```

## Output directory layout

Directory names extract the keywords from avian tree sequence filenames (deposited on
[Dryad](https://datadryad.org)):

```
null/alt-branch_change_(rate-time)-index.ts
```
File structure:

```
simulations/
  null/                              ts/, simulated/, concat/
  10X/{up,down}/<branch>/            ts/, simulated/, concat/
  recombination/{up,down}/<branch>/  ts/, simulated/, concat/
  admixture/{low,high}/<branch>/     ts/, simulated/, concat/
```

- `ts/` — raw tree sequence files (`<chunk_idx>.ts`, chunks 1-12)
- `simulated/` — sequences simulated to FASTA (`<chunk_idx>.fa`), Neoaves-filtered for null/node chunks
- `concat/` — final concatenated 6Mbp FASTA files (null + inserted alt proportion)

`admixture` experiments are split into `low`/`high` by divergence time (`> 4,000,000` generations
is `high`); `10X` and `recombination` experiments are split by `up`/`down` direction. Classification
is derived solely from `parse_filename()`'s structured metadata — there is no path-substring
guessing or fallback location.

## Scripts

### `simulate.py`
Main simulation pipeline entry point (`simulate` console script via `pyproject.toml`). Exactly one
of four modes is selected by the positional `target` argument:

```bash
# 1. No target: complete every organized experiment directory (simulate, map, concat sweep). Skips null.
simulate

# 2. -r/--run-all: organize loose .ts files into experiment dirs, then complete every experiment, including null.
simulate -r

# 3. 'null': simulate/filter/concat only the null experiment.
simulate null

# 4. <relative_dir>: simulate/map/concat one existing experiment directory against null.
simulate 10X/up/N554
```

Optional flags (modes 3 and 4 only):

| Flag | Description |
|---|---|
| `interval` (2nd positional) | Explicit insertion interval, e.g. `40-60`. Omit for the default six-proportion centered sweep (`0.02, 0.05, 0.1, 0.15, 0.20, 0.25`). |
| `-p`, `--prop`, `--proportion` | Proportion of alternative sequence (e.g. `0.25`, `3/12`). |
| `--np`, `--null-pattern` | Null chunk order pattern (default `n1n2n3n4n5n6n7n8n9n10n11n12`). |
| `--ap`, `--alt-pattern` | Alt chunk order pattern (default `a1a2a3`). |
| `-mrate`, `--mutation-rate` | Mutation rate for sequence simulation (default `1.9086e-08`). |

Species-tree mapping (`P1`/`P2`/`P3`/`Po` for `dstar`) is generated automatically.
The Neoaves tree (`speciestree/main-neoaves-numgen.nwk`) is tried first, then the full Avian tree
(`speciestree/main-avian-numgen.nwk`).

Run `simulate -h` for the full flag reference.

## Environment Configuration

All scripts require `$CONNECTION_DIR`. This repository reads from and writes to `$CONNECTION_DIR/simulations/`; Phlag reads from `$CONNECTION_DIR/simulations/` and writes to `$CONNECTION_DIR/phlag/`. Set it via a `.env` file in the repo root:

```env
CONNECTION_DIR=/path/to/shared_directory
```