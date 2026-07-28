import os
import sys
import argparse
import collections

from pathlib import Path


def load_env_vars() -> dict[str, str]:
    env_vars = {}
    current = Path(__file__).resolve().parent
    for _ in range(3):
        env_path = current / ".env"
        if env_path.exists():
            try:
                with open(env_path, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith("#") and "=" in line:
                            key, val = line.split("=", 1)
                            env_vars[key.strip()] = val.strip().strip("'").strip('"')
            except Exception:
                pass
            break
        current = current.parent
    return env_vars


def get_large_dir() -> Path:
    val = os.environ.get("LARGE_DIR") or os.environ.get("LARGE_FILES_DIR")
    if not val:
        env_vars = load_env_vars()
        val = env_vars.get("LARGE_DIR") or env_vars.get("LARGE_FILES_DIR") or env_vars.get("PHLAG_DIR") or os.environ.get("PHLAG_DIR")
    if val:
        return Path(val)
    sys.exit("Error: LARGE_DIR environment variable is not set and could not be found in .env file.")


def get_large_files_dir() -> Path:
    return get_large_dir()


def get_output_dir() -> Path:
    return get_large_dir() / "simulations"


def get_phlag_dir() -> Path:
    val = os.environ.get("PHLAG_DIR")
    if not val:
        env_vars = load_env_vars()
        val = env_vars.get("PHLAG_DIR")
    if val and val != str(get_large_files_dir()):
        return Path(val)
    cand = Path("/home/iang/phlag")
    if cand.exists():
        return cand
    return Path(__file__).resolve().parent.parent

import tskit
import msprime
import dendropy
import daiquiri
import numpy as np
import treeswift as ts


daiquiri.setup(level="DEBUG")
RS=110

def simulate_ts_to_fasta(input_ts: Path, output_path: Path, mutation_rate: float = 1.9086e-08):
    np.random.seed(RS)
    mrate = mutation_rate

    ts = tskit.load(input_ts)
    seq_len = int(ts.sequence_length)

    rate_gamma_shape = np.random.lognormal(1.5, 1)
    lrate_gamma_shape = rate_gamma_shape * 20

    rate_l = mrate * np.random.gamma(
        lrate_gamma_shape, 1 / lrate_gamma_shape, size=seq_len
    )
    lastp = 0
    while lastp < seq_len:
        rlen = np.random.geometric(1 / 1e4)
        curp = min(lastp + rlen + 1, seq_len)
        rate_l[lastp:curp] *= np.random.gamma(rate_gamma_shape, 1 / rate_gamma_shape)
        lastp = curp
    rate_map = msprime.RateMap(position=range(seq_len + 1), rate=rate_l)

    eq_freq = np.random.dirichlet([36, 26, 28, 32])
    model = msprime.GTR(
        relative_rates=np.random.dirichlet([5, 15, 3, 6, 16, 5]),
        equilibrium_frequencies=eq_freq,
    )
    mts = msprime.sim_mutations(ts, rate=rate_map, model=model, random_seed=RS)

    dseq = np.random.choice(["A", "C", "G", "T"], seq_len, p=eq_freq)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for seq, pi in zip(mts.alignments(), mts.samples()):
            pop = mts.populations()[mts.nodes()[pi].population].metadata["name"]
            f.write(f">{pop}\n")
            nseq_l = [dseq[ix] if c == "N" else c for ix, c in enumerate(seq)]
            f.write("".join(nseq_l) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_ts", type=Path, nargs="+", help="Input tree sequence object(s)."
    )
    parser.add_argument(
        "-o",
        "--output-fasta",
        type=Path,
        required=False,
        help="Output FASTA alignment (only applied if a single input_ts is provided).",
    )
    parser.add_argument(
        "-m",
        "--mutation-rate",
        type=float,
        required=False,
        default=1.9086e-08,
        help="Mutation rate.",
    )
    args = parser.parse_args()

    output_dir = get_output_dir()
    from concat_sequences import get_experiment_folder_name

    for ts_path in args.input_ts:
        exp_name = get_experiment_folder_name(ts_path)
        output_subdir = output_dir / exp_name / "simulated"
        output_subdir.mkdir(parents=True, exist_ok=True)

        if args.output_fasta is None or len(args.input_ts) > 1:
            output_path = output_subdir / f"{ts_path.stem}.fa"
        elif not args.output_fasta.is_absolute():
            output_path = output_subdir / args.output_fasta.name
        else:
            output_path = args.output_fasta

        simulate_ts_to_fasta(ts_path, output_path, args.mutation_rate)

