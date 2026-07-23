import os
import sys
import argparse
import collections

from pathlib import Path


def get_output_dir() -> Path:
    phlag_dir = os.environ.get("PHLAG_DIR")
    if not phlag_dir:
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
                                if key.strip() == "PHLAG_DIR":
                                    phlag_dir = val.strip().strip("'").strip('"')
                                    break
                except Exception:
                    pass
            if phlag_dir:
                break
            current = current.parent

    if phlag_dir:
        return Path(phlag_dir)
    else:
        repo_root = Path(__file__).resolve().parent.parent
        return repo_root / "data"

import tskit
import msprime
import dendropy
import daiquiri
import numpy as np
import treeswift as ts


daiquiri.setup(level="DEBUG")
RS=110

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input-ts", type=Path, required=True, help="Input tree sequence object."
    )
    parser.add_argument(
        "-o",
        "--output-fasta",
        type=Path,
        required=False,
        help="Output FASTA alignment.",
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
    (output_dir / "msa").mkdir(parents=True, exist_ok=True)

    if args.output_fasta is None:
        output_path = output_dir / "msa" / f"{args.input_ts.stem}.fa"
    elif not args.output_fasta.is_absolute():
        output_path = output_dir / "msa" / args.output_fasta.name
    else:
        output_path = args.output_fasta

    np.random.seed(RS)

    mrate = args.mutation_rate  # 1.908590331970381e-08  # mean mutation rate

    ts = tskit.load(args.input_ts)
    seq_len = int(ts.sequence_length)

    rate_gamma_shape = np.random.lognormal(1.5, 1)
    lrate_gamma_shape = rate_gamma_shape * 20

    rate_l = mrate * np.random.gamma(
        lrate_gamma_shape, 1 / lrate_gamma_shape, size=seq_len
    )
    lastp = 0
    while lastp < seq_len:
        rlen = np.random.geometric(1 / 1e4)  # mean local length is 1e4
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
    with open(output_path, "w") as f:
        for seq, pi in zip(mts.alignments(), mts.samples()):
            pop = mts.populations()[mts.nodes()[pi].population].metadata["name"]
            f.write(f">{pop}\n")
            nseq_l = [dseq[ix] if c == "N" else c for ix, c in enumerate(seq)]
            f.write("".join(nseq_l) + "\n")
