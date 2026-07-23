import sys
import collections
import argparse
from pathlib import Path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input-fasta", type=Path, required=True, help="Input tree sequence object."
    )

    parser.add_argument(
       "--group-size", type=int, required=False, default=8, help="Number of consecutive genes to group and choose one from."
    )
    parser.add_argument(
       "--group-index", type=int, required=False, default=1, help="Index of the group to extract genes from."
    )
    parser.add_argument(
       "--gene-length", type=int, required=False, default=500, help="Length of each gene."
    )
    parser.add_argument(
       "--gene-index", type=int, required=False, default=2, help="Index of the gene that will be chosen from each group."
    )
    args = parser.parse_args()
    assert(args.group_size >= args.gene_index)
    group_length = args.gene_length * args.group_size

    with open(args.input_fasta, 'r') as f:
        for ix, line in enumerate(f):
            if line[0] == ">":
                print(line.strip(), end="\n")
            else:
                ix_s = group_length * (args.group_index - 1) + (args.gene_index - 1) * args.gene_length
                assert((ix_s + args.gene_length) <= len(line))
                print(line.strip()[ix_s : ix_s +  args.gene_length], end="\n")
