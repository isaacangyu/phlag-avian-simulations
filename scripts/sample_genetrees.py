import sys
import random
import tskit


def sample_genetrees(inputdir, rep, ns, nt):
    sampled_gtrees_l = []
    for i in range(ns):
        ts = tskit.load(f"{inputdir}/{rep}-{i+1}.ts")
        table = ts.dump_tables()
        nd_to_pop = {}
        for name, nd in enumerate(table.nodes):
            nd_to_pop[name] = table.populations[nd.population].metadata["name"]
        for j in range(0, ts.get_num_trees(), int(ts.get_num_trees() / nt * ns)):
            sampled_gtrees_l.append(ts.at(j).as_newick(node_labels=nd_to_pop))
    return sampled_gtrees_l


def save_genetrees(gtrees_l, filepath):
    with open(filepath, "w") as f:
        f.write("\n".join(gtrees_l))


if __name__ == "__main__":
    inputdir = sys.argv[1]
    name = sys.argv[2]
    outputpath = sys.argv[3]
    num_blocks = sys.argv[4]
    num_genetrees = sys.argv[5]

    save_genetrees(
        sample_genetrees(
            inputdir,
            name,
            num_blocks,
            num_genetrees,
        ),
        outputpath
    )
