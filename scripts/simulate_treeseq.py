import sys
import argparse
import collections
from pathlib import Path

import tskit
import demes
import msprime
import dendropy
import daiquiri
import treeswift as ts
from dendropy.simulate import treesim


daiquiri.setup(level="DEBUG")


def scale_time(popn, time):
    return (
        popn_to_stime[popn]
        + (time - popn_to_time[popn]) / popn_to_dur[popn] * popn_to_sdur[popn]
    )
    # return popn_to_stime[popn] + (time - popn_to_time[popn]) * pop_to_rmult[popn]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-i", "--input-tree", type=Path, required=True, help="Input species tree."
    )
    parser.add_argument(
        "-d",
        "--input-demography",
        type=Path,
        help="Input demes demography, overrides the species tree.",
    )
    parser.add_argument(
        "-p",
        "--pop-sizes",
        required=True,
        type=Path,
        help="Population to size mapping.",
    )
    parser.add_argument(
        "-s",
        "--sub-units",
        required=True,
        type=Path,
        help="Population to substitution units mapping.",
    )
    parser.add_argument(
        "-g",
        "--num-generations",
        required=True,
        type=Path,
        help="Population to number of generations mapping.",
    )
    parser.add_argument(
        "-o", "--output-path", required=True, type=Path, help="Output path."
    )
    parser.add_argument(
        "-l", "--seq-len", required=True, type=int, help="Sequence length."
    )
    parser.add_argument(
        "-r",
        "--recombination-rate",
        required=False,
        default=5e-9,
        type=float,
        help="Recombination rate.",
    )
    parser.add_argument(
        "--recombination-rate-change",
        required=False,
        default=1e-5,
        type=float,
        help="Recombination rate change.",
    )
    parser.add_argument(
        "--parameter-change-population",
        required=False,
        type=str,
        help="Population for which the parameter changes will take effect.",
    )
    parser.add_argument(
        "--parameter-change-time",
        required=False,
        default=0,
        type=float,
        help="Parameter change time.",
    )
    parser.add_argument(
        "--parameter-change-duration",
        required=False,
        default=0,
        type=float,
        help="Parameter change duration.",
    )
    parser.add_argument("--seed", required=True, type=int, help="Random seed.")
    parser.add_argument(
        "--root-population",
        required=False,
        default="root",
        type=str,
        help="Name of the root population.",
    )
    args = parser.parse_args()

    if args.parameter_change_time and not args.parameter_change_duration:
        parser.error("--parameter-change-time requires --parameter-change-duration to be specified.")
    if args.parameter_change_time and args.parameter_change_population:
        parser.error("--parameter-change-time and --parameter-change-population cannot be specified together.")
    if args.parameter_change_duration and args.parameter_change_population:
        parser.error("--parameter-change-duration and --parameter-change-population cannot be specified together.")

    print("--input-tree", args.input_tree, file=sys.stderr)
    print("--input-demography", args.input_demography, file=sys.stderr)
    print("--pop-sizes", args.pop_sizes, file=sys.stderr)
    print("--sub-units", args.sub_units, file=sys.stderr)
    print("--num-generations", args.num_generations, file=sys.stderr)
    print("--output-path", args.output_path, file=sys.stderr)
    print("--seq-len", args.seq_len, file=sys.stderr)
    print("--recombination-rate", args.recombination_rate, file=sys.stderr)
    print("--recombination-rate-change", args.recombination_rate, file=sys.stderr)
    print("--parameter-change-population", args.parameter_change_population, file=sys.stderr)
    print("--parameter-change-time", args.parameter_change_time, file=sys.stderr)
    print(
        "--parameter-change-duration", args.parameter_change_duration, file=sys.stderr
    )
    print("--seed", args.seed, file=sys.stderr)
    print("--root-population", args.root_population, file=sys.stderr)

    t = dendropy.Tree.get(
        path=args.input_tree, schema="newick", preserve_underscores=True
    )

    samples = []
    height = t.max_distance_from_root()
    for ix, node in enumerate(t.leaf_nodes()):
        # st = height - node.distance_from_root()
        # node.edge_length += st
        st = 0
        samples.append(
            msprime.SampleSet(1, ploidy=1, time=st, population=node.taxon.label)
        )

    tmain = ts.read_tree_newick(args.input_tree)
    label_to_node = tmain.label_to_node(selection="all")

    initial_size = {}
    with open(args.pop_sizes, "r") as f:
        initial_size = dict(map(lambda x: x.strip().split("\t"), f.readlines()[1:]))
    initial_size = {k: float(v) / 2 for k, v in initial_size.items()}

    if args.input_demography is not None:
        graph = demes.load(args.input_demography)
        demography = msprime.Demography.from_demes(graph)
    else:
        tree = t.as_string(
            schema="newick", suppress_rooting=True, unquoted_underscores=True
        )
        demography = msprime.Demography.from_species_tree(tree, initial_size)

    events = [e for e in demography.events]
    rmult_ignore = set()
    for e in events:
        edict = e.asdict()
        if edict.get("source", "") != "":
            nd_event = label_to_node[edict["source"]]
            for nd in nd_event.traverse_bfs(include_self=True):
                rmult_ignore.add(nd[0].get_label())
        if edict.get("dest", "") != "":
            nd_event = label_to_node[edict["dest"]]
            for nd in nd_event.traverse_bfs(include_self=True):
                rmult_ignore.add(nd[0].get_label())

        if (isinstance(edict.get("ancestral", ""), list)):
            for lbl_ancestral in e.ancestral:
                nd_ancestral = label_to_node[lbl_ancestral]
                for nd in nd_ancestral.traverse_bfs(include_self=True):
                    rmult_ignore.add(nd[0].get_label())
            nd_derived = label_to_node[e.derived]
            for nd in nd_derived.traverse_bfs(include_self=True):
                rmult_ignore.add(nd[0].get_label())
    demography.sort_events()

    pop_to_popn = {}
    for ix, pop in enumerate(demography.populations):
        pop_to_popn[ix] = pop.name

    if args.parameter_change_time > 0 and args.parameter_change_duration > 0:
        tseq = msprime.sim_ancestry(
            samples=samples,
            ploidy=2,
            demography=demography,
            recombination_rate=args.recombination_rate,
            sequence_length=args.seq_len,
            random_seed=args.seed,
            end_time=args.parameter_change_time,
        )
        tseq = msprime.sim_ancestry(
            initial_state=tseq,
            ploidy=2,
            demography=demography,
            recombination_rate=args.recombination_rate * args.recombination_rate_change,
            sequence_length=args.seq_len,
            random_seed=args.seed,
            end_time=args.parameter_change_time + args.parameter_change_duration,
        )
        tseq = msprime.sim_ancestry(
            initial_state=tseq,
            ploidy=2,
            demography=demography,
            recombination_rate=args.recombination_rate,
            sequence_length=args.seq_len,
            random_seed=args.seed,
        )
    elif args.parameter_change_population:
        nd_interest = label_to_node[args.parameter_change_population]
        population_interest = [nd.get_label() for nd in tmain.extract_subtree(nd_interest).traverse_leaves()]
        samples_interest = []
        samples_outside = []
        for sample in samples:
            if sample.population in population_interest:
                samples_interest.append(sample)
            else:
                samples_outside.append(sample)
        time_start = int(tmain.extract_subtree(nd_interest).height())
        time_stop = int(tmain.extract_subtree(nd_interest.get_parent()).height())

        tseqP = msprime.sim_ancestry(
            samples=samples_interest,
            ploidy=2,
            demography=demography,
            recombination_rate=args.recombination_rate,
            sequence_length=args.seq_len,
            random_seed=args.seed,
            end_time = time_start
        )
        tseqA = msprime.sim_ancestry(
            initial_state=tseqP,
            ploidy=2,
            demography=demography,
            model=msprime.StandardCoalescent(duration=time_stop-time_start),
            recombination_rate=args.recombination_rate*args.recombination_rate_change,
            sequence_length=args.seq_len,
            random_seed=args.seed,
            # start_time = time_start,
            # end_time = time_stop
        )
        tseqB = msprime.sim_ancestry(
            samples=samples_outside,
            ploidy=2,
            demography=demography,
            recombination_rate=args.recombination_rate,
            sequence_length=args.seq_len,
            random_seed=args.seed,
            end_time = time_stop
        )
        tablesB = tseqB.dump_tables()
        tablesA = tseqA.dump_tables()

        tablesB.edges.child += tseqA.num_nodes
        tablesB.edges.parent += tseqA.num_nodes
        tablesB.nodes.individual += tseqA.num_individuals

        individuals_dict = tablesB.individuals.asdict()
        del individuals_dict["metadata_schema"]
        tablesA.individuals.append_columns(**individuals_dict)
        nodes_dict = tablesB.nodes.asdict()
        del nodes_dict["metadata_schema"]
        tablesA.nodes.append_columns(**nodes_dict)
        edges_dict = tablesB.edges.asdict()
        del edges_dict["metadata_schema"]
        tablesA.edges.append_columns(**edges_dict)
        tablesA.sort()
        tablesA.build_index()
        tseq_joint = tablesA.tree_sequence()

        tseq = msprime.sim_ancestry(
            initial_state=tseq_joint,
            ploidy=2,
            demography=demography,
            recombination_rate=args.recombination_rate,
            sequence_length=args.seq_len,
            random_seed=args.seed,
            start_time = time_stop,
        )
    else:
        tseq = msprime.sim_ancestry(
            samples=samples,
            ploidy=2,
            demography=demography,
            recombination_rate=args.recombination_rate,
            sequence_length=args.seq_len,
            random_seed=args.seed,
        )
    tseq.simplify()
    tables = tseq.dump_tables()

    node_to_popn = {}
    for node, i in enumerate(tables.nodes):
        node_to_popn[node] = tables.populations[i.population].metadata["name"]
    retained_popn_s = set(node_to_popn.values())

    popn_to_sunit = {}
    with open(args.sub_units, "r") as f:
        for ix, l in enumerate(f):
            if ix > 0:
                pop, sunit = l.strip().split("\t")
                popn_to_sunit[pop] = float(sunit)

    popn_to_ngen = {}
    with open(args.num_generations, "r") as f:
        for ix, l in enumerate(f):
            if ix > 0:
                pop, ngen = l.strip().split("\t")
                ngen = eval(ngen)
                if ngen is not None:
                    popn_to_ngen[pop] = float(ngen)

    popn_to_rmult = {}
    popn_to_srate = {}
    mean_srate = 0
    ix = 0
    for pop, ngen in popn_to_ngen.items():
        popn_to_srate[pop] = popn_to_sunit[pop] / ngen
        if pop in retained_popn_s:
            mean_srate += popn_to_srate[pop]
            ix = ix + 1

    mean_srate /= ix
    for pop, srate in popn_to_srate.items():
        popn_to_rmult[pop] = srate / mean_srate
    popn_to_rmult["root"] = 1
    popn_to_srate["root"] = mean_srate
    for label in rmult_ignore:
        popn_to_rmult[label] = 1
        popn_to_srate[label] = mean_srate


    popn_to_stime = {}
    popn_to_sdur = {}
    popn_to_time = {}
    popn_to_dur = {}

    for nd in tmain.traverse_levelorder(internal=True, leaves=True):
        l = nd.get_edge_length()
        if l is None:
            l = 0
        popn_to_time[nd.get_label()] = tmain.extract_subtree(nd).height() - l
        popn_to_dur[nd.get_label()] = l

    for nd in tmain.traverse_postorder(internal=True, leaves=False):
        for nd_child in nd.child_nodes():
            nd_child.set_edge_length(
                nd_child.get_edge_length() * popn_to_rmult[nd_child.get_label()]
            )

    h = tmain.height()
    for nd, dist in tmain.distances_from_root(
        leaves=True, internal=True, unlabeled=True, weighted=True
    ):
        hp = h - dist
        popn_to_stime[nd.get_label()] = hp
        popn_to_sdur[nd.get_label()] = nd.get_edge_length()
    popn_to_sdur[args.root_population] = 1
    popn_to_dur[args.root_population] = 1

    tables.nodes.time = [scale_time(pop_to_popn[i.population], i.time) for node, i in enumerate(tables.nodes)]

    tables.sort()
    tables.simplify(keep_unary=False)
    ntseq = tables.tree_sequence()
    ntseq.dump(args.output_path)
