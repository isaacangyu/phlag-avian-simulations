import os
import sys
import argparse
from pathlib import Path
import treeswift as ts


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

def generate_mapping(tree_path: Path, node_label: str, output_path: Path):
    tree = ts.read_tree_newick(tree_path)
    label_to_node = tree.label_to_node(selection="all")
    
    if node_label not in label_to_node:
        sys.exit(f"Error: Internal node '{node_label}' not found in species tree '{tree_path}'.")
        
    nd_target = label_to_node[node_label]
    
    # P2 = Leaves under target internal node
    p2_leaves = set(nd.get_label() for nd in tree.extract_subtree(nd_target).traverse_leaves())
    
    # P1 = Leaves under sister clade(s) of target internal node
    parent_node = nd_target.get_parent()
    p1_leaves = set()
    if parent_node:
        for child in parent_node.child_nodes():
            if child != nd_target:
                p1_leaves.update(nd.get_label() for nd in tree.extract_subtree(child).traverse_leaves())
                
    # P3 = Leaves under sister clade(s) of parent node
    p3_leaves = set()
    grandparent_node = parent_node.get_parent() if parent_node else None
    if grandparent_node:
        for child in grandparent_node.child_nodes():
            if child != parent_node:
                p3_leaves.update(nd.get_label() for nd in tree.extract_subtree(child).traverse_leaves())
                
    # Po = All remaining leaves in species tree
    all_leaves = set(nd.get_label() for nd in tree.traverse_leaves())
    po_leaves = all_leaves - p1_leaves - p2_leaves - p3_leaves
    
    with open(output_path, "w") as f:
        for leaf in sorted(p1_leaves):
            f.write(f"{leaf}\tP1\n")
        for leaf in sorted(p2_leaves):
            f.write(f"{leaf}\tP2\n")
        for leaf in sorted(p3_leaves):
            f.write(f"{leaf}\tP3\n")
        for leaf in sorted(po_leaves):
            f.write(f"{leaf}\tPo\n")
            
    print(f"Successfully generated mapping at '{output_path}':")
    print(f"  P1: {len(p1_leaves)} leaves")
    print(f"  P2: {len(p2_leaves)} leaves")
    print(f"  P3: {len(p3_leaves)} leaves")
    print(f"  Po: {len(po_leaves)} leaves")

def main():
    parser = argparse.ArgumentParser(
        description="Generate 4-taxon population mapping TSV for dstar based on an internal node in a species tree."
    )
    parser.add_argument(
        "-t", "--tree", "--input-tree",
        dest="input_tree",
        type=Path,
        required=True,
        help="Path to species tree file in Newick format."
    )
    parser.add_argument(
        "-n", "--node", "--internal-node",
        dest="internal_node",
        type=str,
        required=True,
        help="Label of internal node (e.g. N151)."
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        type=Path,
        default=None,
        help="Output mapping TSV path (default: <tree_stem>_<node>_mapping.tsv)."
    )
    
    args = parser.parse_args()
    
    if args.output is None:
        output_dir = get_output_dir()
        (output_dir / "mapping").mkdir(parents=True, exist_ok=True)
        args.output = output_dir / "mapping" / f"{args.input_tree.stem}_{args.internal_node}_mapping.tsv"
    elif not args.output.is_absolute():
        output_dir = get_output_dir()
        (output_dir / "mapping").mkdir(parents=True, exist_ok=True)
        args.output = output_dir / "mapping" / args.output.name
        
    generate_mapping(args.input_tree, args.internal_node, args.output)

if __name__ == "__main__":
    main()
