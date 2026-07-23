import sys

if __name__ == "__main__":
    input_filepath = sys.argv[1]
    output_filepath = sys.argv[2]
    edge_label = sys.argv[3]
    change_coef = float(sys.argv[4])

    label_to_size = {}
    with open(input_filepath, 'r') as f:
        for ix, line in enumerate(f):
            if ix == 0:
                continue
            label, size = line.strip().split()
            size = float(size)
            if label == edge_label:
                label_to_size[label] = size * change_coef
            else:
                label_to_size[label] = size

    with open(output_filepath, 'w') as f:
        for k, v in label_to_size.items():
            f.write(f"{k}\t{v}\n")
