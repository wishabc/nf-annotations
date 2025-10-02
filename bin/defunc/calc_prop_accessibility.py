import pandas as pd
import numpy as np
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Calculate proportion of accessibility')
    parser.add_argument('binary_matrix', help='Path to binary matrix to calculate proportion of accessibility')
    parser.add_argument('dhs_meta', help='Path to DHS metadata (index file)')
    parser.add_argument('output', help='Path to the output file')
    parser.add_argument('--samples_weights', help='Path to samples weights (to avoid class imbalanced)', default=None)
    args = parser.parse_args()

    dhs_annotations = pd.read_table(args.dhs_meta)
    binary_matrix = np.load(args.binary_matrix).astype(np.float16)
    assert binary_matrix.shape[0] == len(dhs_annotations)
    assert binary_matrix.shape[0] == dhs_annotations.shape[0]

    if args.samples_weights is not None:
        sample_weights = pd.read_table(args.samples_weights)['weight'].to_numpy()
    else:
        sample_weights = np.ones(binary_matrix.shape[1])

    
    acc_proportions = np.average(binary_matrix, axis=1, weights=sample_weights)
    dhs_annotations['mean_acc'] = acc_proportions
    dhs_annotations.to_csv(args.output, index=False, sep="\t")