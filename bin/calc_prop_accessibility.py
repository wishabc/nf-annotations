import pandas as pd
import numpy as np
import argparse
from scipy.sparse import csc_matrix
from tqdm import tqdm

def sample_with_weights(weights, n=1000, n_samples=10000, seed=None):
    prob_distribution = weights / weights.sum()
    
    # Pre-allocate a matrix to store boolean masks
    masks = np.zeros((n_samples, len(weights)), dtype=bool)
    
    # Generate 10,000 samples
    dat = np.arange(len(weights))
    for i in range(n_samples):
        if seed is not None:
            seed += n_samples
            np.random.seed(seed)
        # Draw n unique sample indices using the probability distribution
        sampled_indices = np.random.choice(
            dat, size=n,
            replace=False,
            p=prob_distribution
        )
        
        # Create a boolean mask for the sampled indices
        masks[i, sampled_indices] = True
    
    return masks


def main(binary_matrix, sample_weights, n=1000, n_samples=1000):

    sampled_masks = sample_with_weights(sample_weights, n=n, n_samples=n_samples, seed=0)

    acc_counts = np.zeros((n_samples, binary_matrix.shape[0]), dtype=int)

    # Apply each mask and calculate the sums directly in a sparse-efficient way
    for i, mask in tqdm(enumerate(sampled_masks), total=len(sampled_masks)):
        acc_counts[i, :] = binary_matrix[:, mask].mean(axis=1)

    return acc_counts



if __name__ == '__main__':
    parser = argparse.ArgumentParser('Calculate proportion of accessibility')
    parser.add_argument('binary_matrix', help='Path to binary matrix to calculate proportion of accessibility')
    parser.add_argument('dhs_meta', help='Path to DHS metadata (index file)')
    parser.add_argument('dhs_annotations', help='DHS annotations file (BED format)')
    parser.add_argument('ouput', help='Path to the output file')
    parser.add_argument('--samples_weights', help='Path to samples weights (to avoid class imbalanced)', default=None)
    parser.add_argument('--sampling_n', type=int, help='Number of sampling iterations to average proportion of accessibility', default=1000)
    parser.add_argument('--n', type=int, help='Number of samples to choose', default=None)
    args = parser.parse_args()

    dhs_ids = pd.read_table(args.dhs_meta, header=None)[3]
    dhs_annotations = pd.read_table(args.dhs_annotations)
    dhs_annotations = dhs_annotations[dhs_annotations['dhs_id'].isin(dhs_ids)].reset_index(drop=True)
    binary_matrix = np.load(args.binary_matrix).astype(int)
    
    print(binary_matrix.shape)

    if args.samples_weights is not None:
        sample_weights = pd.read_table(args.samples_weights)['weight'].to_numpy()
    else:
        sample_weights = np.ones(binary_matrix.shape[1])
    
    if args.n is None:
        n = binary_matrix.shape[1]
    else:
        n = args.n

    assert binary_matrix.shape[0] == dhs_annotations.shape[0]
    acc_proportions = main(binary_matrix, sample_weights, n=n, n_samples=args.sampling_n)
    dhs_annotations['mean_acc'] = acc_proportions.mean(axis=0)
    dhs_annotations['std_acc'] = acc_proportions.std(axis=0)
    dhs_annotations.to_csv(args.output, index=False, sep="\t")