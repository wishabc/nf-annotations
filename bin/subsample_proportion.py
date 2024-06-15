#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import norm
import argparse
from genome_tools.utils.sampling import stratified_sampling


def calculate_zscore(result_array, motif_counts):
    mu = result_array.mean(axis=0)
    sd = result_array.std(axis=0)
    
    z_score = (motif_counts - mu) / sd

    p_val = norm.sf(z_score)
    
    return mu, sd, z_score, p_val

def sparse_dot_product(sample_arr, binary_mat, motif_indicator):
    # sparse the matrix
    X_sparse = sparse.csr_matrix(sample_arr.T) # 1000 x dhs
    Y_sparse = sparse.csc_matrix(binary_mat) # dhs x sample

    # dot product
    result_arr = X_sparse.dot(Y_sparse) # 1000 x sample

    motifs_per_sample = np.dot(motif_indicator.T, binary_mat).squeeze() # sample

    return result_arr.toarray(), motifs_per_sample


def transform_to_bins(data, n_quantiles=100):
    bins = np.quantile(data, np.linspace(0, 1, n_quantiles + 1))
    return pd.cut(data, np.unique(bins), include_lowest=True, labels=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='All_coeffs file to components')
    parser.add_argument('motif_id', help='ID of the motif') # motif_id
    parser.add_argument('indicator', help='Path to a file with boolean indicator of overlap with significant motif hits') # indicator_file
    parser.add_argument('output', help='Path to output file') # name
    parser.add_argument('matrix_file', help='Path to matrix file') # params.nmf_matrix
    parser.add_argument('dhs_meta', help='Path to dhs annotations')
    parser.add_argument('--sample_names', help='File with sample names in the same order as the matrix', default=None)
    parser.add_argument('--n_bins', type=int, help='File with sample names in the same order as the matrix', default=100)
    args = parser.parse_args()

    indicator_file = pd.read_table(args.indicator, header=None)
    binary_matrix = np.load(args.matrix_file).astype(int)

    combined_masterlist = pd.read_table(args.dhs_meta)
    combined_masterlist['overlaps_motif'] = indicator_file
    combined_masterlist['gc_bin'] = transform_to_bins(combined_masterlist['percent_gc'], n_quantiles=args.n_bins)
    combined_masterlist['acc_bin'] = transform_to_bins(combined_masterlist['mean_acc'], n_quantiles=args.n_bins)


    sampled_indices = stratified_sampling(
        combined_masterlist,
        combined_masterlist.query('overlaps_motif == 1'),
        matching_fields=['gc_bin', 'acc_bin'],
        num_samples=1000,
        starting_seed=0,
        return_indicators=True
    )
    
    sampled_peaks, motif_hits_per_sample = sparse_dot_product(sampled_indices, binary_matrix, combined_masterlist['overlaps_motif'])
    
    mu_np, sd_np, z_score_np, pvalue = calculate_zscore(sampled_peaks, motif_hits_per_sample)

    print("Done Z-score")
    if args.sample_names is not None:
        sample_names = np.loadtxt(args.sample_names, dtype=str)
    else:
        sample_names = np.arange(binary_matrix.shape[1])

    output_df = pd.DataFrame({
        'names': sample_names,
        'mu': mu_np,
        'sd': sd_np,
        'z_score': z_score_np,
        'motif_hits_per_sample': motif_hits_per_sample,
        'p_value': pvalue
    })
 
    output_df['motif_id'] = args.motif_id

    output_df.to_csv(args.output, sep='\t', index=False)
