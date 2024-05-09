#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import norm
from genome_tools.utils.sampling import stratified_sampling
import argparse # using argparse for args
print('import packages completed')




def calculate_zscore(result_array, motif_agid):
    mu = result_array.mean(axis=0)
    sd = result_array.std(axis=0)
    
    z_score = (motif_agid - mu) / sd

    # Calculate p-value (one-tailed)
    p_val = norm.sf(np.abs(z_score))
    
    return mu, sd, z_score, p_val

def sparse_dot_product(resample_arr, binary_mat, motif_indicator):
    # sparse the matrix
    X_sparse = sparse.csr_matrix(resample_arr.T)
    Y_sparse = sparse.csc_matrix(binary_mat)

    # dot product
    result_arr = X_sparse.dot(Y_sparse)

    # motif in sample agid
    motif_agid_np = np.dot(motif_indicator.T, binary_mat).squeeze()

    return result_arr.toarray(), motif_agid_np


def transform_to_bins(data, n_quantiles=100):
    bins = np.quantile(data, np.linspace(0, 1, n_quantiles + 1))
    return pd.cut(data, np.unique(bins), include_lowest=True, labels=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='All_coeffs file to components')
    parser.add_argument('motif_id', help='ID of the motif') # motif_id
    parser.add_argument('indicator', help='Path to a file with boolean indicator of overlap with significant motif hits') # indicator_file
    parser.add_argument('output', help='Path to output file') # name
    # parser.add_argument('index_meta', help='Path to sample metadata')
    parser.add_argument('matrix_file', help='Path to matrix file') # params.nmf_matrix
    # parser.add_argument('meta_data', help='Path to metadata to named the sample or components')
    parser.add_argument('sample_meta', help='Path to sample metadata')
    parser.add_argument('dhs_meta', help='Path to dhs annotations')
    args = parser.parse_args()

    motif_id_name = args.motif_id
    indicator_file = pd.read_table(args.indicator, header=None)
    # output name
    # index_masterlist = pd.read_table(args.index_meta)
    binary_matrix = np.load(args.matrix_file) #.T
    # motifs_meta = pd.read_table(args.meta_data, header=None, names=['#chr', 'start', 'end', 'dhs_id'])
    sample_meta = pd.read_table(args.sample_meta)
    # acc_prop_df = pd.read_table(args.acc_proportion)

    combined_masterlist = pd.read_table(args.dhs_meta)
    combined_masterlist['overlaps_motif'] = indicator_file
    combined_masterlist['gc_bin'] = transform_to_bins(combined_masterlist['percent_gc'], n_quantiles=100)
    combined_masterlist['acc_bin'] = transform_to_bins(combined_masterlist['mean_acc'], n_quantiles=100)
    matching_fields = ['gc_bin', 'acc_bin']

    indices = stratified_sampling(
        combined_masterlist,
        combined_masterlist.query('overlaps_motif == 1'),
        matching_fields,
        num_samples=1000,
        starting_seed=0,
        return_indicators=True
    )
    
    result_array, motif_agid = sparse_dot_product(indices, binary_matrix, combined_masterlist['overlaps_motif'])
    
    # call function for zscore
    mu_np, sd_np, z_score_np, pvalue = calculate_zscore(result_array, motif_agid)

    print("Done Z-score")

    output_df = pd.DataFrame({
        'component_number': np.arange(binary_matrix.shape[1]),
        'mu': mu_np,
        'sd': sd_np,
        'z_score': z_score_np,
        'motif_agid': motif_agid,
        'p_value': pvalue
    })
 

    output_df['motif_id'] = motif_id_name

    output_df.to_csv(args.output, sep='\t', index=False)
