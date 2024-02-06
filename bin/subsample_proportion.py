#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy import sparse
import argparse # using argparse for args
print('import packages completed')

def shuffle_along_axis(a, axis):
    idx = np.random.rand(*a.shape).argsort(axis=axis)
    return np.take_along_axis(a,idx,axis=axis)

def apply_funct(df_for_gc_bin):
    matrix = shuffle_along_axis(np.tile(df_for_gc_bin['motif_indicator_file'].values, (1000, 1)).T, axis=0)
    return pd.DataFrame(matrix, index=df_for_gc_bin.index)

def bin_resample(df_masterlist):

    # Obtained s_cut and bins which is the edge for our DHSs with motif
    s_cut, bins = pd.cut(df_masterlist[df_masterlist['motif_indicator_file'] == 1]['percent_gc'],
                         bins=50, retbins=True)

    # Add two empty bins such that 0 to min and max to 1
    new_bin_edges = sorted([0] + list(bins) + [1])
    df_masterlist['bins'] = pd.cut(df_masterlist['percent_gc'], new_bin_edges, include_lowest=True)

    # Group by bins and apply function:
    # shuffle_along_axis, apply_funct: Shuffle the sample based on bins given and sample 1000 times
    result_df = df_masterlist.groupby('bins').apply(apply_funct)

    return result_df

def calculate_zscore(result_array, motif_agid):
    mu = result_array.mean(axis=0)
    sd = result_array.std(axis=0)
    
    z_score = (motif_agid - mu) / sd
    
    return mu, sd, z_score

def sparse_dot_product(resample_arr, binary_mat, indicator_df):
    # sparse the matrix
    X_sparse = sparse.csr_matrix(resample_arr.T)
    Y_sparse = sparse.csc_matrix(binary_mat)

    # dot product
    result_arr = X_sparse.dot(Y_sparse)

    # motif in sample agid
    motif_agid_np = np.dot(indicator_df.to_numpy().T, binary_mat)[0]

    return result_arr, motif_agid_np


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='All_coeffs file to components')
    parser.add_argument('motif_id', help='value of motif id') # motif_id
    parser.add_argument('indicator', help='Path to indicator file') # indicator_file
    parser.add_argument('output', help='Path to output file') # name
    parser.add_argument('matrix_file', help='Path to matrix file') # params.nmf_matrix
    parser.add_argument('meta_data', help='Path to metadata to named the sample or components')
    parser.add_argument('sample_meta', help='Path to sample metadata')
    args = parser.parse_args()

    motif_id_name = args.motif_id
    index_masterlist = pd.read_table(args.index_m)
    motifs_meta = pd.read_table(args.meta_data, header=None, names=['#chr', 'start', 'end', 'dhs_id'])
    indicator_file = pd.read_table(args.indicator, header=None)
    binary_matrix = np.load(args.matrix_file)
    sample_meta = pd.read_table(args.sample_meta)

    # Merge so that it filter out only Index DHSs
    combined_masterlist = motifs_meta.merge(index_masterlist)

    # Save the total number of motif
    total_DHSmotif = sum(indicator_file[0])

    # Combined the motif indicator
    combined_masterlist['motif_indicator_file'] = indicator_file

    # run bins and resample by bins function
    resample_array = bin_resample(combined_masterlist)


    result_array, motif_agid = sparse_dot_product(resample_array, binary_matrix, indicator_file)
    
    # call function for zscore
    mu_np, sd_np, z_score_np = calculate_zscore(result_array.toarray(), motif_agid)

    d = {'ag_id': sample_meta['ag_id'].values, 'mu': mu_np, 'sd': sd_np,
         'z_score': z_score_np, 'motif_agid': motif_agid}
    output_df = pd.DataFrame(data = d)

    output_df['motif_id'] = motif_id_name

    output_df.to_csv(args.output, sep='\t', index=False)
