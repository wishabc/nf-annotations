#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import norm
from numba import jit, prange
from numpy.lib.stride_tricks import as_strided
from genome_tools.utils.sampling import stratified_sampling
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
    motif_agid_np = np.dot(motif_indicator.T, binary_mat)[0]

    return result_arr, motif_agid_np

def get_probs(bin_distances, w):
    transformed = np.exp(-bin_distances / (2*w**2))
    return transformed / transformed.sum(axis=1)

@jit(nopython=True, parallel=True)
def sample_multinomial(n_array, p_matrix, num_samples):
    """Sample for each n and row of ps with parallel execution."""
    num_categories = p_matrix.shape[1]
    results = np.zeros((num_samples, num_categories), dtype=np.int32)

    for i in prange(num_samples):
        for j in range(len(n_array)):
            n = n_array[j]
            p = p_matrix[j]
            counts = np.zeros(len(p), dtype=np.int32)

            # Inlined categorical sampling
            for _ in range(n):
                r = np.random.random()
                cumulative = 0.0
                for k, p_val in enumerate(p):
                    cumulative += p_val
                    if r < cumulative:
                        counts[k] += 1
                        break

            results[i] += counts
    return results

def perturb_bin_counts(bin_counts, w=0.01, num_samples=1000):
    bin_distances = squareform(pdist(np.array(bin_counts.index.tolist()), 'euclidean'))
    if w > 0:
        probs = get_probs(bin_distances, w)
        return sample_multinomial(bin_counts.values, probs, num_samples)
    else:
        return np.tile(bin_counts.values, (num_samples, 1))

# Changed to parallel=False
@jit(nopython=True, parallel=False)
def get_sample_indicators(counts_to_sample, all_counts, seed=0):
    np.random.seed(seed)
    total_rows = all_counts.sum()
    num_samples, num_columns = counts_to_sample.shape
    print((total_rows, num_samples))
    res = np.zeros((total_rows, num_samples), dtype=np.bool_)

    cums = 0
    for i in range(num_columns):
        all_c = all_counts[i]
        for j in range(num_samples):
            c_to_sample = counts_to_sample[j, i]
            binary_matrix = np.arange(all_c) < c_to_sample
            shuffled_indices = np.random.permutation(all_c)
            res[cums: cums + all_c, j] = binary_matrix[shuffled_indices]
            seed += 10000
            np.random.seed(seed)
        cums += all_c

    return res


def transform_to_bins(data, n_quantiles=100):
    bins = np.quantile(data, np.linspace(0, 1, n_quantiles + 1))
    return pd.cut(data, np.unique(bins), include_lowest=True, labels=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='All_coeffs file to components')
    parser.add_argument('motif_id', help='value of motif id') # motif_id
    parser.add_argument('matrix_type', help='Path to type of matrix file') # New addition
    parser.add_argument('indicator', help='Path to indicator file') # indicator_file
    parser.add_argument('output', help='Path to output file') # name
    # parser.add_argument('index_meta', help='Path to sample metadata')
    parser.add_argument('matrix_file', help='Path to matrix file') # params.nmf_matrix
    # parser.add_argument('meta_data', help='Path to metadata to named the sample or components')
    parser.add_argument('sample_meta', help='Path to sample metadata')
    parser.add_argument('dhs_meta', help='Path to dhs annotations')
    args = parser.parse_args()

    motif_id_name = args.motif_id
    matrix_type = args.matrix_type
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
    combined_masterlist['acc_bin'] = transform_to_bins(combined_masterlist['acc_prop'], n_quantiles=100)
    matching_fields = ['gc_bin', 'acc_bin']

    indices = stratified_sampling(
        combined_masterlist,
        combined_masterlist.query('overlaps_motif == 1'),
        matching_fields,
        num_samples=1000,
        starting_seed=0,
    )
    
    result_array, motif_agid = sparse_dot_product(indices, binary_matrix, indicator_file)

    print("Done Dot Product")
    
    # call function for zscore
    mu_np, sd_np, z_score_np, pvalue = calculate_zscore(result_array.toarray(), motif_agid)

    print("Done Z-score")

    if matrix_type == "NMF":
        d = {'component_number': [i for i in range(1,binary_matrix.shape[1] + 1)], 'mu': mu_np, 'sd': sd_np,
         'z_score': z_score_np, 'motif_agid': motif_agid, 'p_value': pvalue}
        output_df = pd.DataFrame(data = d)
    elif matrix_type == "DHS_Binary":
        d = {'ag_id': sample_meta['ag_id'].values, 'mu': mu_np, 'sd': sd_np,
         'z_score': z_score_np, 'motif_agid': motif_agid, 'p_value': pvalue}
        output_df = pd.DataFrame(data = d)
    else:
        print("Type of matrix not specified in workflow")    

    output_df['motif_id'] = motif_id_name
    output_df['matrix_type'] = matrix_type

    output_df.to_csv(args.output, sep='\t', index=False)
