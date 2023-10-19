#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import hypergeom
import argparse # using argparse for args
print('past the import stuff')
def log_odd_ratio(ori_df, indicator_df):
    
    log_odds_ratios = []
    pval_fischer = []
    
    for col in ori_df.columns:
        # motif=1, sample=1
        a = (ori_df[col] & indicator_df.iloc[:, 0]).sum()
        # motif=1, sample=0
        b = indicator_df.iloc[:, 0].sum() - a
        # motif=0, sample=1
        c = ori_df[col].sum() - a
        # motif=0, sample=0
        d = (len(ori_df[col]) - ori_df[col].sum()) - b
        
        # Calculate the log odds ratio
        log_odds_ratios.append(np.log((a * d) / (b * c)))
        
        # Compute p-value
        try:
            _, p_value = fisher_exact(np.array([[a, b], [c, d]]))
        except ValueError:
            print(a, b, c, d)
            raise

        pval_fischer.append(p_value)
        

    return pd.DataFrame({'log_odd_ratio' : log_odds_ratios,
                        'fischer_pval' : pval_fischer},
                        index=ori_df.columns)

def log_odd_ratio_np(ori_matrix_np, indicator_df):
    print('inside function')
    a = (ori_matrix_np * indicator_df).sum(axis=0) # (3883,)
    b = indicator_df.sum() - a # (3883,)
    c = ori_matrix_np.sum(axis=0) - a # (3883,)
    d = len(ori_matrix_np) - (a + b + c)
    
    #log_odd_rat = np.log((a * d) / (b * c))
    log_odd_rat = np.where((b * c) != 0, np.log((a * d) / (b * c)), 0)
    
    M = a + b + c + d
    n = a + c
    N = a + b
    
    p_values = np.zeros_like(log_odd_rat)
    
    # later check for case where a == expected, whether this is taken care by python
    for i in range(len(a)):
        p_values[i] = 1 - hypergeom.cdf(a[i]-1, M[i], n[i], N[i])
    
    return log_odd_rat, p_values

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='All_coeffs file to components')
    parser.add_argument('motif_id', help='value of motif id')
    parser.add_argument('indicator', help='Path to indicator file')
    args = parser.parse_args()

    # matrix_original_df = pd.DataFrame(np.load('/net/seq/data2/projects/afathul/motif_enrichment/odd_ratio/nmf_32.H.npy').T)
    # matrix_original_df.columns = ['comp_' + str(i) for i in list(range(1, len(matrix_original_df.columns) + 1))]
    # Find maximum value in each row and mark as 1, others as 0
    # binary_nmf_matrix = matrix_original_df.eq(matrix_original_df.max(axis=1), axis=0).astype(int)

    matrix_original = np.load('/net/seq/data2/projects/afathul/motif_enrichment/odd_ratio/original_binary_matrix.npy').T
    # sample_metadata = pd.read_table('/net/seq/data2/projects/afathul/motif_enrichment/odd_ratio/sample_metadata_new.tsv')
    # matrix_original_df.columns = sample_metadata['ag_id'].values
    indicator_matrix = pd.read_csv(args.indicator, header=None).values

    print('load data done')
    # result_df = log_odd_ratio(binary_nmf_matrix, indicator_df).reset_index()
    logodd, pval = log_odd_ratio_np(matrix_original, indicator_matrix)
    print('log_odd_np done')
    motif_id_name = args.motif_id
    # result_df['motif_id'] = motif_id_name

    # will output one file only
    #result_df.to_csv(motif_id_name + ".coeff.tsv", sep='\t', index=False)
    np.save(motif_id_name + '.logodd.npy', np.append(logodd, motif_id_name))
    np.save(motif_id_name + '.pval.npy', np.append(pval, motif_id_name))
