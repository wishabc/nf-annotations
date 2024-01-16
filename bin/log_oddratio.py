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
    
    #M = a + b + c + d
    #n = a + c
    #N = a + b
    
    p_values = np.zeros_like(log_odd_rat)
    
    # later check for case where a == expected, whether this is taken care by python
    #for i in range(len(a)):
    #    p_values[i] = 1 - hypergeom.cdf(a[i]-1, M[i], n[i], N[i])
    
    return log_odd_rat, p_values

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='All_coeffs file to components')
    parser.add_argument('motif_id', help='value of motif id') # motif_id
    parser.add_argument('indicator', help='Path to indicator file') # indicator_file
    parser.add_argument('output', help='Path to output file') # name
    parser.add_argument('matrix_file', help='Path to matrix file') # params.nmf_matrix
    parser.add_argument('meta_data', help='Path to metadata to named the sample or components') # 
    args = parser.parse_args()

    #matrix_original = np.load('/net/seq/data2/projects/aabisheva/Encode/nextflow_results/nmf_results/november_3517_samples_727k_interesting_peaks.24.H.npy').T
    #matrix_original = np.load('/net/seq/data2/projects/aabisheva/Encode/nextflow_results/nmf_results/november/november_3517_samples_727k_interesting_peaks.24.H_new.npy').T
    matrix_original = np.load(args.matrix_file) # .T
    
    indicator_matrix = pd.read_csv(args.indicator, header=None).values

    print('load data done')

    logodd, pval = log_odd_ratio_np(matrix_original, indicator_matrix)
    print('log_odd_np done')
    motif_id_name = args.motif_id
    id_name = pd.read_table(args.meta_data, usecols=['ag_id'])

    result = pd.DataFrame({'logodds': logodd, 'pvalue': pval, 'ag_id': id_name.values}) # ag_id for sample and comp for number of component
    result['motif_id'] = motif_id_name

    result.to_csv(args.output, sep='\t', index=False)
    
    # will output one file only
    #result_df.to_csv(motif_id_name + ".coeff.tsv", sep='\t', index=False)
    #np.save(motif_id_name + '.logodd.npy', np.append(logodd, motif_id_name))
    #np.save(motif_id_name + '.pval.npy', np.append(pval, motif_id_name))
