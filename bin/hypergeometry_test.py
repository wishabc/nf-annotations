#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import hypergeom
from scipy.stats import fisher_exact
import argparse # using argparse for args


def calculate_pval(col):
    a = ((col == 1) & (binary_nmf_matrix['motif_indicator'] == 1)).sum()
    b = ((col == 1) & (binary_nmf_matrix['motif_indicator'] == 0)).sum()
    c = ((col == 0) & (binary_nmf_matrix['motif_indicator'] == 1)).sum()
    d = ((col == 0) & (binary_nmf_matrix['motif_indicator'] == 0)).sum()
    

    
    M = a + b + c + d   # total number of observations
    n = a + c           # all rows where motif_indicator == 1
    N = a + b           # all rows where col == 1

    # Compute the hypergeometric CDF
    p_value = hypergeom.cdf(a - 1, M, n, N)

    return 1 - p_value

def calculate_logodd(col):
    a = ((col == 1) & (binary_nmf_matrix['motif_indicator'] == 1)).sum()
    b = ((col == 1) & (binary_nmf_matrix['motif_indicator'] == 0)).sum()
    c = ((col == 0) & (binary_nmf_matrix['motif_indicator'] == 1)).sum()
    d = ((col == 0) & (binary_nmf_matrix['motif_indicator'] == 0)).sum()
    
    mat = np.array([[a, b],
                    [c, d]])
    
    # Compute odds ratio
    _, p_value = fisher_exact(mat)
    OR = (mat[0, 0] * mat[1, 1]) / (mat[0, 1] * mat[1, 0])

    # Compute log odds ratio
    log_or = np.log(OR)
    
    return log_or

# Function to aggregate p-value
def cauchy_combination(p_val):
    return stats.cauchy.sf(np.sum(np.tan(np.multiply(np.subtract(0.5, p_val), np.pi))) / len(p_val))

# Function to generate plot
def generate_top(modal, cluster, component, **kwargs):
    X_components = 'X' + str(component)
    vis = modal[modal["Components"] == X_components].head(10).sort_values(by="Estimation")
    # vis.plot(x="tf_name", y="Estimation", kind="bar")
    colors_order = get_colors_order(component)
    plt.barh(data=vis, y=cluster, width='Estimation', color = colors_order[component - 1] ,  **kwargs)
    plt.title(X_components)
    file_name = X_components + '.png'
    plt.savefig(file_name, bbox_inches='tight')
    plt.show()

def merge_aggregate_fdr(all_coeffs, motif_meta_df):
    # combine coeffs with some columns with motifs_metadata
    new_df = pd.merge(all_coeffs, motif_meta_df[['motif_id','tf_name']], on='motif_id', how='left')
    # Cauchy aggregation
    qval_df = new_df.pivot_table(index='tf_name', columns='Unnamed: 0',
                             values='Pr(>|z|)', aggfunc=cauchy_combination)
    # Do fdr corrected on p-val, fdr_by
    qval_df_fdr = qval_df.applymap(lambda x: multipletests(x, alpha=0.05,
                                                           method='fdr_by', maxiter=1)[1][0])
    # Average estimate based on tf_name
    tf_name_est_df = new_df.pivot_table(index='tf_name',
                                        columns='Unnamed: 0',values='estimate', aggfunc=np.mean)
    # Indexing based on estimation
    pvaltest_df = pd.DataFrame({'index': np.argmax(tf_name_est_df.to_numpy(), axis=1),
                                'value': tf_name_est_df.max(axis=1)})
    max_df = pvaltest_df.sort_values(by=['index', 'value'], ascending=[True, False])
    # Lock intended index on df for masking
    est = tf_name_est_df.loc[max_df.index, :]
    qvals = qval_df_fdr.loc[max_df.index, :]
    # set masking and return long format of sorted estimation by components
    filtered_est = est.where(~((est<0) | (qvals>0.05)))
    filtering = filtered_est.reset_index()
    tf_result_cluster_name = filtering.melt('tf_name',
                                            var_name='Components',
                                            value_name='Estimation').sort_values('Estimation', ascending=False)
    return tf_result_cluster_name

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='All_coeffs file to components')
    parser.add_argument('motif_id', help='value of motif id')
    parser.add_argument('indicator', help='Path to indicator file')
    args = parser.parse_args()

    nmf_matrix = pd.DataFrame(np.load("/net/seq/data2/projects/aabisheva/Encode/nextflow_results/nmf_results/bin_new_unweight_full.16.H.npy").T)

    indicator_file_df = pd.read_table(args.indicator, header=None)

    rename_dict = {i: f'comp_{i+1}' for i in range(16)}
    new_df = nmf_matrix.rename(columns=rename_dict)
    # Normalize each row by its sum
    normalized_df = new_df.div(new_df.sum(axis=1), axis=0).fillna(0)

    # Find maximum value in each row and mark as 1, others as 0
    binary_nmf_matrix = normalized_df.eq(normalized_df.max(axis=1), axis=0).astype(int)

    binary_nmf_matrix['motif_indicator'] = indicator_file_df

    results_pval = binary_nmf_matrix.drop('motif_indicator', axis=1).apply(calculate_pval)

    results_logodd = binary_nmf_matrix.drop('motif_indicator', axis=1).apply(calculate_logodd)

    pval_result = pd.DataFrame(results_pval).T
    logodd_result = pd.DataFrame(results_logodd).T

    motif_id_name = args.motif_id

    pval_result['motif_id'] = motif_id_name
    logodd_result['motif_id'] = motif_id_name

    # Write to a tab-separated file 
    name1 = motif_id_name + ".coeff.tsv"
    name2 = motif_id_name + "pval.tsv"
    logodd_result.to_csv(name1, sep='\t', index=False)
    pval_result.to_csv(name2, sep='\t', index=False)

