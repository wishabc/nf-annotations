#!/bin/which python3

import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
from tqdm import tqdm
import dask.dataframe as dd
tqdm.pandas()

# Params
_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
flank_width = 20
n_shuffles = 1000
es_fld = 'es_weighted_mean'


columns = [
    'motifs_name',
    "log_odds",
    "pval",
    "total_inside",
    "imbalanced_inside",
    "imbalanced_inside_median",
    "n_imbalanced_more_7",
    'r2',
    'concordance',
    'ref_bias'
]

class NoDataException(Exception):
    pass


def complement(x):
    return np.vectorize(_complement.get)(x)


def logtr(arr):
    return 1 - 1 / (1 + np.power(2, np.array(arr)))


def get_concordant(x, y, expected_es=0, x_mar=0):
    diff = x - expected_es
    unmasked_values = (np.abs(diff) >= x_mar) & (y != 0)
    return ((y * diff > 0) & unmasked_values).sum() / unmasked_values.sum()


def join_columns(df):
    df['variant_id'] = df[['#chr', 'start', 'end', 'ref', 'alt']].astype(str).agg('@'.join, axis=1)
    return df



def set_index(df):
    if len(df.index) == 0:
        exit(0)
    df = df.map_partitions(join_columns)
    return df.set_index('variant_id')


def get_annotations(snvs):

    snvs_in = snvs[(snvs['within'] == 1)]

    try:
        if len(snvs_in) == 0:
            raise NoDataException()

        predictor_array = snvs_in[es_fld].to_numpy()
        expected_es = 0

        if (predictor_array == 0).all() == 0:
            raise NoDataException()
        
        X = predictor_array
        Y = snvs_in['ddg']

        X = sm.add_constant(X)

        model = sm.OLS(Y, X).fit()

        outliers = model.get_influence().cooks_distance[0] > (4/X.shape[0])

        model = sm.OLS(Y[~outliers], X[~outliers]).fit()
        
        r2 = model.rsquared

        concordance = get_concordant(
            predictor_array, 
            snvs_in["ddg"], 
            expected_es=expected_es
        )
        
        ref_bias = (predictor_array > 0).sum() / (predictor_array != 0).sum()
    except NoDataException:
        r2, concordance, ref_bias = np.nan, np.nan, np.nan

    return r2, concordance, ref_bias


#Enrichment code
def calc_enrichment(motifs_df, imbalanced):
    bins = np.arange(motifs_df['offset'].min(), motifs_df['offset'].max() + 2) # add 2 to length: one for '0' and one for last element
    n_all = np.histogram(motifs_df['offset'], bins=bins)[0]
    n_imbalanced = np.histogram(motifs_df['offset'][imbalanced], bins=bins)[0]

    total_inside = np.nansum(n_all[flank_width:-flank_width])
    total_imbalanced_inside = np.nansum(n_imbalanced[flank_width:-flank_width])
    if total_imbalanced_inside == 0 or total_inside - total_imbalanced_inside == 0:
        return np.nan, np.nan, total_inside, 0, 0, 0

    log_odds = np.log2(total_imbalanced_inside) - np.log2(total_inside - total_imbalanced_inside) - \
        np.log2(np.nansum(n_imbalanced) - total_imbalanced_inside) + \
            np.log2(np.nansum(n_all) - np.nansum(n_imbalanced) - total_inside + total_imbalanced_inside)
    pval = -stats.hypergeom.logsf(
        total_imbalanced_inside,
        np.nansum(n_all),
        np.nansum(n_imbalanced),
        total_inside
    ) / np.log(10)

    return [
        log_odds,
        pval,
        total_inside,
        np.nansum(n_imbalanced[flank_width:-flank_width]),
        np.nanmedian(n_all[flank_width:-flank_width]),
        np.nansum(n_imbalanced[flank_width:-flank_width] >= 7)
    ]


def get_stats(motifs_df):
    imbalanced_index = motifs_df['min_fdr'] <= 0.05
    try:
        if imbalanced_index.sum() == 0:
            raise NoDataException()

        r2, concordance, ref_bias = get_annotations(motifs_df[imbalanced_index])

        # log_odds, pval, n_inside, n_imb_inside, n_median_inside, n_inside_more_7
        enrichment_stats = calc_enrichment(motifs_df, imbalanced_index)
    except NoDataException:
        return pd.Series([], columns=columns)

    data = pd.Series([
        motifs_df.name,
        *enrichment_stats,
        r2,
        concordance,
        ref_bias
    ], columns=columns)

    return data


def main(variants_df_path, counts_df_path):
    # Load variant imbalance file
    names = "#chr    start   end     ID      ref     alt     ref_counts      alt_counts      sample_id       AAF     RAF     FMR #chr.1  start.1 end.1   BAD     SNP_count       SNP_ID_count    sum_cover       Q1.00   Q1.50   Q2.00   Q2.50   Q3.00   Q4.00   Q5.00   Q6.00   coverage        w       es      pval_ref        pval_alt        is_tested       footprints      hotspots        group_id".split()

    print('Reading variants df')
    variants_df = set_index(
        dd.read_table(variants_df_path,
            header=None,
            dtype={'AAF': 'object', 'RAF': 'object'},
            names=[x for x in names if x != ' '])
    )
    if variants_df.empty:
        return
    # Load motifs dataframe
    print('Reading motifs df')
    motifs_df = set_index(
        dd.read_table(
            counts_df_path,
            header=None,
            names=['#chr', 'start', 'end', 'rsid', 'ref', 'alt',
            'motif', 'offset', 'within', 'strand',
            'ref_score', 'alt_score', 'seq']
        )
    )
    for key in ('ref', 'alt'):
        motifs_df[key] = np.where(motifs_df['strand'] == '-', complement(motifs_df[key]), motifs_df[key])
   
    # Add imbalance data
    df = variants_df[[es_fld, 'min_fdr']].join(motifs_df)

    # Compute preferred allele
    df["prefered_allele"] = np.where(df[es_fld] >= 0, df["ref"], df["alt"])
    df['ddg'] = df['ref_score'] - df['alt_score']
    return df.groupby('motif').apply(get_stats).compute()


if __name__ == '__main__':
    res_df = main(sys.argv[1], sys.argv[2])
    res_df.to_csv(sys.argv[3], index=False, sep='\t')