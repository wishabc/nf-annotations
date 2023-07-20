#!/bin/which python3

import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm


# Params
_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
flank_width = 20
n_shuffles = 1000
es_fld = 'es_weighted_mean'

class CustomException(Exception):
    pass

def complement(x):
    return _complement[x]


def logtr(arr):
    return 1 - 1 / (1 + np.power(2, np.array(arr)))


def get_concordant(x, y, x0=0, x_mar=0):
    return ((y * (x - x0) > 0) & (np.abs(x - x0) >= x_mar)).sum() / ((np.abs(x - x0) >= x_mar) & (y != 0)).sum()


def reverse_complement(x):    
    bases = list(x) 
    bases = reversed([_complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def flip_by_strand(x):
    if x['strand'] == '-':
        x['ref'] = complement(x['ref'])
        x['alt'] = complement(x['alt'])
        #x['aa'] =  complement(x['aa'])
    return x

def get_prefered_allele(x):
    x["prefered_allele"] = x["ref"] if x[es_fld] >= 0 else x["alt"]
    return x

def set_index(df):
    if df.empty:
        exit(0)
    print(df.columns)
    df['variant_id'] = df.apply(lambda x: 
        '@'.join(map(str, 
        [x[field] for field in ['#chr', 'start', 'end', 'ref', 'alt']])), axis=1)
    df.set_index('variant_id', inplace=True)
    return df

def get_annotations(snvs, imbalanced):
    #snvs = snvs.apply(lambda x: prefered_allele(x, 'es_weighted_mean'), axis=1)
    inside = snvs['within'] == 1
    
    snvs_in = snvs[inside & imbalanced]

    try:
        if len(snvs_in) == 0:
            raise CustomException
        
        pred = 'log_es'

        if pred == 'log_es':
            predictor_array = snvs_in[es_fld]
            x_0 = 0
        elif pred == 'raw_es':
            predictor_array = logtr(snvs_in[es_fld])
            x_0 = 0.5
        elif pred == 'signed_fdr':
            predictor_array = -np.log10(snvs_in['min_fdr']) * np.sign(snvs_in['es_weighted_mean'])
            x_0 = 0

        if (predictor_array != 0).sum() == 0:
            raise CustomException

        snvs_in['ddg'] = snvs_in['ref_score'] - snvs_in['alt_score']
        
        X = predictor_array
        Y = snvs_in['ddg']

        X = sm.add_constant(X)

        model = sm.OLS(Y, X).fit()

        outliers = model.get_influence().cooks_distance[0] > (4/X.shape[0])

        model = sm.OLS(Y[~outliers], X[~outliers]).fit()
        
        r2 = model.rsquared
        concordance = get_concordant(predictor_array, snvs_in["ddg"], x0=x_0)
        
        ref_bias = (predictor_array > 0).sum() / (predictor_array != 0).sum()
    except CustomException:
        r2, concordance, ref_bias = np.nan, np.nan, np.nan

    return r2, concordance, ref_bias


#Enrichment code
def calc_enrichment(motifs_df, imbalanced):
    bins = np.arange(motifs_df['offset'].min(), motifs_df['offset'].max() + 2) # add 2 to length: one for '0' and one for last element
    n_all = np.histogram(motifs_df['offset'], bins=bins)[0]
    n_imbalanced = np.histogram(motifs_df['offset'][imbalanced], bins=bins)[0]
    n_not_imbalanced = n_all - n_imbalanced

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
        np.nansum(n_imbalanced[flank_width:-flank_width]>=7)
    ]


def get_stats(motifs_df):
    imbalanced = (motifs_df['min_fdr'] <= 0.05) # & ((motifs_df['ard'] > 0.65))# | (df['ard'] < 0.3))
    log_odds, pval, n_inside, n_imb_inside, n_median_inside, n_inside_more_7  = calc_enrichment(motifs_df, imbalanced)
    if imbalanced.sum() == 0:
        r2, concordance, ref_bias = np.nan, np.nan, np.nan
    else:
        r2, concordance, ref_bias = get_annotations(motifs_df, imbalanced)

    fields = [
        motifs_df.name,
        log_odds,
        pval,
        n_inside,
        n_imb_inside, 
        n_median_inside,
        n_inside_more_7,
        r2,
        concordance,
        ref_bias
    ]

    print('\t'.join(map(str, fields)))


def main(variants_df_path, counts_df_path):
    # Load variant imbalance file
    names = "#chr    start   end     ID      ref     alt     ref_counts      alt_counts      sample_id       AAF     RAF     FMR #chr.1  start.1 end.1   BAD     SNP_count       SNP_ID_count    sum_cover       Q1.00   Q1.50   Q2.00   Q2.50   Q3.00   Q4.00   Q5.00   Q6.00   coverage        w       es      pval_ref        pval_alt        is_tested       footprints      hotspots        group_id".split()

    variants_df = set_index(
        pd.read_table(variants_df_path, header=None, names=[x for x in names if x != ' '])
    )
    if variants_df.empty:
        return
    # Load motifs dataframe
    motifs_df = set_index(pd.read_table(counts_df_path,
        header=None,
        names=['#chr', 'start', 'end', 'rsid', 'ref', 'alt',
         'motif', 'offset', 'within', 'strand', 'ref_score', 'alt_score', 'seq']))

    # df_ct_all = pd.crosstab(snvs['offset'], snvs['ref'])
    # df_ct_imbalanced = pd.crosstab(snvs.loc[imbalanced]['offset'],
    #         snvs.loc[imbalanced]['prefered_allele'])
    # Flip data with motif is on '-' strand
    motifs_df = motifs_df.apply(flip_by_strand, axis=1)

    # Add imbalance data
    df = variants_df[[es_fld, 'min_fdr']].join(motifs_df)

    # Compute preferred allele
    df = df.apply(get_prefered_allele, axis=1)
    df.groupby('motif').apply(get_stats)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])