#!/bin/which python3

import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import multiprocessing as mp


# Params
_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
flank_width = 20
n_shuffles = 1000
es_fld = 'es_weighted_mean'

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
    df['variant_id'] = df.apply(lambda x: 
        '@'.join(map(str, 
        [x[field] for field in ['#chr', 'start', 'end', 'ref', 'alt']])), axis=1)
    df.set_index('variant_id', inplace=True)
    return df

def get_annotations(snvs, imbalanced):
    #snvs = snvs.apply(lambda x: prefered_allele(x, 'es_weighted_mean'), axis=1)
    inside = snvs['within'] == 1
    
    snvs_in = snvs[inside & imbalanced]
    if len(snvs_in) > 0:
        pred = 'log_es'

        if pred == 'log_es':
            lim = (-3, 3)
            nbins = lim[1] * 8 + 1
            fit_lim = (-2, 2)
            predictor_array = snvs_in[es_fld]
            x_0 = 0
        elif pred == 'raw_es':
            lim = (0, 1)
            nbins = 25
            predictor_array = logtr(snvs_in[es_fld])
            fit_lim = (0.2, 0.8)
            x_0 = 0.5
        elif pred == 'signed_fdr':
            lim = (-25, 25)
            nbins = 25
            predictor_array = -np.log10(snvs_in['min_fdr']) * np.sign(snvs_in['es_weighted_mean'])
            fit_lim = (-15, 15)
            x_0 = 0

        snvs_in['ddg'] = snvs_in['ref_score'] - snvs_in['alt_score']

        bins=np.linspace(*lim, nbins)

        xdiff = (bins[1]-bins[0])/2.0
        x = bins[:-1]+xdiff

        snvs_in["bin"] = pd.cut(predictor_array, bins)
        
        
        X = predictor_array
        Y = snvs_in['ddg']

        X = sm.add_constant(X)


        model = sm.OLS(Y, X).fit()

        outliers = model.get_influence().cooks_distance[0] > (4/X.shape[0])

        model = sm.OLS(Y[~outliers], X[~outliers]).fit()

        xx = np.array(fit_lim)
        yy = xx * model.params.iloc[1] + model.params.iloc[0]
        
        r2 = model.rsquared
        concordance = get_concordant(predictor_array, snvs_in["ddg"], x0=x_0)
        
        log_ref_bias = np.log2((predictor_array > 0).sum()) - np.log2((predictor_array < 0).sum())
    else:
         r2, concordance, log_ref_bias = np.nan, np.nan, np.nan

    return r2, concordance, log_ref_bias


#Enrichment code
def calc_enrichment(motifs_df, imbalanced):
    bins = np.arange(motifs_df['offset'].min(), motifs_df['offset'].max() + 2) # add 2 to length: one for '0' and one for last element
    n_all = np.histogram(motifs_df['offset'], bins=bins)[0]
    n_imbalanced = np.histogram(motifs_df['offset'][imbalanced], bins=bins)[0]
    n_not_imbalanced = n_all - n_imbalanced
    n_inside = np.nansum(n_all[flank_width:-flank_width])
    log_odds = np.log2((n_imbalanced[flank_width:-flank_width].sum() / n_imbalanced.sum()) / (n_not_imbalanced[flank_width:-flank_width].sum() / n_not_imbalanced.sum()) )
    # log_odds_per_nt = np.log2((n_imbalanced / n_imbalanced.sum()) / (n_not_imbalanced / n_not_imbalanced.sum()))

    perm = np.zeros(n_shuffles)
    # perm_per_nt = np.zeros((n_shuffles, len(bins) - 1))

    for i in range(n_shuffles):
        n_exp_imbalanced = np.histogram(motifs_df['offset'][np.random.permutation(imbalanced)], bins=bins)[0] + 1
        n_exp_not_imbalanced = n_all - n_exp_imbalanced + 1 

        perm[i] = np.log2( (n_exp_imbalanced[flank_width:-flank_width].sum() / n_exp_imbalanced.sum()) / (n_exp_not_imbalanced[flank_width:-flank_width].sum() / n_exp_not_imbalanced.sum()) )
        # perm_per_nt[i,:] = np.log2( (n_exp_imbalanced / np.sum(n_exp_imbalanced)) / (n_exp_not_imbalanced / np.sum(n_exp_not_imbalanced)))

    pval = -1 * stats.norm.logsf(log_odds,
                    loc=np.nanmean(perm, axis=0),
                    scale=np.nanstd(perm, axis=0)) / np.log(10)
    # pvals_per_nt = -1 * stats.norm.logsf(log_odds_per_nt,
    #             loc=np.nanmean(perm_per_nt, axis=0), 
    #             scale=np.nanstd(perm_per_nt, axis=0)) / np.log(10)

    return [
        log_odds,
        pval,
        n_inside,
        np.nansum(n_imbalanced[flank_width:-flank_width]),
        np.nanmedian(n_all[flank_width:-flank_width]),
        np.nansum(n_imbalanced[flank_width:-flank_width]>=7)
    ]


def get_stats(motifs_df):
    imbalanced = (motifs_df['min_fdr'] <= 0.05) # & ((motifs_df['ard'] > 0.65))# | (df['ard'] < 0.3))
    if imbalanced.sum() == 0:
        fields = [motifs_df.name, np.nan, np.nan, n_inside,
                     0, 0, 0, np.nan, np.nan, np.nan]
    else:
        log_odds, pval, n_inside, n_imb_inside, n_median_inside, n_inside_more_7  = calc_enrichment(motifs_df, imbalanced)
        r2, concordance, log_ref_bias = get_annotations(motifs_df, imbalanced)

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
            log_ref_bias
        ]

    print('\t'.join(map(str, fields)))


def main(variants_df_path, counts_df_path):
    # Load variant imbalance file
    variants_df = set_index(pd.read_table(variants_df_path))
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