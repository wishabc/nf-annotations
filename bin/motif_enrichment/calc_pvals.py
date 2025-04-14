import scipy.stats as st
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import sys


def sampling_to_signif(df):
    df['frac'] = df.eval('motif_hits / total')
    ref = df.query('sampling == "reference"')
    sampled = df.query('sampling != "reference"')
    ref_set = ref.set_index(['prefix', 'motif_id']).sort_values('frac', ascending=False)
    stats = sampled.groupby(['prefix', 'motif_id']).agg(
        sampled_frac=('frac', 'mean'),
        sampled_frac_var=('frac', 'var')
    ).join(
        ref_set
    ).sort_values(
        'frac',
        ascending=False
    )
    stats['z'] = stats.eval('(frac - sampled_frac) / sqrt(sampled_frac_var)')
    stats['logodds'] = stats.eval('log(sampled_frac) - log(1 - sampled_frac) - log(frac) + log(1 - frac)') / np.log(2)
    stats['pval'] = st.norm.sf(stats['z'])
    stats['fdr'] = multipletests(stats['pval'], method='fdr_bh')[1]
    return stats.reset_index().drop(columns=['sampling'])


if __name__ == "__main__":
    sampling_stats = pd.read_table(sys.argv[1])
    enrichments = sampling_to_signif(sampling_stats)
    enrichments.to_csv(
        sys.argv[2],
        sep='\t',
        index=False,
    )