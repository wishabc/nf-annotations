import pandas as pd
import sys
import numpy as np


def main(df):
    df['P'] = np.power(10, -df['neglog10_pval_meta_hq'])
    df['#chr'] = "chr" + df['chr']
    df['start'] = df['pos'] - 1
    df.rename(
        columns={
            'pos': 'end',
            'beta_meta_hq': 'Beta',
            'se_meta_hq': 'Beta_se',
            'neglog10_pval_meta_hq': 'neglog10_p'
        }, inplace=True
    )
    return df[['#chr', 'start', 'end', 'ref', 'alt', 'Beta', 'Beta_se', 'P', 'neglog10_p']]

if __name__ == '__main__':
    phen_df = pd.read_table(sys.argv[1], compression='gzip')
    main(phen_df).to_csv(sys.argv[2], sep='\t', index=False, header=False)