import pandas as pd
import sys
import numpy as np

result_columns = ['#chr', 'start', 'end', 'ref', 'alt', 'Beta', 'Beta_se', 'P', 'neglog10_p']

def main(df):
    try:
        df['P'] = np.power(10, -df['neglog10_pval_EUR'])
    except KeyError:
        return pd.DataFrame([], columns=result_columns)
    df['#chr'] = "chr" + df['chr']
    df['start'] = df['pos'] - 1
    df.rename(
        columns={
            'pos': 'end',
            'beta_EUR': 'Beta',
            'se_EUR': 'Beta_se',
            'neglog10_pval_EUR': 'neglog10_p'
        }, 
        inplace=True
    )
    return df[result_columns]

if __name__ == '__main__':
    phen_df = pd.read_table(sys.argv[1], compression='gzip')
    main(phen_df).to_csv(sys.argv[2], sep='\t', index=False, header=False)