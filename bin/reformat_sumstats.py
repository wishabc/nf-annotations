import pandas as pd
import sys
import numpy as np

result_columns = ['#chr', 'start', 'end', 'ref', 'alt', 'Beta', 'Beta_se', 'P', 'neglog10_p']

def main(df, population):
    try:
        df['P'] = np.power(10, -df[f'neglog10_pval_{population}'])
    except KeyError:
        return pd.DataFrame([], columns=result_columns)
    df['#chr'] = "chr" + df['chr'].astype(str)
    df['start'] = df['pos'] - 1
    df.rename(
        columns={
            'pos': 'end',
            f'beta_{population}': 'Beta',
            f'se_{population}': 'Beta_se',
            f'neglog10_pval_{population}': 'neglog10_p'
        }, 
        inplace=True
    )
    return df[result_columns]

if __name__ == '__main__':
    phen_df = pd.read_table(sys.argv[1], compression='gzip')
    main(phen_df, sys.argv[3]).to_csv(sys.argv[2], sep='\t', index=False, header=False)