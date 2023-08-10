import pandas as pd
import sys
import numpy as np

result_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'Beta', 'Beta_se', 'P', 'neglog10_p', 'INFO']

def main(df, population):
    try:
        df['P'] = np.power(10, -df[f'neglog10_pval_{population}'])
    except KeyError:
        return pd.DataFrame([], columns=result_columns)
    df['#chr'] = "chr" + df['chr']
    df['start'] = df['pos'] - 1
    assert df[['chrom', 'pos', 'ref', 'alt']].equals(df[['chr', 'pos.1', 'ref.1', 'alt.1']])
    df.rename(
        columns={
            'pos': 'end',
            f'beta_{population}': 'Beta',
            f'se_{population}': 'Beta_se',
            f'neglog10_pval_{population}': 'neglog10_p',
            'info': 'INFO',
            'rsid': 'ID'
        }, 
        inplace=True
    )
    return df[result_columns].dropna(subset=['Beta', 'neglog10_p'])

if __name__ == '__main__':
    phen_df = pd.read_table(sys.stdin, dtype={'chr': str, 'chrom': str})
    print(phen_df.columns)
    main(phen_df, sys.argv[2]).to_csv(sys.argv[1], sep='\t', index=False, header=False)