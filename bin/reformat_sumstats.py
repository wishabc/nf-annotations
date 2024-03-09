import pandas as pd
import sys
import numpy as np

result_columns = ['#chr', 'start', 'end', 'SNP', 'ref', 'alt', 'Beta', 'Beta_se', 'P', 'neglog10_p', 'INFO', 'phen_id', 'N']

def main(df, population):
    assert np.all(df['varid'] == df['chr'] + ':' + df['pos.1'].astype(str) + '_' + df['ref.1'] + '_' + df['alt.1'])
    assert np.all(df['varid'] == df['chrom'] + ':' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt'])
    try:
        df['P'] = np.power(10, -df[f'neglog10_pval_{population}'])
    except KeyError:
        return pd.DataFrame([], columns=result_columns)
    df.rename(
        columns={

        },
        inplace=True
    )

    df.rename(
        columns={
            'hg38_chr': '#chr',
            'hg38_start': 'start',
            'hg38_end': 'end',
            f'beta_{population}': 'Beta',
            f'se_{population}': 'Beta_se',
            f'neglog10_pval_{population}': 'neglog10_p',
            'info': 'INFO',
            'rsid': 'SNP'
        }, 
        inplace=True
    )
    df['#chr'] = 'chr' + df['#chr']
    return df.dropna(subset=['Beta', 'neglog10_p'])

if __name__ == '__main__':
    phen_df = pd.read_table(sys.stdin, dtype={'chr': str, 'chrom': str})
    n_cases = int(sys.argv[3])
    n_controls = None if sys.argv[4] == "N/A" else sys.argv[4]
    phen_id = sys.argv[5]
    N_eff = n_cases if n_controls is None else 4/(1/n_cases + 1/n_controls)
    main(phen_df, sys.argv[2]).assign(
        phen_id=phen_id,
        N=N_eff
    )[result_columns].to_csv(
        sys.argv[1], sep='\t',
        index=False,
    )