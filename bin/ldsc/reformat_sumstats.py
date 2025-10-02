import pandas as pd
import sys
import numpy as np

result_columns = ['#chr', 'start', 'end', 'SNP', 'ref', 'alt', 'Beta', 'Beta_se', 'P', 'neglog10_p', 'INFO', 'phen_id', 'N']

def get_req_columns(population):
    return [
        'hg38_chr', 'hg38_start', 'hg38_end', 'chr', 'chrom', 'pos', 'pos.1', 'ref', 'ref.1', 'alt', 'alt.1',
        f'beta_{population}', f'se_{population}', f'neglog10_pval_{population}', 
        'info', 'rsid'
    ]

def main(df, population):  
    try:
        df['P'] = np.power(10, -df[f'neglog10_pval_{population}'])
    except KeyError:
        return pd.DataFrame([], columns=result_columns)

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
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    return df

if __name__ == '__main__':
    population = sys.argv[2]
    columns = get_req_columns(population)
    phen_df = pd.read_table(sys.stdin, dtype={'chr': str, 'chrom': str}, usecols=columns)
    assert np.all(
        (phen_df['chrom'] == phen_df['chr']) & 
        (phen_df['pos'] == phen_df['pos.1']) & 
        (phen_df['ref'] == phen_df['ref.1']) & 
        (phen_df['alt'] == phen_df['alt.1'])
    )
    phen_df.dropna(inplace=True)
    n_samples = pd.to_numeric(sys.argv[3], errors='coerce').astype(int)
    phen_id = sys.argv[4]
    main(phen_df, population).assign(
        phen_id=phen_id,
        N=n_samples
    )[result_columns].to_csv(
        sys.argv[1], sep='\t',
        index=False,
    )