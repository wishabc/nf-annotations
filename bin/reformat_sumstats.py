import pandas as pd
import sys
import numpy as np

result_columns = ['#chr', 'start', 'end', 'SNP', 'ref', 'alt', 'Beta', 'Beta_se', 'P', 'neglog10_p', 'INFO', 'phen_id', 'N']


def process_chunk(df, population):
    assert np.all(
        (df['chrom'] == df['chr']) & 
        (df['pos'] == df['pos.1']) & 
        (df['ref'] == df['ref.1']) & 
        (df['alt'] == df['alt.1'])
    )    
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
    df.dropna(subset=['Beta', 'neglog10_p', 'start'], inplace=True)
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    return df

def main(input_file, population):
    chunk_size = 5e5  # Adjust chunk size according to your memory constraints
    chunks = pd.read_table(input_file, dtype={'chr': str, 'chrom': str}, chunksize=chunk_size)

    processed_chunks = []
    for chunk in chunks:
        processed_chunk = process_chunk(chunk, population)
        processed_chunks.append(processed_chunk)
    
    return pd.concat(processed_chunks)


if __name__ == '__main__':
    phen_df = pd.read_table(sys.stdin, dtype={'chr': str, 'chrom': str})
    n_samples = pd.to_numeric(sys.argv[3], errors='coerce').astype(int)
    phen_id = sys.argv[4]
    main(phen_df, sys.argv[2]).assign(
        phen_id=phen_id,
        N=n_samples
    )[result_columns].to_csv(
        sys.argv[1], sep='\t',
        index=False,
    )