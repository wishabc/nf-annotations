import pandas as pd
import sys


if __name__ == '__main__':
    header = pd.read_table(sys.argv[1], nrows=0).columns.tolist()
    df = pd.read_table(sys.stdin, header=None, names=header).dropna(
        subset=[
            'tss_bin', 'ld_bin', 'maf_bin'
        ]
    )
    df['sampling_type'] = 'ref'

    df.to_csv(sys.stdout, sep='\t', index=False, header=True)
