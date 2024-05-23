import pandas as pd
import sys


if __name__ == '__main__':
    initial_variants = pd.read_table(sys.stdin).set_index('varid')
    index = initial_variants.index.copy()

    hg38_variants = pd.read_table(
        sys.argv[1], header=None,
        names=['hg38_chr', 'hg38_start', 'hg38_end', 'varid']
    ).set_index('varid')
    
    initial_variants = initial_variants.join(hg38_variants, how='left').loc[index].reset_index()[
        ['hg38_chr', 'hg38_start', 'hg38_end', 'varid', *initial_variants.columns]
    ]
    
    initial_variants[['hg38_start', 'hg38_end']] = initial_variants[['hg38_start', 'hg38_end']].astype('Int64')
    initial_variants.to_csv(
        sys.argv[2], sep='\t', index=False
    )