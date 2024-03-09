import pandas as pd
import sys


if __name__ == '__main__':
    initial_variants = pd.read_table(sys.stdin).set_index('varid')
    hg38_variants = pd.read_table(
        sys.argv[1], 
        names=['hg38_chr', 'hg38_start', 'hg38_end', 'varid']
    ).set_index('varid')
    index = hg38_variants.index.copy()
    initial_variants.join(hg38_variants, how='left')[index].reset_index()[
        ['hg38_chr', 'hg38_start', 'hg38_end', 'varid', *initial_variants.columns]
    ].to_csv(
        sys.argv[2], sep='\t', index=False
    )