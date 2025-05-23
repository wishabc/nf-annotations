import pandas as pd
import sys

print('Reading data...')
df = pd.read_table(sys.argv[1], names=['#chr', 'start', 'end', 'length', 'n_gc', 'gc', 'gc_bin', 'length_bin']).drop(columns=['n_gc', 'length', 'gc'])

df['#chr'] = df['#chr'].astype('category')
df['gc_bin'] = df['gc_bin'].astype('category')
df['length_bin'] = df['length_bin'].astype('category')
df = df.sort_values(["gc_bin", "length_bin"])

print('Writing data...')
df.to_parquet(
    sys.argv[2],
    index=False,
    engine='pyarrow',
    compression=None,
    partition_cols=['gc_bin', 'length_bin'],
)

