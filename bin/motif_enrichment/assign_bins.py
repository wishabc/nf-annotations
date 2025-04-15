import pandas as pd
import sys
import numpy as np



gc_bins = np.linspace(0, 1.01, 20) # 1.01 to include 1.0 in bins

df = pd.read_table(sys.argv[1], names=['#chr', 'start', 'end', 'length', 'n_gc', 'gc'])
arange_params = [int(x) for x in sys.argv[2].split(',')]

length_bins = [*np.arange(*arange_params), np.inf]
df['gc_bin'] = pd.cut(df['gc'], bins=gc_bins, labels=False, include_lowest=True)
df['length_bin'] = pd.cut(df['length'], bins=length_bins, labels=False)

df.to_csv(
    sys.argv[3],
    sep='\t',
    index=False,
    header=False,
)