import pandas as pd
import sys
import numpy as np



gc_bins = np.linspace(0, 1.01, 20) # 1.01 to include 1.0 in bins
length_bins = [*np.linspace(0, 1000, 50), np.inf]
df = pd.read_table(sys.argv[1], names=['#chr', 'start', 'end', 'length', 'n_gc', 'gc'])

df['gc_bin'] = pd.cut(df['gc'], bins=gc_bins, labels=False, include_lowest=True)
df['length_bin'] = pd.cut(df['length'], bins=length_bins, labels=False)

df.to_csv(
    sys.argv[2],
    sep='\t',
    index=False,
    header=False,
)