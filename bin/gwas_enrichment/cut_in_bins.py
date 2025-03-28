import numpy as np
import pandas as pd
import sys


histogram_data = pd.read_table(sys.argv[1])

val_counts = histogram_data[['tss_bin', 'ld_bin', 'maf_bin']].value_counts().reset_index()

val_counts.to_csv(sys.argv[2], sep='\t', index=False)