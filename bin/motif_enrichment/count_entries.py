import pandas as pd
import sys
import numpy as np


# python3 $moduleDir/bin/motif_enrichment/count_entries.py \
#     ${bed_file} \
#     mask.txt \
#     ${motif_id} \
#     ${prefix} \
#     ${name}

bed_file = pd.read_table(sys.argv[1], names=['#chr', 'start', 'end', 'sampling'])
mask = np.loadtxt(sys.argv[2], dtype=bool)
bed_file['motif_hits'] = mask
stats = bed_file.groupby('sampling').agg(
    total=('motif_hits', 'count'),
    motif_hits=('motif_hits', 'sum')
)
stats = stats.reset_index()

stats['motif_id'] = sys.argv[3]
stats['prefix'] = sys.argv[4]
stats.to_csv(
    sys.argv[5],
    sep='\t',
    header=True,
    index=False,
)