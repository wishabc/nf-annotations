import pandas as pd
import numpy as np
import sys
import os


def main(W, densitity_files, samples_order, topX, suffix):
    assert W.shape[1] == densitity_files.shape[0] == samples_order.shape[0]
    top_samples = []
    for i in range(W.shape[0]):
        major_samples = np.where(W[i, :] == np.max(W, axis=0))[0]
        sorted_samples = major_samples[np.argsort(W[i, major_samples])[::-1]]
        top_samples.append(sorted_samples[:topX])
        for sample in sorted_samples:
            data = densitity_files.iloc[sample]
            ag_id = densitity_files.index[sample]
            os.symlink(
                data,
                f'{i}.{ag_id}.component_{suffix}.bw'
            )


if __name__ == "__main__":
    W = np.load(sys.argv[1])
    samples_order = np.loadtxt(sys.argv[2], dtype=str)
    density_tracks = pd.read_table(
        sys.argv[3]
    ).set_index(
        'ag_id'
    )['normalized_density_bw'].loc[samples_order]

    unique_suffix = sys.argv[4] 
    top = int(sys.argv[5])
    main(W, density_tracks, samples_order, topX=top, suffix=unique_suffix)