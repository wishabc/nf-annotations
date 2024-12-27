import pandas as pd
import numpy as np
import sys
import os

from genome_tools.data.anndata import read_zarr_backed


def main(W, densitity_files, samples_order, topX, suffix):
    assert W.shape[1] == densitity_files.shape[0] == samples_order.shape[0]
    top_samples = []
    for i in range(W.shape[0]):
        major_samples = np.where(W[i, :] == np.max(W, axis=0))[0]
        sorted_samples = major_samples[np.argsort(W[i, major_samples])[::-1]]
        for sample in sorted_samples[:topX]:
            data = densitity_files.iloc[sample]
            ag_id = densitity_files.index[sample]
            os.symlink(
                data,
                f'{i}.{ag_id}.component_{suffix}.bw'
            )
            top_samples.append([ag_id, i])
    return pd.DataFrame.from_records(top_samples, columns=['ag_id', 'component'])

if __name__ == "__main__":
    W = np.load(sys.argv[1]).T
    samples_order = np.loadtxt(sys.argv[2], dtype=str)
    density_tracks = read_zarr_backed(
        sys.argv[3]
    ).obs.loc[samples_order]['normalized_density_bw']

    unique_suffix = sys.argv[4] 
    top = int(sys.argv[5])
    top_samples = main(W, density_tracks, samples_order, topX=top, suffix=unique_suffix)
    top_samples.to_csv(f"{unique_suffix}.top_samples.tsv", index=False, sep="\t")
    
    print(len(sys.argv))
    if len(sys.argv) == 7:
        basepath = f"{sys.argv[6]}/{unique_suffix}"
        tracks_paths = pd.DataFrame({
            'component': np.arange(W.shape[0]),
            'aggregated_bw': [f'{basepath}.{i}.top_samples.bw' for i in range(W.shape[0])],
            'aggregated_bg': [f'{basepath}.{i}.top_samples.bg' for i in range(W.shape[0])],
        })
        tracks_paths['n_samples'] = top_samples['component'].value_counts()
        
        tracks_paths.to_csv(f'{unique_suffix}.density_tracks_meta.tsv', sep='\t', index=False)