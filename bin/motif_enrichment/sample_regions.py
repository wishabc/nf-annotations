import pandas as pd
import numpy as np
from tqdm import tqdm
import sys


matching_cols = ['gc_bin', 'length_bin']


def sample_random(df: pd.DataFrame, count_to_sample: pd.Series, seed):
    if df.name not in count_to_sample.index:
        return []
    n = count_to_sample.loc[df.name]
    return df.sample(n=n, random_state=seed).index.tolist()


def main(masterlist_df, regions_pool):
    count_to_sample = masterlist_df[matching_cols].value_counts()

    sampled_indices = {}
    grouped = regions_pool.groupby(matching_cols)

    for name, group in tqdm(grouped, total=len(grouped)):
        if name in count_to_sample.index:
            n = count_to_sample.loc[name]
            for random_state in np.arange(n_samples):
                sampled = group.sample(n=n, random_state=random_state, replace=False).index
                sampled_indices.setdefault(random_state, []).extend(sampled)

    return sampled_indices


if __name__ == "__main__":
    header = ['#chr', 'start', 'end', 'length', 'n_gc', 'gc', 'gc_bin', 'length_bin']
    masterlist_df = pd.read_table(sys.argv[1], names=header)
    regions_pool = pd.read_table(sys.argv[2], names=header)
    motif_indicator = np.loadtxt(sys.argv[3], dtype=bool)
    annotation_indicator = np.loadtxt(sys.argv[4], dtype=bool)

    n_samples = 100
    masterlist_df = masterlist_df[motif_indicator & annotation_indicator]
    sampled_indices = main(masterlist_df, regions_pool, n_samples)
    samples = []
    for random_state, indices in tqdm(sampled_indices.items(), total=len(sampled_indices)):
        samples.append(
            regions_pool.loc[indices, :].eval(f'sample = {random_state}')
        )
    pd.concat(samples).sort_index().to_csv(
        sys.argv[5],
        sep='\t',
        index=False,
        header=False,
    )
    
        
