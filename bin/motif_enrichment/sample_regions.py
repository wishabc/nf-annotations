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


def cut_in_bins(df: pd.DataFrame):
    gc_bins = np.linspace(0, 1, 30)
    length_bins = np.linspace(20, 1000, 30)
    df['gc_bin'] = pd.cut(df['gc'], bins=gc_bins)
    df['length_bin'] = pd.cut(df.eval('end - start'), bins=length_bins, include_lowest=True)
    return df

def main(masterlist_df, regions_pool, random_state):
    masterlist_df = cut_in_bins(masterlist_df)
    count_to_sample = masterlist_df[matching_cols].value_counts()

    regions_pool = cut_in_bins(regions_pool)
    sampled_indices = []
    grouped = regions_pool.groupby(matching_cols)

    for name, group in tqdm(grouped, total=len(grouped)):
        if name in count_to_sample.index:
            n = count_to_sample.loc[name]
            sampled = group.sample(n=n, random_state=random_state, replace=False).index
            sampled_indices.extend(sampled)

    return regions_pool.loc[sampled_indices, ['#chr', 'start', 'end', *matching_cols]].sort_values(['#chr', 'start'])


if __name__ == "__main__":
    masterlist_df = pd.read_table(sys.argv[1])
    regions_pool = pd.read_table(sys.argv[2])
    motif_indicator = np.loadtxt(sys.argv[3], dtype=bool)
    seed = int(sys.argv[4])
    masterlist_df = masterlist_df[motif_indicator]
    sampled_df = main(masterlist_df, regions_pool, random_state=seed)
    sampled_df.to_csv(sys.argv[5], sep='\t', index=False, header=False)