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
    n_samples = 100
    masterlist_df = masterlist_df[motif_indicator]
    sampled_indices = main(masterlist_df, regions_pool, n_samples)
    for random_state, indices in sampled_indices.items():
        sampled_df = regions_pool.loc[indices, ['#chr', 'start', 'end', *matching_cols]].sort_values(['#chr', 'start'])
        sampled_df = sampled_df.reset_index()[[*sampled_df.columns, 'rs_id']]
        # Save the sampled dataframe to a file
        # The output file name is constructed by appending the random state to the original output file name
        # This way, each sampled dataframe will be saved in a separate file
        # The output file name is passed as the 5th argument to the script
        sampled_df.to_csv(f"{sys.argv[4]}_{random_state}.txt", sep='\t', index=False)
    sampled_df.to_csv(sys.argv[5], sep='\t', index=False, header=False)