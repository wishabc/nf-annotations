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


def main(masterlist_df, regions_pool_path):
    count_to_sample = masterlist_df[matching_cols].value_counts()

    sampled_data = []

    for index, n in tqdm(count_to_sample, total=len(count_to_sample)):
        data = pd.read_parquet(regions_pool_path, filters=[('gc_bin', '==', index[0]), ('length_bin', '==', index[1])])
        for random_state in np.arange(n_samples):
            sampled = data.sample(n=n, random_state=random_state, replace=False).eval(f'sample_id = {random_state}')
            sampled_data.append(sampled)
            
    return sampled_data


if __name__ == "__main__":
    header = ['#chr', 'start', 'end', 'length', 'n_gc', 'gc', 'gc_bin', 'length_bin']
    masterlist_df = pd.read_table(sys.argv[1], names=header)
    masterlist_df['start'] -= 1
    masterlist_df = masterlist_df.set_index(['#chr', 'start', 'end'])
    regions_pool_path = sys.argv[2]
    motif_indicator = np.loadtxt(sys.argv[3], dtype=bool)

    annotation_coordinates = pd.read_table(sys.argv[5]).set_index(['#chr', 'start', 'end'])

    annotation_indicator = np.loadtxt(sys.argv[4], dtype=bool)
    annotation_mask = masterlist_df.index.isin(annotation_coordinates.index)
    print(annotation_mask.sum())
    annotation_indicator_mask = np.zeros(len(masterlist_df), dtype=bool)
    annotation_indicator_mask[annotation_mask] = annotation_indicator
                                                     
    print('Filtering data...')
    n_samples = 100
    masterlist_df = masterlist_df[motif_indicator & annotation_indicator_mask]
    sampled_data = main(masterlist_df, regions_pool_path, n_samples)
    print('Writing data...')
    sampled_data = pd.concat(sampled_data)
    sampled_data.sort_values(['#chr', 'start']).to_csv(
        sys.argv[6],
        sep='\t',
        index=False,
        header=False,
    )
    
        
