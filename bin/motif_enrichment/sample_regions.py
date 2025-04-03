import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import time

def sample_random(df: pd.DataFrame, count_to_sample: pd.Series, seed):
    if df.name not in count_to_sample.index:
        return []
    n = count_to_sample.loc[df.name]
    return df.sample(n=n, random_state=seed).index.tolist()


def main(masterlist_df, regions_pool_path, n_samples):
    count_to_sample = masterlist_df[['gc_bin', 'length_bin']].value_counts()

    sampled_data = []

    for index, n in tqdm(count_to_sample.items(), total=len(count_to_sample)):
        data = pd.read_parquet(regions_pool_path, filters=[('gc_bin', '==', index[0]), ('length_bin', '==', index[1])])
        sampled = data.sample(n=n * n_samples, random_state=random_state, replace=True)
        sampled['sample_id'] = np.repeat(np.arange(n_samples), n)
        sampled_data.append(sampled)
            
    return pd.concat(sampled_data, ignore_index=True)


if __name__ == "__main__":
    header = ['#chr', 'start', 'end', 'length', 'n_gc', 'gc', 'gc_bin', 'length_bin']
    masterlist_df = pd.read_table(sys.argv[1], names=header)
    masterlist_df['start'] -= 1 
    masterlist_df = masterlist_df.set_index(['#chr', 'start', 'end'])
    regions_pool_path = sys.argv[2]

    annotation_coordinates = pd.read_table(sys.argv[4]).set_index(['#chr', 'start', 'end'])

    annotation_indicator = np.loadtxt(sys.argv[3], dtype=bool)
    annotation_mask = masterlist_df.index.isin(annotation_coordinates.index)
    print(annotation_mask.sum())
    annotation_indicator_mask = np.zeros(len(masterlist_df), dtype=bool)
    annotation_indicator_mask[annotation_mask] = annotation_indicator

    print(f'DHSs with annotation: {annotation_indicator_mask.sum()}')                                               
    print('Filtering data...')
    random_state = 42
    n_samples = 100
    masterlist_df = masterlist_df[annotation_indicator_mask]
    sampled_data = main(masterlist_df, regions_pool_path, n_samples)
    sampled_data['start'] -= 1 # FIXME from top level script
    print('Writing data...')
    sampled_data.sort_values(['#chr', 'start']).to_csv(
        sys.argv[6],
        sep='\t',
        index=False,
        header=False,
    )
    
        
