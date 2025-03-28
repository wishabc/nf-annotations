import sys
import pandas as pd


matching_cols = ['tss_bin', 'ld_bin', 'maf_bin']


def sample_random(df: pd.DataFrame, count_to_sample: pd.Series, seed):
    if df.name not in count_to_sample.index:
        return []
    n = count_to_sample.loc[df.name]
    return df.sample(n=n, random_state=seed).index


def main(ref_pop_df: pd.DataFrame, count_to_sample, seed):
    sampled = ref_pop_df[matching_cols].groupby(matching_cols).apply(
        sample_random,
        count_to_sample=count_to_sample,
        seed=seed
    )
    return ref_pop_df.loc[sampled.values, ['chrom', 'start', 'end', *matching_cols]].sort_values(['chrom', 'start'])



if __name__ == "__main__":
    ref_pop = pd.read_table(sys.argv[1]).set_index('rs_id')
    sampling_counts = pd.read_table(sys.argv[2]).set_index(matching_cols)['count']
    sampling_seed = int(sys.argv[3])

    sampled_df = main(ref_pop, count_to_sample=sampling_counts, seed=sampling_seed)
    sampled_df = sampled_df.reset_index()[[*sampled_df.columns, 'rs_id']]
    sampled_df.to_csv(
        sys.argv[4],
        sep='\t',
        index=False
    )