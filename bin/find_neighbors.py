import pandas as pd
import sys


# assumes df is sorted
def get_pairs(df, group_name):
    # returns indeces of the pairs in dataframe + distance
    combs = list(zip(range(len(df.index) - 1), range(1, len(df.index))))
    for comb in combs:
        assert comb[0] < comb[1]
    if len(combs) == 0:
        return []
    starts, ends = map(list, zip(*combs))
    positions = df.end.to_numpy()
    es = df.es.to_numpy()
    result = [(positions - 1)[starts], positions[starts], positions[ends], positions[ends] - positions[starts], es[starts], es[ends]]
    columns=['pair_first_start', 'pair_first', 'pair_second', 'distance', 'start_es', 'end_es']
    result = pd.DataFrame({y: x for x, y in zip(result, columns)})
    result['chr'] = group_name[0]
    result['sample_id'] = group_name[1]
    return result[['chr', *columns, 'sample_id']]




def main(sample_df):
    groups = sample_df.groupby(['#chr', 'sample_id'])

    pairs = []
    for group_name in list(groups.groups):
        p = get_pairs(groups.get_group(group_name).reset_index(drop=True), group_name)
        if len(p) > 0:
            pairs.append(p)
    return pd.concat(pairs).reset_index(drop=True)

if __name__ == '__main__':
    sample_df = pd.read_table(sys.argv[1])
    result_df = main(sample_df)
    result_df.to_csv(sys.stdout, index=False, sep='\t', header=False)