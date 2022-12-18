import pandas as pd
import sys
import os


def main(base_annotations, custom_annotation, output):
    base_df = pd.read_table(base_annotations)[['CHR', 'BP', 'SNP', 'CM']]
    base_df['index'] = "chr" + base_df['CHR'].astype(str) + '_' + base_df['BP'].astype(str)
    key = os.path.basename(os.path.splitext(custom_annotation)[0])
    base_df[key] = 0
    
    custom_df = pd.read_table(custom_annotation)[['#chr', 'end']]
    custom_df.drop_duplicates()
    custom_df['index'] = custom_df['#chr'] + '_' + custom_df['end'].astype(str)
    
    base_df.loc[base_df['index'].isin(custom_df['index']), key] = 1

    base_df.to_csv(output, sep='\t', index=False)



if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    