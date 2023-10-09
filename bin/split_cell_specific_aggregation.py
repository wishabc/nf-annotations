import pandas as pd
import argparse

def filter_and_write_bed(pval_file: str, fdr_tr: float, suffix: str = ""):
    # Read the file, considering the first line starting with '#' as header
    df = pd.read_table(pval_file)

    # Filter rows
    filtered_df = df[df['min_fdr'] <= fdr_tr]

    # Group by the 'group_id' column and write to individual .bed files
    for group_id, group_data in filtered_df.groupby('group_id'):
        filename = f"{group_id}{suffix}"
        group_data.to_csv(filename, sep='\t', index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter data and write to .bed files based on group_id.')
    parser.add_argument('pval_file', type=str, help='Path to the pval file.')
    parser.add_argument('--fdr_tr', type=float, default=0.05, help='FDR threshold. Default is 0.05.')
    parser.add_argument('--suffix', type=str, default="", help='Suffix for the output .bed files.')

    args = parser.parse_args()
    filter_and_write_bed(args.pval_file, args.fdr_tr, args.suffix)
