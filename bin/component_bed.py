#!/usr/bin/env python3
import pandas as pd
import argparse
import numpy as np

def convert_to_bed(comb_df):
    for i in range(1, 25):
        comb_df[comb_df[i] == 1].reset_index().iloc[:, 1:5].to_csv(f'comp_{i}.bed',
                                                                       sep='\t', header=False, index=False)

def create_binary_matrix(nmf_mat):
     # Assuming nmf_24 is your existing array
    rows, cols = nmf_mat.shape

    # Initialize a zero matrix of the same shape
    binary_matrix = np.zeros((rows, cols))

    # Iterate over each row and set the max value's index to 1
    for i in range(rows):
        max_index = np.argmax(nmf_mat[i, :])
        binary_matrix[i, max_index] = 1
    
    return binary_matrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate bed files for each component')
    parser.add_argument('metadata', help='Metadata df for the matrix')
    parser.add_argument('nmf_matrix', help='Path to nmf matrix')
    args = parser.parse_args()

    meta_df = pd.read_table(args.metadata, header=None, names=['chr','str','end','chunk_id'])
    binary_nmf = create_binary_matrix(np.load(args.nmf_matrix).T)

    nmf_df = pd.DataFrame(binary_nmf, columns=[i for i in range(1, 25)])

    combined_df = pd.concat([meta_df, nmf_df], axis=1)

    convert_to_bed(combined_df)
