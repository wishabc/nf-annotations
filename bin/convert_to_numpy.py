#!/bin/usr/python
import datatable as dt
import argparse
import numpy as np



def read_matrix(input_file, dtype):
    # Make sure to say there is NO header in this file. Otherwise will be off by 1
    df = dt.fread(input_file, header=False)
    return df.to_numpy().astype(dtype)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts matrix in txt format to npy')
    parser.add_argument('matrix', help='Path to matrix file')
    parser.add_argument('outpath', help='Path to output binary matrix file with .npy extension')
    args = parser.parse_args()
    input_path = args.matrix
    out_path = args.outpath
    matrix_dense = read_matrix(input_path, dtype=bool)
    np.save(out_path, matrix_dense)
