#!/bin/usr/python

import sys
from scipy.sparse import coo_matrix, csr_matrix, save_npz
import logging
import datatable as dt
import argparse
import numpy as np

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
logger.setLevel('INFO')
formatter = logging.Formatter('%(asctime)s  %(levelname)s  %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def read_matrix(input_file, dtype):
    # Make sure to say there is NO header in this file. Otherwise will be off by 1
    df = dt.fread(input_file, header=False)
    logger.info('Converting to numpy array')
    return df.to_numpy().astype(dtype)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts matrix in txt format to csr')
    parser.add_argument('matrix', help='Path to matrix file')
    parser.add_argument('outpath', help='Path to output binary matrix file with .npz extension')
    args = parser.parse_args()
    input_path = args.matrix
    out_path = args.outpath
    logger.info('Starting processing')
    matrix_dense = read_matrix(input_path, dtype=bool)
    np.save(out_path, matrix_dense)
    # logger.info(f'Matrix size: {matrix_dense.shape}. '
    #             f'Density: {matrix_dense.sum() / matrix_dense.size}'
    #             )
    # convert_to_sparse(matrix_dense, out_path)
