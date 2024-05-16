import sys
import numpy as np


def split_matrix(matrix, sample_names, prefix):
    for i in range(matrix.shape[1]):
        sample_name = sample_names[i]
        np.savetxt(f'{prefix}.{sample_name}.txt', matrix[:, i], fmt='%d')


if __name__ == '__main__':
    matrix = np.load(sys.argv[1])
    sample_names = np.loadtxt(sys.argv[2], dtype=str)
    prefix = sys.argv[3]
    split_matrix(matrix, sample_names, prefix)
