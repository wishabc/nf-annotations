import sys
import numpy as np
import scipy.sparse as sp


def split_matrix(matrix, sample_names, prefix):
    for i in range(matrix.shape[1]):
        sample_name = sample_names[i]
        column = matrix.getcol(i) 
        np.savetxt(f'{prefix}.{sample_name}.txt', column.toarray(), fmt='%d')


if __name__ == '__main__':
    name = sys.argv[1]
    if name.endswith('.npy'):
        matrix = sp.csr_matrix(np.load(name))
    elif name.endswith('.npz'):
        matrix = sp.load_npz(name)
    else:
        raise ValueError("Unsupported file format. Expected formats are .npy or .npz")
    sample_names = np.loadtxt(sys.argv[2], dtype=str)
    prefix = sys.argv[3]
    split_matrix(matrix, sample_names, prefix)
