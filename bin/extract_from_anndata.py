from genome_tools.data.anndata import read_zarr_backed
import sys
import numpy as np

if __name__ == "__main__":
    anndata = read_zarr_backed(sys.argv[1])
    binary_matrix = anndata.layers['binary'][:, :].todense().T
    np.save(sys.argv[2], binary_matrix)
    np.savetxt(sys.argv[3], anndata.obs_names, fmt='%s')