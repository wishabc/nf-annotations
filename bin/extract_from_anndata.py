from genome_tools.data.anndata import read_zarr_backed
import sys
import numpy as np

if __name__ == "__main__":
    peaks_mask = np.loadtxt(sys.argv[2], dtype=bool)
    anndata = read_zarr_backed(sys.argv[1])[:, peaks_mask]
    binary_matrix = anndata.layers['binary'].todense().T
    np.save(sys.argv[3], binary_matrix)
    np.savetxt(sys.argv[4], anndata.obs_names, fmt='%s')
    anndata.var.to_csv(sys.argv[5], sep='\t', index=False, header=False)