import pandas as pd
import sys
import numpy as np
from scipy.stats import hypergeom


def read_as_np(binary_matrix_path):
    return np.load(binary_matrix_path)


def main(binary_matrix_path, motif_hits_path, motif, sample_ids_path, sample_id, step):
    binary_matrix = read_as_np(binary_matrix_path)
    sample_ids = pd.read_table(sample_ids_path, header=None).iloc[:, 0].tolist()
    sample_id = int(sample_id)
    step = int(step)
    sample_ids = sample_ids[sample_id - 1: min(sample_id + step, len(sample_ids))]
    motif_hits = pd.read_table(motif_hits_path, header=None).to_numpy()
    assert motif_hits.shape[0] == binary_matrix.shape[0]
    M = motif_hits.shape[0]
    N = binary_matrix.sum(axis=0)
    n = motif_hits.sum()
    s = (binary_matrix * motif_hits).sum(axis=0)
    enrichment = -hypergeom.logsf(s - 1, M, n, N)
    logodds = np.log2(s) + np.log2(M + s - n - N) - np.log2(n - s) - np.log2(N - s)
    for sample_id, odds, enrich, N_i, s_i in zip(sample_ids, logodds, enrichment, N, s):
        print(sample_id, motif, s_i, N_i, n, M, enrich, odds, sep='\t')

if __name__ == '__main__':
    main(*sys.argv[1:])