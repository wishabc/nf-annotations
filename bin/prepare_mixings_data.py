import numpy as np
import sys

from nmf_tools.component_labels.labels_at_tr import compute_labels_absolute, create_label_matrix



def main(H):
    H_labels, _ = compute_labels_absolute(H, absolute_threshold=0.05, purity_threshold=0.5)

    matrix_for_enrichment, uq_labels = create_label_matrix(H_labels)
    uq_labels, counts = np.unique(H_labels, return_counts=True)

    abundant_labels_mask = (counts >= 2000) & (uq_labels != '')

    matrix_for_enrichment = matrix_for_enrichment[abundant_labels_mask, :].todense().T
    return matrix_for_enrichment, uq_labels[abundant_labels_mask]


if __name__ == "__main__":
    H = np.load(sys.argv[1]).T # comp x DHSs
    matrix_name = sys.argv[2]
    order_name = sys.argv[3]
    annotation, names = main(H)
    assert annotation.shape[0] == H.shape[1]
    np.save(matrix_name, annotation)
    np.savetxt(order_name, names, fmt='%s')
