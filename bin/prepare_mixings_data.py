import numpy as np
import sys
from tqdm import tqdm


def binary_labels_to_strings(binary_matrix, func=str, sep='_'):
    result = []
    for col in tqdm(range(binary_matrix.shape[1]), total=binary_matrix.shape[1]):
        indices = np.where(binary_matrix[:, col] == 1)[0]
        result.append(sep.join(map(func, indices)))
    return np.array(result)

    
def construct_binary_matrix_for_enrichment(binary_matrix_of_labels, n_dhs, min_count=5_000):
    index_labels = binary_labels_to_strings(binary_matrix_of_labels)
    uq_labels, label_counts = np.unique(index_labels, return_counts=True)
    labels_for_enrichment = uq_labels[label_counts >= min_count]

    mat_for_enrichment = np.zeros((n_dhs, len(labels_for_enrichment)), dtype=bool)
    for i, l in enumerate(labels_for_enrichment):
        mat_for_enrichment[:, i] = index_labels == l

    return mat_for_enrichment, labels_for_enrichment


def construct_binary_matrix_of_labels(H, threshold=0.8):
    k, n = H.shape
    total_sums = H.sum(axis=0)
    sorted_indices = np.argsort(-H, axis=0)
    sorted_H = np.take_along_axis(H, sorted_indices, axis=0)
    cum_sums = np.cumsum(sorted_H, axis=0)
    mask = cum_sums >= (threshold * total_sums)
    binary_matrix = np.zeros((k, n), dtype=int)
    
    for col in tqdm(range(n), total=n):
        # Find the first index where cumulative sum exceeds the threshold
        first_true_index = np.argmax(mask[:, col])
        # Collect the indices up to that point
        selected_indices = sorted_indices[:first_true_index + 1, col]
        # Create binary label
        binary_label = np.zeros(k, dtype=int)
        binary_label[selected_indices] = 1
        binary_matrix[:, col] = binary_label

    return binary_matrix


def main(H):
    clean_annotations = H >= 0.5
    binary_matrix_of_labels = construct_binary_matrix_of_labels(H, threshold=0.8)
    mixing_annotations = construct_binary_matrix_for_enrichment(binary_matrix_of_labels, n_dhs=H.shape[1], min_count=5_000)
    return (clean_annotations.T, np.arange(H.shape[0])), mixing_annotations

if __name__ == "__main__":
    H = np.load(sys.argv[1]).T # comp x DHSs
    H =  H / H.sum(axis=0, keepdims=True) # normalize per DHS
    prefix = sys.argv[2]
    (pure_annotation, pure_order), (mixing_annotation, mixing_order) = main(H)
    assert pure_annotation.shape[0] == H.shape[1] == mixing_annotation.shape[0]
    np.save(f'{prefix}.pure.50pr.npy', pure_annotation)
    np.savetxt(f'{prefix}.pure.50pr.order.txt', pure_order, fmt='%s')

    np.save(f'{prefix}.mixing.80pr.npy', mixing_annotation)
    np.savetxt(f'{prefix}.mixing.80pr.order.txt', mixing_order, fmt='%s')
