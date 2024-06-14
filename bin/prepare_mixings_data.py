import numpy as np
import sys
from tqdm import tqdm


def construct_binary_labels(H, threshold=0.8, min_count=5_000):
    k, n = H.shape
    labels = np.zeros((n,), dtype=int)
    total_sums = H.sum(axis=0)
    sorted_indices = np.argsort(-H, axis=0)
    sorted_H = np.take_along_axis(H, sorted_indices, axis=0)
    cum_sums = np.cumsum(sorted_H, axis=0)
    mask = cum_sums >= (threshold * total_sums)
    
    for col in tqdm(range(n), total=n):
        # Find the first index where cumulative sum exceeds the threshold
        first_true_index = np.argmax(mask[:, col])
        # Collect the indices up to that point
        selected_indices = sorted_indices[:first_true_index + 1, col]
        # Create binary label
        binary_label = np.zeros(k, dtype=int)
        binary_label[selected_indices] = 1
        # Convert binary label to an integer
        labels[col] = int("".join(binary_label.astype(str)), 2)
    
    unique_labels, indices, counts = np.unique(labels, axis=0, return_inverse=True, return_counts=True)
    
    sufficient_dhs_labels_mask = counts >= min_count
    filtered_unique_labels = unique_labels[sufficient_dhs_labels_mask]
    
    # Create the binary matrix with unique labels as columns
    filtered_indices = np.array([i for i, index in enumerate(indices) if sufficient_dhs_labels_mask[index]])
    binary_matrix = np.zeros((n, filtered_unique_labels.shape[0]), dtype=bool)
    for i, index in enumerate(filtered_indices):
        binary_matrix[i, index] = 1
    
    return unique_labels, binary_matrix


def main(H):
    clean_annotations = H >= 0.5
    mixing_annotations = construct_binary_labels(H, 0.8)
    return (np.arange(H.shape[0]), clean_annotations), mixing_annotations

if __name__ == "__main__":
    H = np.load(sys.argv[1]).T
    H =  H / H.sum(axis=0, keepdims=True)
    prefix = sys.argv[2]
    clean_ann, mixing_ann = main(H)
    np.save(f'{prefix}.clean.50pr.npy', clean_ann[1])
    np.savetxt(f'{prefix}.clean.50pr.order.txt', clean_ann[0], fmt='%s')

    np.save(f'{prefix}.mixing.80pr.npy', mixing_ann[1])
    np.savetxt(f'{prefix}.mixing.80pr.order.txt', mixing_ann[0], fmt='%s')
