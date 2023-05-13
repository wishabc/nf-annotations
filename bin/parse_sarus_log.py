import sys
import pandas as pd
import os
from parse_variants_motifs import read_pfm
from tqdm import tqdm

def choose_best(data, motif_id):
    key, value = data
    chrom, end, ref, alt = key.split('@')
    end = int(end)
    best_elem = max(value, key=lambda x: x[1]) # allele, p-value, position, orientation
    return [[chrom, end - 1, end, ref, alt, motif_id, *x] for x in value if (x[2] == best_elem[2]) and (x[3] == best_elem[3])]


def main(sarus_log, motif_length, window_size, out_file, motif_path):
    result = {}
    
    motif_id, pfm = read_pfm(motif_path)
    print('Processing started')
    with open(sarus_log) as f:
        for line in tqdm(f.readlines()):
            if line.startswith('>'):
                key, allele = line.strip()[1:].rsplit('@', 1)
            else:
                try:
                    motif_score, motif_pos, motif_orient = line.strip('\n').split()
                except:
                    print(line)
                    raise
                difference = window_size - motif_length
                if difference < 0:
                    raise AssertionError
                motif_pos = int(motif_pos)
                # max motif pos

                motif_pos = motif_pos if line[2] == '-' else window_size - 1 - motif_pos 
                if motif_pos > motif_length:
                    continue
                motif_pos -= window_size
                motif_score = float(motif_score)
                result.setdefault(key, []).append([allele, motif_score, motif_pos, motif_orient])
    
 
    pd.DataFrame([choose_best(x, motif_id) for x in tqdm(result.items())], 
        columns=['#chr', 'start', 'end', 'ref', 'alt', 'allele', 
                'motif_logpval', 'motif_pos', 'motif_orient', 'motif']
        ).to_csv(out_file, sep='\t', index=False)
            


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5])