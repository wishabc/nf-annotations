import sys
import pandas as pd
import os
from parse_variants_motifs import read_pfm
from tqdm import tqdm

def choose_best(data, motif_id, tr):
    key, value = data
    chrom, end, ref, alt = key.split('@')
    end = int(end)
    best_allele, best_pval, best_pos, best_orient = max(value, key=lambda x: x[1]) # allele, p-value, position, orientation
    res = [x for x in value
        if (x[0] != best_allele) and (x[2] == best_pos) and (x[3] == best_orient)]
    assert len(res) == 1
    _, other_pval, *_ = res[0]
    stats = [best_pos, 1 if best_pval >=tr else 0, best_orient]
    if best_allele == 'ref':
        stats.extend([best_pval, other_pval])
    else:
        stats.extend([other_pval, best_pval])

    return [chrom, end - 1, end, ref, alt, motif_id, *stats]


def main(sarus_log, motif_length, window_size, out_file, motif_path, tr=4):
    result = {}
    
    motif_id, pfm = read_pfm(motif_path)
    print('Processing started')
    with open(sarus_log) as f:
        for line in tqdm(f.readlines()):
            if line.startswith('>'):
                key, allele = line.strip('\n')[1:].rsplit('@', 1)
            else:
                motif_score, motif_pos, motif_orient = line.strip('\n').split()
                difference = window_size + 1 - motif_length
                if difference < 0:
                    raise AssertionError
                motif_pos = int(motif_pos)
                # max motif pos

                motif_pos = motif_pos - difference if line[2] == '-' else window_size - motif_pos 
                if motif_pos < 0 or motif_pos >= motif_length:
                    continue
                motif_pos = window_size - motif_pos
                motif_score = float(motif_score)
                result.setdefault(key, []).append([allele, motif_score, motif_pos, motif_orient])
    
 
    pd.DataFrame([choose_best(x, motif_id, tr) for x in tqdm(result.items())], 
        columns=['#chr', 'start', 'end', 'ref', 'alt', 'motif',
                'motif_pos', 'signif_hit', 'motif_orient', 'motif_logpval_ref', 'motif_logpval_alt',]
        ).to_csv(out_file, sep='\t', index=False)
            


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5])