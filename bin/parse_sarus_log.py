import sys
import pandas as pd
from parse_variants_motifs import read_pfm
from tqdm import tqdm
import pyfaidx


def choose_best(data, motif_id, motif_length, fasta, tr):
    key, value = data
    chrom, end, rs_id, ref, alt = key.split('@')
    end = int(end)
    best_allele, best_pval, best_pos, best_orient = max(value, key=lambda x: x[1]) # allele, p-value, position, orientation
    res = [x for x in value
        if (x[0] != best_allele) and (x[2] == best_pos) and (x[3] == best_orient)]
    assert len(res) == 1
    _, other_pval, *_ = res[0]
    stats = [best_pos, 1 if best_pval >= tr else 0, best_orient]
    if best_allele == 'ref':
        stats.extend([best_pval, other_pval])
    else:
        stats.extend([other_pval, best_pval])
    start = end - 1
    seq = fasta[chrom][start - best_pos:end + motif_length - best_pos]
    return [chrom, start, end, rs_id, ref, alt, motif_id, *stats, seq]


def main(sarus_log, fasta_path, motif_length, window_size, out_file, motif_path, tr=4):
    result = {}
    
    motif_id, pfm = read_pfm(motif_path)
    print('Processing started')
    with open(sarus_log) as f:
        for line in f:
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
                motif_score = float(motif_score)
                result.setdefault(key, []).append([allele, motif_score, motif_pos, motif_orient])
    with pyfaidx.Fasta(fasta_path, sequence_always_upper=True) as fasta:
        pd.DataFrame([choose_best(x, motif_id, motif_length, fasta, tr=tr) for x in tqdm(result.items())], 
            columns=['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'motif',
                    'motif_pos', 'signif_hit', 'motif_orient',
                    'motif_logpval_ref', 'motif_logpval_alt', 'seq']
            ).to_csv(out_file, sep='\t', index=False)
            


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]),
        int(sys.argv[4]), sys.argv[5], sys.argv[6])