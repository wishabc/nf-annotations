import sys
import pandas as pd
import os

def main(log_iterator, motif_length, window_size, out_file, motif):
    result = []
    chrom = end = ref = alt = suffix = None
    for line in log_iterator:
        if line.startswith('>'):
            chrom, end, ref, alt, suffix = line.strip()[2:].split('@')
            end = int(end)
        else:
            motif_score, motif_pos, motif_orient = line.strip('\n').split()
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

            result.append(
                [chrom, end - 1, end, ref, alt, suffix, motif_score, motif_pos, motif_orient, motif]
            )
    
    pd.DataFrame(result, 
        columns=['#chr', 'start', 'end', 'ref', 'alt', 'suffix', 
                'motif_score', 'motif_pos', 'motif_orient', 'motif']
    ).to_csv(out_file, sep='\t', index=False)
            


if __name__ == '__main__':
    motif_name = os.path.splitext(os.path.basename(sys.argv[4]))[0]
    main(sys.stdin, int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], motif_name)