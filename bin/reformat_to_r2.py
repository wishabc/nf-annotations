import sys

prev_chrom = prev_start = prev_end = None
for variant_line in sys.stdin:
    chrom, start, end = variant_line.strip().split('\t')
    if chrom == prev_chrom:
        print(chrom, prev_start, prev_end, end)
    prev_chrom, prev_start, prev_end = chrom, start, end 