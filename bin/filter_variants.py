import pandas as pd
import sys


def main():


if __name__ == '__main__':
    bed_file = pd.read_table(sys.argv[1], header=None,
        names=['chr', 'start', 'end', 'ref', 'alt',
            'mut_rates_roulette', 'mut_rates_gnomad',
            'chr1', 'start1', 'end1', 'ref', 'alt'
            ])
    main()