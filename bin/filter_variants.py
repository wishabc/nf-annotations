import pandas as pd
import sys


def main(bed_file):
    pass

if __name__ == '__main__':
    bed_file = pd.read_table(sys.argv[1])
    print(bed_file.columns)
    exit(1)
    main(bed_file)