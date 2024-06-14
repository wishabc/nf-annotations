import pandas as pd
import numpy as np
import sys


def main():
    pass


if __name__ == "__main__":
    W = np.load(sys.argv[1])
    H = np.load(sys.argv[2])
    samples_order = np.loadtxt(sys.argv[3], dtype=str)
    metadata = pd.read_table(sys.argv[4])
    main()