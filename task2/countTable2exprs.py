import pandas as pd
import sys


def main():
    if len(sys.argv) != 3:
        print("Usage: python " + sys.argv[0]
              + " source_file output_file")
        exit()
    source_filename = sys.argv[1]
    output_filename = sys.argv[2]

    data = pd.read_csv('raw_countTable.csv', sep=',', index_col=0)
    colsum = data.sum(0)

    for col in data:
        data[col] /= colsum[col]

    data = data[(data.T != 0).any()]
    data = data.transpose()
    data.to_csv(output_filename, index=False)


if __name__ == '__main__':
    main()
