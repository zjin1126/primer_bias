#!/usr/bin/env python
"""
convert position according to reference
"""
import argparse
import pickle

# import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pos_info', help='Position info generate by extract_posision.py')
    parser.add_argument('usearch_result', help='USEARCH search_oligodb userout result')
    args = parser.parse_args()

    # read position information
    with open(args.pos_info, mode='rb') as fh:
        pos_info = pickle.load(fh)

    # read search_oligodb result
    with open(args.usearch_result) as fh:
        for line in fh:
            fields = line.rstrip().split('\t')
            name = fields[0].split()[0]
            primer_name = fields[1]
            start = int(fields[6])
            end = int(fields[7])
            if start >= len(pos_info[name]):
                start = -1
            else:
                start = pos_info[name][start]
            if end >= len(pos_info[name]):
                end = -1
            else:
                end = pos_info[name][end]
            print('{}\t{}\t{}\t{}'.format(name, primer_name, start, end))


if __name__ == '__main__':
    main()
