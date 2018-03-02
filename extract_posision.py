#!/usr/bin/env python
"""
Extract position information from aligned sequences
"""
import argparse
import array
import pickle

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref', help='Aligned reference sequence (e.g. E. coli 16S rRNA)')
    parser.add_argument('query', help='Aligned query sequences in fasta format')
    parser.add_argument('out', help='Output file name')
    args = parser.parse_args()

    # read reference sequence and convert to position
    ref_seq = str(SeqIO.read(args.ref, 'fasta').seq)
    pos = -1
    ref_pos = []
    for base in ref_seq:
        if base != '.' and base != '-':
            pos += 1
        ref_pos.append(pos)

    # read query sequences and build position list
    progress_count = 0
    result = {}
    for record in SeqIO.parse(args.query, 'fasta'):
        query_pos = []
        query_seq = str(record.seq)
        for i in range(len(query_seq)):
            if query_seq[i] != '.' and query_seq[i] != '-':
                query_pos.append(ref_pos[i])
        result[record.id] = array.array('h', query_pos)
        progress_count += 1
        if progress_count % 1000 == 0:
            print(progress_count, 'sequences proceded...')

    # write result to disk
    with open(args.out, mode='wb') as fh:
        pickle.dump(result, fh)

if __name__ == '__main__':
    main()
