#!/usr/bin/env python
"""
Extract reference sequence by ID (require fai index)
"""

import argparse
import os
import sys

from pyfaidx import Fasta

parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Reference sequences in fasta format')
parser.add_argument('reference_id', help='Reference ID to be extract')
args = parser.parse_args()

# check if fai index exists
if not os.path.exists(args.fasta + '.fai'):
    print('fai index does not exist', file=sys.stderr)
    exit(1)

# load index
reference_seqs = Fasta(args.fasta)
print('>' + reference_seqs[args.reference_id].name)
print(reference_seqs[args.reference_id])
