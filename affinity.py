#!/usr/bin/env python

import sys
from pynucleic.thermo import *
from Bio import Seq
from itertools import product


def extend_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


with open(sys.argv[1]) as fh:
    print('\t'.join(['Name', 'dH', 'dS', 'dG']))
    for line in fh:
        seq_id, primer_seq = line.rstrip().split()
        dmodel = DNAModel()
        resolved_seq = extend_ambiguous_dna(primer_seq)
        sum_dH = sum_dS = sum_dG = 0
        for seq in resolved_seq:
            dH, dS, dG = dmodel.calc_thermo(Duplex(seq))
            sum_dH += dH
            sum_dS += dS
            sum_dG += dG
        resolved_seq_num = len(resolved_seq)
        print('{}\t{:.4}\t{:.4}\t{:.4}'.format(
            seq_id,
            sum_dH / resolved_seq_num,
            sum_dS / resolved_seq_num,
            sum_dG / resolved_seq_num
        ))
