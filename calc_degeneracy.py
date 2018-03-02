#!/usr/bin/env python
from Bio import SeqIO

IUPAC_CODE = {
    "A": 1,
    "C": 1,
    "G": 1,
    "T": 1,
    "U": 1,
    "R": 2,
    "Y": 2,
    "S": 2,
    "W": 2,
    "K": 2,
    "M": 2,
    "B": 3,
    "D": 3,
    "H": 3,
    "V": 3,
    "N": 4
}

with open("primers.fa") as fh:
    sequences = SeqIO.parse(fh, "fasta")
    for seq in sequences:
        deg_degree = 1
        for base in seq:
            deg_degree *= IUPAC_CODE[base]
        print("{}\t{}".format(seq.id, deg_degree))
