#!/usr/bin/env python3

import sys
import math
import argparse
from Bio import SeqIO
from Bio.Seq import reverse_complement
import wordcount as wc

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Compute dinucleotide force for regions specified '
        'in the coordinate file')
parser.add_argument('fasta_infile', help='FASTA file with the reference')
parser.add_argument('coordinate_file', help='coordinate file')
parser.add_argument('-d', '--dimer', default = "CG",
    help='dimer to compute the force for')
parser.add_argument('-L', '--min-length', type=int, default=100,
     help='minimal length of the region to compute the force for')
parser.add_argument('-s', '--strand', action='store_true',
     help='take the strand of the alignment into account')


args = parser.parse_args()

sd = {rec.id:str(rec.seq) for rec in  SeqIO.parse(args.fasta_infile, "fasta")}

for line in open(args.coordinate_file):
    S = line.split()
    contig = S[0]
    try:
        seq = sd[contig]
    except KeyError:
        continue
    start = max(0, int(S[1]) - 1)
    end = min(len(seq), int(S[2]))
    if end - start < args.min_length:
        continue
    subseq = seq[start:end]
    if args.strand:
        strand = S[3]
        if strand == "-":
            subseq = reverse_complement(subseq)
        elif strand == "+":
            pass
        else:
            raise ValueError("strand needs to be + or - exclusively")
    N_Valid = wc.count_overlapping_words(subseq, 2, normalize=False).sum()
    dForce = wc.DimerForce(subseq, args.dimer) if N_Valid > 0 else math.nan
    sys.stdout.write("{}\t{}\t{}\n".format(line.strip("\n\r"), N_Valid,
            dForce))

