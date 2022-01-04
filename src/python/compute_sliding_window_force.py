#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
import wordcount as wc

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Compute dinucleotide force using sliding window')
parser.add_argument('fasta_infile', help='FASTA file with the reference')
parser.add_argument('-L', '--window', type=int, default=3000,
     help='length of the sliding window')
parser.add_argument('-c', '--contig', 
    help='only use this conting in the reference file')
parser.add_argument('-d', '--dimer', default = "CG",
    help='dimer to compute the force for')

args = parser.parse_args()

nt2int = {"A":0, "C":1, "G":2, "T":3, "U":3}
nt1 = nt2int[args.dimer[0].upper()]
nt2 = nt2int[args.dimer[1].upper()]

for rec in SeqIO.parse(args.fasta_infile, "fasta"):
    if args.contig and rec.id != args.contig:
        continue
    seq = str(rec.seq).upper().replace("U", "T")
    for start in range(len(seq) - args.window + 1):
        subseq = seq[start:start+args.window]
        dForce = wc.DimerForce(subseq, nt1, nt2)
        N_ACGT = wc.count_words(subseq, 1, normalize=False).sum()
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(rec.id, start + 1, N_ACGT,
            dForce))
