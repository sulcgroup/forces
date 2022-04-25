#!/usr/bin/env python3

import sys
import math
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
    help='only use this contig in the reference file')
parser.add_argument('-d', '--dimer', default = "CG",
    help='dimer to compute the force for')
parser.add_argument('-e', '--end', type=int,
    help='end at this coordinate in the --contig, 1-based')
parser.add_argument('-s', '--start', type=int,
    help='start at this coordinate in the --contig, 1-based')

args = parser.parse_args()


for rec in SeqIO.parse(args.fasta_infile, "fasta"):
    if args.contig and rec.id != args.contig:
        continue
    seq = str(rec.seq).upper().replace("U", "T")
    if args.start is None or args.end is None:
        start = 0
        end = len(seq) - args.window + 1
    else:
        start = max(0, args.start -1)
        end = min(len(seq) - args.window + 1, args.end)
    for pos in range(start, end):
        subseq = seq[pos:pos+args.window]
        N_Valid = wc.count_overlapping_words(subseq, 2, normalize=False).sum()
        dForce = wc.DimerForce(subseq, args.dimer) if N_Valid > 0 else math.nan
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(rec.id, pos + 1, N_Valid,
            dForce))
