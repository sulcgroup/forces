#!/usr/bin/env python3

"""
A tool that checks sequence complexity
"""

from io import BytesIO

import gzip
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Estimate sequence complexity using gzip; '
    'output: complexity of the input sequence')
parser.add_argument('sequence', help='sequence to estimate complexity of')

args = parser.parse_args()

def get_complexity(seq):
    out = BytesIO()
    bseq = bytes(str(seq.upper()), encoding="ASCII")
    with gzip.GzipFile(fileobj=out, mode="wb") as out_file:
        out_file.write(bseq)
    #obtained empirically for ran seqs
    return len(out.getvalue())/ ( len(seq) * 0.641)

print(get_complexity(args.sequence))
