#!/usr/bin/env python

#A tool that checks sequence complexity 

import StringIO ## for Python 2

import gzip
import sys
import Bio

def get_complexity(s):
    
    out = StringIO.StringIO()
    
    with gzip.GzipFile(fileobj=out, mode="w") as f:
        f.write(str(s).upper() )
    return len(out.getvalue())/ ( len(s) * 0.641)  #obtained empirically for ran seqs

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print ('Usage: %s sequence' % (sys.argv[0]) )
        print ('Output: complexity of the sequence specified as an input'  )
        sys.exit(1)

    sequence = sys.argv[1]
    print(get_complexity(sequence))