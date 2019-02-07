#! /usr/bin/env python
import os
import sys
import argparse
from collections import namedtuple
import time
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='remove strand from frequency file')
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i','--input',type=argparse.FileType('r'),
            required=False,default=sys.stdin,help="input methylation frequency file")
    args = parser.parse_args()
    return args

methylFreq = namedtuple('methylFreq',['chrom','pos','strand','methylated','unmethylated','di','tri'])
if __name__=="__main__":
    args=parseArgs()
    for line in args.input :
        fields = line.strip().split("\t")
        new = methylFreq(fields[0],int(fields[1]),
                fields[2],int(fields[3]),int(fields[4]),fields[5],fields[6])
        try : entry.chrom
        except ( NameError,AttributeError ):
            entry = new
            continue
        if ( new.pos-1 == entry.pos and
                entry.strand == "+" and
                new.strand == "-" ):
            entry = entry._replace(methylated = entry.methylated + new.methylated,
                    unmethylated =entry.unmethylated + new.unmethylated,
                    strand = "." )
            new = None
        print("\t".join([str(x) for x in list(entry)]))
        entry = new
