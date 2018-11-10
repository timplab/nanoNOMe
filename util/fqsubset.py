#! /usr/bin/env python
import sys
import os
import argparse
import gzip
from Bio import SeqIO

def parseArgs():
    parser = argparse.ArgumentParser( description='match read names with file names' )
    parser.add_argument('-i','--input',type=os.path.abspath,required=True,
            help="input fq")
    parser.add_argument('-n','--names',type=argparse.FileType('r'),required=False,
            default=sys.stdin,help="read names")
    parser.add_argument('-o','--out',type=argparse.FileType('w'),required=False,
            default=sys.stdout,help="output fq")
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args=parseArgs()
    rnames = [ x.strip() for x in args.names ]
    with gzip.open(args.input,'rt') as fh :
        for record in SeqIO.parse(fh,"fastq") :
            if record.id in rnames :
                print(record.format("fastq").strip(),file=args.out)
