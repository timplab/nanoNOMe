#! /usr/bin/env python
import sys
import argparse
import os
from collections import namedtuple
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='merge overlapping features in bed')
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i','--input',type=argparse.FileType('r'),required=False,
            default=sys.stdin, help="input bed")
    parser.add_argument('--maxsize',type=int,required=False,default=1000,
            help="max size to consider")
    args = parser.parse_args()
    return args

BedEntry = namedtuple('BedEntry',['chrom','start','end','extcols'])
if __name__=="__main__":
    args = parseArgs()
    if args.verbose : 
        print("cleaning...",file=sys.stderr)
    line = args.input.readline()
    fields = line.strip().split("\t")
    oldentry = BedEntry(fields[0],int(fields[1]),int(fields[2]),
            fields[3:])
    for line in args.input :
        fields = line.strip().split("\t")
        newentry = BedEntry(fields[0],int(fields[1]),int(fields[2]),
                fields[3:])
        if oldentry.end - oldentry.start > args.maxsize :
            # skip oldentry if this is too big
            oldentry = newentry
        elif (newentry.chrom == oldentry.chrom and 
                newentry.start <= oldentry.end ) :
            # merge if overlap
            oldentry = BedEntry(newentry.chrom,oldentry.start,
                    max(oldentry.end,newentry.end),oldentry.extcols)
        else :
            # print otherwise
            print("\t".join([oldentry.chrom,
                str(oldentry.start),
                str(oldentry.end)]+oldentry.extcols),file=sys.stdout)
            oldentry = newentry
    print("\t".join([oldentry.chrom,
        str(oldentry.start),
        str(oldentry.end)]+oldentry.extcols),file=sys.stdout)
