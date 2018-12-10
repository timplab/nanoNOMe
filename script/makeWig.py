#! /usr/bin/env python
import os
import math
import sys
import argparse
import gzip
import numpy as np
from collections import namedtuple
import pysam
import re
import multiprocessing as mp
import time
from multiprocess_utils import listener,init_mp,close_mp
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='make wig file')
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i','--input',type=str,required=False,default="stdin",
            help="input file (methylation freq, default stdin)")
    parser.add_argument('-o','--output',type=argparse.FileType('w'),required=False, 
            default=sys.stdout,help="output file (default stout)")
    parser.add_argument('-t','--threshold',type=int,required=False,
            default=3,help="coverage thresdhold (default 3)")
    parser.add_argument('-c','--coverage',type=str,required=False,
            help="output path for coverage file (optional)")
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

if __name__=="__main__":
    args=parseArgs()
    if args.coverage :
        covout = open(args.coverage,'w')
        print("track type=wiggle_0",file=covout)
        def printfxn(prechrom,chrom,pos,meth,cov,out,covout) :
            if chrom != prechrom :
                prechrom = fields[0]
                print("variableStep chrom={}".format(chrom),file=out)
                print("variableStep chrom={}".format(chrom),file=covout)
            print("{}\t{}".format(pos,meth/cov),file=out)
            print("{}\t{}".format(pos,cov),file=covout)
            return prechrom
    else :
        covout = None
        def printfxn(prechrom,chrom,pos,meth,cov,out,covout) :
            if chrom != prechrom :
                prechrom = fields[0]
                print("variableStep chrom={}".format(chrom),file=out)
            print("{}\t{}".format(pos,meth/cov),file=out)
            return prechrom
    if args.input == "stdin" :
        if args.verbose : print("reading from stdin",file=sys.stderr)
        in_fh = sys.stdin
    elif args.input.split(".")[-1] == "gz" :
        if args.verbose : 
            print("reading gzipped input {}".format(args.input),file=sys.stderr)
        in_fh = gzip.open(args.input,'rt')
    else :
        if args.verbose : print("reading {}".format(args.input),file=sys.stderr)
        in_fh = open(args.input)
    i = 0
    chrom=""
    # header
    print("track type=wiggle_0",file=args.output)
    for line in in_fh :
        fields = line.strip().split("\t")
        meth = int(fields[3])
        unmeth = int(fields[4])
        cov = meth+unmeth
        if cov < args.threshold : continue
        chrom = printfxn(chrom,fields[0],fields[1],meth,cov,args.output,covout)
        i += 1
        if args.verbose : 
            if i % 100000 == 0 :
                print("processed {} lines".format(i),file=sys.stderr)
            
    in_fh.close()



