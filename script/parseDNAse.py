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
from multiprocess_utils import init_mp,close_mp
from methylbed_utils import bed_to_coord,tabix
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='parse methylation frequency file')
    subparsers = parser.add_subparsers()
    # parent
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    # methylation by distance to center
    parser_region_meth = subparsers.add_parser('by-distance',parents=[parent_parser],
            help = 'get methylation summary by distance to center of features')
    parser_region_meth.add_argument('-i','--input',type=argparse.FileType('r'),
            required=False,default=sys.stdin,
            help="input DNAse-seq signal-to-region distance file (done by bedtools closest) (default stdin)")
    parser_region_meth.add_argument('-w','--window',type=int,
            required=False,default=2000,help="window")
    parser_region_meth.set_defaults(func=DNAseByDistance_main)
    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

def DNAseByDistance_main(args) :
    if args.verbose : print("parsing",file=sys.stderr)
    agg_dict = dict()
    side = args.window/2
    for line in args.input :
        fields = line.strip().split("\t")
        if fields[4] == "." : continue
        start,end=int(fields[1]),int(fields[2])
        regstart,regend=int(fields[5])-side,int(fields[6])+side
        if start < regstart :
            if end < regend : 
                dist_list = list(range(0,int(end-regstart)))
            else : 
                dist_list = list(range(0,int(args.window)))
        elif end >= regend : 
            if start > regstart : 
                dist_list = list(range(int(args.window-regend+start),int(args.window)))
        else : 
            dist_list = list(range(int(start-regstart),int(args.window-regend+end)))
        dist_signed = [ x-side for x in dist_list ]
        if fields[9] == "-" : 
            dist_signed = [ -x for x in dist_signed ]
        for pos in dist_signed : 
            try : 
                agg_dict[pos].append(float(fields[3]))
            except KeyError : 
                agg_dict[pos] = [float(fields[3])]
    avg_list = [ "\t".join([str(x),str(np.mean(agg_dict[x]))]) for x in agg_dict.keys() ]
    for line in avg_list :
        print(line,file=sys.stdout)
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)

if __name__=="__main__":
    args=parseArgs()
    args.func(args)
