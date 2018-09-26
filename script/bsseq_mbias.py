#! /usr/bin/env python
import os
import math
import sys
import argparse
import numpy as np
import pysam
import re
import multiprocessing as mp
from multiprocess_utils import listener,init_mp,close_mp
import time
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='parse methylation around SVs')
    parser.add_argument('-t','--threads',type=int,required=False,default=1, 
            help="number of parallel processes (default : 1 )")
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help="bam file")
    parser.add_argument('-o','--output',type=str,required=False, 
            default = "stdout",help="output path (default : stdout)")
    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

def main() :
    args=parseArgs()
    if args.verbose : print("using {} parallel processes".format(args.threads),file=sys.stderr)
    manager,q,pool = init_mp(args.threads)
    # watcher for output
    watcher = pool.apply_async(listener,(q,args.output,args.verbose))
    jobs = list()
    with pysam.AlignmentFile(args.bam,'rb') as bam :
        for entry in bam.fetch() :
            print(entry.to_string())
#            jobs = pool.apply_async(nomeseq_mbias,
#                    args = (entry,args.bam,args.cpg,args.gpc,args.window,args.verbose,q))
#    output = [ p.get() for p in jobs ]
    close_mp(q,pool)
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)

if __name__=="__main__":
    main()

