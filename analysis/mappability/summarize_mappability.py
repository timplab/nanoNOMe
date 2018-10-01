#! /usr/bin/env python
import os
import sys
import argparse
import numpy as np
from collections import namedtuple
import time
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='parse methylation frequency file')
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i','--input',type=argparse.FileType('r'),
            required=False,default=sys.stdin,help="input bedgraph intersect output")
    parser.add_argument('-b','--bins',type=int,
            required=False,default=4,help="number of bins to separate by (default: 4)")
    parser.add_argument('-c','--coverage',type=int,
            required=False,default=2,help="coverage filter thrshold (default: 2)")
    args = parser.parse_args()
    return args

methylFreq = namedtuple('methylFreq',['chrom','pos','strand','methylated','unmethylated'])
if __name__=="__main__":
    args=parseArgs()
    bins = np.linspace(0,1,args.bins+1)
    qdict = { x:0 for x in range(len(bins)) }
    nomap = 0
    i = 0 
    for line in args.input :
        fields = line.strip().split("\t")
        mfreq = methylFreq(fields[0],int(fields[1]),
                fields[2],int(fields[3]),int(fields[4]))
        cov = mfreq.methylated+mfreq.unmethylated 
        if cov < args.coverage : 
            continue
        i+=1
        try : 
            b = np.digitize(float(fields[-1]),bins)
        except ValueError : 
            nomap += 1
            continue
        qdict[b-1] += 1
    print(i)
    print(qdict)
    print(nomap)
#        freq = float(mfreq.methylated)/cov
#        outline = "\t".join([str(x) for x in [fields[-1],cov,freq]])
#        print(outline)

