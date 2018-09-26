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
    parser = argparse.ArgumentParser(description='parse methylation frequency file')
    subparsers = parser.add_subparsers()
    # parent
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-t','--threads',type=int,required=False,default=1, 
            help="number of parallel processes (default : 1 )")
    parent_parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parent_parser.add_argument('-i','--input',type=os.path.abspath,required=False,
            help="input methylation frequency file")
    parent_parser.add_argument('-o','--output',type=str,required=False, 
            default="stdout",help="output file (default stout)")
    # methylation by region
    parser_region_meth = subparsers.add_parser('by-region',parents=[parent_parser],
            help = 'get methylation summary by region')
    parser_region_meth.add_argument('-r','--regions',type=argparse.FileType('r'),
            required=False,default=sys.stdin, help="windows in bed format (default stdin)")
    parser_region_meth.add_argument('-c','--coverage',type=int,required=False,default=2,
            help = "coverage threshold (default : 2)")
    parser_region_meth.add_argument('-e','--exclude',type=str,required=False,default="GCG",
            help = "motif to exclude (default :  GCG)")
    parser_region_meth.set_defaults(func=RegionMeth_main)
    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

methFreq = namedtuple('methFreq',
        ['chrom','pos','strand','methylated','unmethylated','dinuc','trinuc'])
def make_key(bedentry) :
    fields=bedentry.strip().split("\t")
    return fields[0]+":"+fields[1]+"-"+fields[2]

def tabix_mfreq(fpath,window) :
    with pysam.TabixFile(fpath,'r') as tabix :
        entries = [ methFreq._make(x.strip().split('\t')) for x in tabix.fetch(window)]
    return entries

def getRegionMeth(infile,window,regline,coverage,exclude,verbose,q) :
    if verbose : print("reading {}".format(window),file=sys.stderr)
    try : 
        mfreq_list = tabix_mfreq(infile,window)
    except : 
        return 
    if not mfreq_list : return
    # filter out exclude
    meth_list = [ [int(y) for y in 
        [x.pos,x.methylated,x.unmethylated]] for x in 
        mfreq_list if x.trinuc != exclude]
    # coverage filter
    meth_list = [ x for x in meth_list if sum(x[1:3]) >= coverage ]
    meth_array = np.array(meth_list)
    if len(meth_array) == 0 : return
    # accounting for strands
    diffs = np.diff(meth_array[:,0])
    coord_diff = np.where(diffs==1)[0]
    if len(coord_diff) > 0 : 
        double_array = [[meth_array[i,1]+meth_array[i+1,1],
            meth_array[i,2]+meth_array[i+1,2]]
            for i in coord_diff]
        double_ind = np.append(coord_diff,coord_diff+1)
        nodouble = np.delete(meth_array,double_ind,axis=0)[:,1:3]
        meth_array = np.append(nodouble,double_array,axis=0)
    else : 
        meth_array = meth_array[:,1:3]
    numsites = len(meth_array)
    num_meth = sum(meth_array[:,0])
    num_unmeth = sum(meth_array[:,1])
    totcov = sum([num_meth,num_unmeth])
    freq = num_meth/totcov
    out = ( regline + 
            '\t'.join([str(x) for x in ['',freq,totcov,numsites]]))
    q.put((window,out))

def RegionMeth_main(args) :
    if args.verbose : print( "parsing regions",file=sys.stderr)
    windows = { make_key(x):x.strip() for x in args.regions }
    if args.verbose : print("{} regions".format(len(windows)),file=sys.stderr)
    if args.verbose : 
        print("using {} parallel processes to get methylation per region".format(args.threads),
                file=sys.stderr)
    # initialize mp
    manager,q,pool = init_mp()
    # watcher for output
    watcher = pool.apply_async(listener,(q,args.output,args.verbose))
    jobs = list()
#    jobs = [ getRegionMeth(args.input,win,reg,args.coverage,args.exclude,args.verbose,q)
#            for win,reg in windows.items() ] 
    jobs = [ pool.apply_async(getRegionMeth, 
            (args.input,win,reg,args.coverage,args.exclude,args.verbose,q))
            for win,reg in windows.items() ] 
    output = [ p.get() for p in jobs ]
    close_mp(q,pool)
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)

if __name__=="__main__":
    args=parseArgs()
    args.func(args)
