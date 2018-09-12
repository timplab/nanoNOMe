#! /usr/bin/env python
import os
import math
import sys
import argparse
import gzip
import numpy as np
from collections import namedtuple
from methylbed_utils import MethRead
import pysam
import re
import multiprocessing as mp
import time
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='convert bam to be igv compatible')
    parser.add_argument('-t','--threads',type=int,required=False,default=1, 
            help="number of parallel processes (default : 1 )")
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help="bam file - sorted and indexed")
    parser.add_argument('-c','--cpg',type=os.path.abspath,required=True,
            help="gpc methylation bed - sorted, bgzipped, and indexed")
    parser.add_argument('-g','--gpc',type=os.path.abspath,required=False,
            default=None,help="gpc methylation bed - sorted, bgzipped, and indexed")
    parser.add_argument('-w','--window',type=str,required=False, 
            help="window from index file to extract [chrom:start-end]")
    parser.add_argument('-r','--regions',type=argparse.FileType('r'),required=False, 
            default = sys.stdin, help="windows in bed format (default: stdin)")
    parser.add_argument('-o','--out',type=str,required=False,default="stdout",
            help="output bam file (default: stdout)")
    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

def make_key(bedentry) :
    fields=bedentry.strip().split("\t")
    return fields[0]+":"+fields[1]+"-"+fields[2]

def read_tabix(fpath,window) :
    with pysam.TabixFile(fpath) as tabix :
        entries = [x for x in tabix.fetch(window)]
    reads = [MethRead(x) for x in entries]
    rdict = dict((x.qname,x) for x in reads)
    return rdict

# https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def listener(q,inbam,outbam,verbose=False) :
    '''listens for messages on the q, writes to file. '''
    if verbose : print("writing output to {}".format(outbam),file=sys.stderr)
    if outbam == "stdout" :
        # write output to stdout
        f = sys.stdout
        with pysam.AlignmentFile(inbam,'rb') as fh :
            print(fh.header,file=f)
        def printread(m,out_fh) :
            print(m,file=out_fh)
    else : 
        with pysam.AlignmentFile(inbam,'rb') as fh :
            f = pysam.AlignmentFile(outbam, 'wb',template = fh) 
        def printread(m,out_fh) :
            read = pysam.AlignedSegment.fromstring(m,out_fh.header)
            out_fh.write(read)
    qname_list= list()
    while True:
        m = q.get()
        if m == 'kill':
            break
        fields = m.split("\t")
        qname = ':'.join([fields[0],fields[2],fields[3]])
        if qname not in qname_list : 
            if verbose: print(qname,file=sys.stderr)
            qname_list.append(qname)
            printread(m,f)
        else : 
            if verbose: print("{} already written".format(qname),file=sys.stderr)
        q.task_done()
    if verbose : print("total {} reads processed".format(len(qname_list)),file=sys.stderr)
    f.close()
    q.task_done()

def convert_cpg(bam,cpg,gpc) :
    # only cpg
    return change_sequence(bam,cpg,"cpg")

def convert_nome(bam,cpg,gpc) :
    # cpg and gpc
    bam_cpg = change_sequence(bam,cpg,"cpg")
    return change_sequence(bam_cpg,gpc,"gpc")

def reset_bam(bam) :
    refseq = bam.get_reference_sequence()
    bam.query_sequence = refseq
    bam.cigarstring = ''.join([str(len(refseq)),"M"])
    return bam

def change_sequence(bam,methread,mod="cpg") :
    start = bam.reference_start
    pos = bam.get_reference_positions(True)
    calls = methread.callarray
    if bam.is_reverse == True : 
        if mod == "cpg" :
            offset=1
            dinuc = "CN"
        elif mod == "gpc" :
            offset=-1
            dinuc = "NC"
        m="G"
        u="A"
    else : 
        if mod == "cpg" :
            offset=0
            dinuc = "NG"
        elif mod == "gpc" :
            offset = 0
            dinuc = "GN"
        m="C"
        u="T"
    if mod == "cpg" :
        seq = np.array(list(bam.query_sequence.replace("CG",dinuc)))
    elif mod == "gpc" :
        seq = np.array(list(bam.query_sequence.replace("GC",dinuc)))
#    seq[gsites-offset] = g
    # methylated
    meth = calls[np.where(calls[:,0]==1),0]+offset
    seq[np.isin(pos,meth)] = m
    # unmethylated
    meth = calls[np.where(calls[:,1]==0),0]+offset
    seq[np.isin(pos,meth)] = u
    
    bam.query_sequence = ''.join(seq)
    return bam

def convertBam(bampath,cfunc,cpgpath,gpcpath,window,verbose,q) :
    if verbose : print("reading {} from bam file".format(window),file=sys.stderr)
    with pysam.AlignmentFile(bampath,"rb") as bam :
        bam_entries = [x for x in bam.fetch(region=window)]
    bam_dict = dict((x.query_name,x) for x in bam_entries)
    if verbose : print("reading {} from cpg data".format(window),file=sys.stderr)
    cpg_dict = read_tabix(cpgpath,window)
    if verbose : print("reading {} from gpc data".format(window),file=sys.stderr)
    gpc_dict = read_tabix(gpcpath,window)
    if verbose : print("converting bams in {}".format(window),file=sys.stderr)
    i = 0
    for qname in bam_dict.keys() :
        try :
            cpg = cpg_dict[qname]
            gpc = gpc_dict[qname]
        except KeyError : 
            continue
        i += 1
        bam = reset_bam(bam_dict[qname])
        newbam = cfunc(bam,cpg,gpc)
        q.put(newbam.to_string())
    if verbose : print("converted {} bam entries in {}".format(i,window),file=sys.stderr)
    

def main() :
    args=parseArgs()
    windows = [ make_key(x) for x in args.regions ]
    # initialize mp
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    if args.verbose : print("using {} parallel processes".format(args.threads),file=sys.stderr)
    # watcher for output
    watcher = pool.apply_async(listener,(q,args.bam,args.out,args.verbose))
    # which convert function
    if args.gpc is None : converter = convert_cpg
    else : converter = convert_nome
    jobs = [ pool.apply_async(convertBam,
        args = (args.bam,converter,args.cpg,args.gpc,win,args.verbose,q))
        for win in windows ]
    output = [ p.get() for p in jobs ]
    # done
    q.put('kill')
    q.join()
    pool.close()
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)

if __name__=="__main__":
    main()

