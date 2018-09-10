#! /usr/bin/env python
import os
import math
import sys
import argparse
import gzip
import numpy as np
from methylbed_utils import MethRead,SnifflesEntry
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
    parser = argparse.ArgumentParser(description='parse methylation around SVs')
    parser.add_argument('-t','--threads',type=int,required=False,default=1, 
            help="number of parallel processes (default : 1 )")
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-s','--sniffles',type=argparse.FileType('r'),required=False, 
            default=sys.stdin,help="sniffles vcf output")
    parser.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help="bam file")
    parser.add_argument('-c','--cpg',type=os.path.abspath,required=True,
            help="gpc methylation bed - sorted, bgzipped, and indexed")
    parser.add_argument('-g','--gpc',type=os.path.abspath,required=True,
            help="gpc methylation bed - sorted, bgzipped, and indexed")
    parser.add_argument('--type',type=str,required=False,
            default="TRA",help="type of SV to parse (default: TRA)")
    parser.add_argument('-o','--output',type=str,required=False, 
            default = "stdout",help="output path (default : stdout)")
    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

def make_coord(chrom,start,end) :
    if start < 1 : start = 1
    return chrom+":"+str(start)+"-"+str(end)

def read_tabix(fpath,window) :
    with pysam.TabixFile(fpath) as tabix :
        entries = [x for x in tabix.fetch(window)]
    reads = [MethRead(x) for x in entries]
    rdict = dict()
    for read in reads :
        try :
            rdict[read.qname].append(read)
        except :
            rdict[read.qname] = [read]
    return rdict

def read_bam(fpath,window) :
    with pysam.AlignmentFile(fpath,'rb') as bam :
        bam_entries = [ x for x in bam.fetch(region=window) ]
    bamdict = dict()
    for bam in bam_entries :
        try :
            bamdict[bam.query_name].append(bam)
        except :
            bamdict[bam.query_name] = [bam]
    return bamdict


# https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def listener(q,out,verbose=False) :
    '''listens for messages on the q, writes to file. '''
    if verbose : print("writing output to {}".format(out),file=sys.stderr)
    if out == "stdout" :
        out_fh = sys.stdout
    else :
        out_fh = open(out,'w')
    reg_list= list()
    while 1:
        m = q.get()
        if m == 'kill':
            break
        key,line = m
        if key not in reg_list:
            reg_list.append(key)
            print(line,file=out_fh)
        else :
            if verbose: print("{} already written - skipping".format(key),file=sys.stderr)
        q.task_done()
    if verbose : print("total {} lines written".format(len(reg_list)),file=sys.stderr)
    out_fh.close()
    q.task_done()

def getRegMeth(read_list,start,end) :
    callarray = np.concatenate([ x.callarray for x in read_list])
    regidx = np.where(np.logical_and(callarray[:,0]>=start, callarray[:,0]<=end))
    callreg = callarray[regidx,1].flatten()
    sigidx = np.argwhere(callreg != -1)
    sigreg = callreg[sigidx]
    methcount = np.count_nonzero(sigreg)
    return len(sigreg),methcount
    
def parse_methylation(q,sv,cpg,gpc,start,end,tag) :
    qname = cpg[0].qname
    cpgcov,cpgmeth = getRegMeth(cpg,start,end)
    gpccov,gpcmeth = getRegMeth(gpc,start,end)
    if (gpccov == 0 or cpgcov == 0) : return
    line = '\t'.join([str(x) for x in [ sv.chrom,sv.pos,sv.pos,
        sv.info["CHR2"],sv.info["END"],sv.info["END"],
        qname,sv.id,".",".",tag,cpgmeth,gpcmeth,cpgcov,gpccov]])
    q.put((qname+sv.id+tag,line))

def TRA_methylation(sv,bamfn,cpgfn,gpcfn,verbose,q) :
    win = 300
    methwin = 250
    sv.activate()
    # if majority supports reference, no significant
    if sv.allele == "0/0" : return
    win1 = make_coord(sv.chrom,sv.pos-win,sv.pos+win)
    win2 = make_coord(sv.info["CHR2"],sv.info["END"]-win,sv.info["END"]+win)
    #q.put((win1+win2,win1+win2))
    bam_dicts = [ read_bam(bamfn,w) for w in [win1,win2] ]
    cpg_dicts = [ read_tabix(cpgfn,w) for w in [win1,win2] ]
    gpc_dicts = [ read_tabix(gpcfn,w) for w in [win1,win2] ]
    qnames = list(bam_dicts[0].keys()) + list(bam_dicts[1].keys())
    for qname in set(qnames) :
        if ( qname in bam_dicts[0] and qname in bam_dicts[1] ):
            # this read is an SV and has both parts
            try :
                cpg1 = cpg_dicts[0][qname]
                gpc1 = gpc_dicts[0][qname]
                cpg = cpg_dicts[1][qname]
                gpc = gpc_dicts[1][qname]
                parse_methylation(q,sv,cpg1,gpc1,sv.pos-methwin,sv.pos+methwin,"target")
                start,end,tag = (sv.info["END"]-methwin,sv.info["END"]+methwin,"insert")
                parse_methylation(q,sv,cpg,gpc,start,end,tag)
                continue
            except :
                pass
        if ( qname in bam_dicts[1] ) :
            try :
                cpg = cpg_dicts[1][qname]
                gpc = gpc_dicts[1][qname]
            except :
                continue
            coords = [ pos for x in bam_dicts[1][qname] for pos in [x.reference_start,x.reference_end] ]
            start,end = (sv.info["END"]-methwin,sv.info["END"]+methwin)
            bp = [ x for x in coords if x >= sv.info["END"]-win and x <= sv.info["END"]+win ]
            if len(bp) > 0 : 
                # this read is the insert of SV
                tag = "insert"
            else :
                # this read is insert that has not been transposed
                tag = "insert_noSV"
        elif ( qname in bam_dicts[0] ) :
            try :
                cpg = cpg_dicts[0][qname]
                gpc = gpc_dicts[0][qname]
            except :
                continue
            coords = [ pos for x in bam_dicts[0][qname] for pos in [x.reference_start,x.reference_end] ]
            start,end = (sv.pos-methwin,sv.pos+methwin)
            bp = [ x for x in coords if x >= sv.pos-win and x <= sv.pos+win ]
            if len(bp) > 0 : 
                # this read is the target portion of SV
                tag = "target"
            else :
                # this read is inserted region without SV
                tag = "target_noSV"
        parse_methylation(q,sv,cpg,gpc,start,end,tag)

   

def main() :
    args=parseArgs()
    if args.verbose : print("parsing sniffles file",file=sys.stderr)
    svlines = [SnifflesEntry(x) for x in args.sniffles.readlines() if x[0]!="#"]
    sv_type = [x for x in svlines if x.type == args.type]
    if args.verbose : print("{} SVs selected out of {}".format(len(sv_type),len(svlines)),file=sys.stderr)
    if args.verbose : print("using {} parallel processes".format(args.threads),file=sys.stderr)
    # initialize mp
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    # watcher for output
    watcher = pool.apply_async(listener,(q,args.output,args.verbose))
    if args.type == "TRA" :
#        jobs = [ TRA_methylation(entry,args.cpg,args.gpc,args.verbose,q) for entry in sv_type ]
        jobs = [ pool.apply_async(TRA_methylation,
            args = (entry,args.bam,args.cpg,args.gpc,args.verbose,q))
            for entry in sv_type ]
    output = [ p.get() for p in jobs ]
    q.put('kill')
    # done
    q.join()
    pool.close()
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)

if __name__=="__main__":
    main()

