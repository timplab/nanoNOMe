#! /usr/bin/env python
import sys
import os
import argparse
import pysam
import numpy as np
from collections import namedtuple
from methylbed_utils import SnifflesEntry,make_coord,read_bam
import multiprocessing as mp
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='dealing with sniffles vcf ')
    subparsers = parser.add_subparsers()
    # parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-v', '--verbose', action='store_true',default=False,
            help="verbose output")
    parent_parser.add_argument('-t','--threads',type=int,default=1,
            help="threads (default 1)")
    parent_parser.add_argument('-s', '--sniffles', nargs='?', type=argparse.FileType('r'),
            required=False, default=sys.stdin,help="input vcf file path (or stdin)")
    parent_parser.add_argument('-o', '--output', nargs='?', type=str,
            required=False, default="stdout",help="output file path (or stdout)")
    # parser to print predicted bed coordinate ranges of SVs
    parser_bed = subparsers.add_parser('bed',parents=[parent_parser],
            help = "make predictions of bed coordinates of breakpoints")
    parser_bed.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help = "bam file")
    parser_bed.add_argument('-w','--window',type=int,required=False,
            default=200,help="window for flanking regions")
    parser_bed.set_defaults(what="bed",func=snifflesBed)
    # parser to subset bam entries 
    parser_bam = subparsers.add_parser('bam',parents=[parent_parser],
            help = "subset bam entries that fall in SV bps")
    parser_bam.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help = "bam file to subset")
    parser_bam.set_defaults(what="bam",func=snifflesBam)
    args = parser.parse_args()
#    print(args)
    return args

# https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def listener(q,out,inbam=None,what="bam",verbose=False) :
    '''listens for messages on the q, writes to file. '''
    if verbose : print("writing output to {}".format(out),file=sys.stderr)
    if what == "bam" :
        in_fh = pysam.AlignmentFile(inbam,'rb')
        if out == "stdout" :
            out_fh = sys.stdout
            print(in_fh.header,file=out_fh)
            def printline(m,out_fh) :
                print(m,file=out_fh) 
        else :
            out_fh = pysam.AlignmentFile(out,'wb',template=in_fh)
            def printline(m,out_fh) :
                read = pysam.AlignedSegment.fromstring(m,out_fh.header)
                out_fh.write(read)
        in_fh.close()
    else :
        if out == "stdout" :
            out_fh = sys.stdout
        else :
            out_fh = open(out,'w')
        def printline(m,out_fh) :
            print(m,file=out_fh)
    key_list = list()
    while 1:
        m = q.get()
        if m == 'kill':
            break
        key,line = m
        if key not in key_list:
            key_list.append(key)
            printline(line,out_fh)
        else :
            if verbose: print("{} already written - skipping".format(key),file=sys.stderr)
        q.task_done()
    if verbose : print("total {} lines written".format(len(key_list)),file=sys.stderr)
    out_fh.close()
    q.task_done()

def printBed(outlist,q) :
    out_line = '\t'.join(outlist)
    if int(outlist[1])<0 : outlist[1] = str(0)
    q.put((out_line,out_line))

def snifflesBed(bamfn,win,sv,verbose,q) : 
    sv.activate()
    if sv.allele == "0/0" : return
    bamwin = 500
    if sv.type == "TRA" :
        coord = make_coord(sv.info["CHR2"],sv.info["END"]-bamwin,sv.info["END"]+bamwin)
        # fetch reads to determine the right edge
        bam_dicts = read_bam(bamfn,coord)
        edge_list = list()
        for qname in sv.rnames :
            try : 
                bam = bam_dicts[qname][-1]
            except : continue
            edge_list.append(bam.reference_end)
        if len(edge_list) == 0 : return
#        edge = max(edge_list)
        edge = int(np.median(edge_list))
        if ( max(edge_list)-edge > win and 
                edge-min(edge_list) > win or 
                edge < sv.info["END"] ): return 
        printBed([sv.chrom,
            str(sv.pos-1),
            str(sv.pos),
            sv.id,".",".",sv.type,"target","breakpoint"],q)
        printBed([sv.chrom,
            str(sv.pos-1-win),
            str(sv.pos),
            sv.id,".",".",sv.type,"target","upstream"],q)
        printBed([sv.chrom,
            str(sv.pos-1),
            str(sv.pos+win),
            sv.id,".",".",sv.type,"target","downstream"],q)
        printBed([sv.info["CHR2"],
            str(sv.info["END"]-1),
            str(edge),sv.id,".",".",sv.type,"insert","region"],q)
        printBed([sv.info["CHR2"],
            str(sv.info["END"]-1-win),
            str(sv.info["END"]),sv.id,".",".",sv.type,"insert","upstream"],q)
        printBed([sv.info["CHR2"],
            str(edge-1),
            str(edge+win),sv.id,".",".",sv.type,"insert","downstream"],q)
        printBed([sv.info["CHR2"],
            str(sv.info["END"]-1),
            str(sv.info["END"]+win),sv.id,".",".",sv.type,"insert","five-prime"],q)
        printBed([sv.info["CHR2"],
            str(edge-1-win),
            str(edge),sv.id,".",".",sv.type,"insert","three-prime"],q)

def snifflesBam(bamfn,window,sv,verbose,q) :
    sv.activate()
    win = 500
    sv = SnifflesEntry(svline)
    if sv.allele == "0/0" : return
    start = sv.pos-win
    end = sv.pos+win
    coord = make_coord(sv.chrom,start,end)
    with pysam.AlignmentFile(bamfn,"rb") as bam :
        bam_entries = [ x for x in bam.fetch(region=coord) ]
    for bam in bam_entries : 
        if bam.query_name in sv.rnamelist : 
            key = '.'.join([bam.query_name,bam.reference_name,str(bam.reference_start)])
            q.put((key,bam.to_string()))

if __name__=="__main__":
    args=parseArgs()
    if args.verbose : print("reading in SVs",file=sys.stderr)
    svlines = [SnifflesEntry(x) for x in args.sniffles.readlines() if x[0]!="#"]
    # currently working only with TRA
    svlines = [x for x in svlines if x.type == "TRA" ]
    # initialze mp
    if args.verbose : print("using {} parallel processes".format(args.threads),file=sys.stderr)
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    # watcher for output
    watcher = pool.apply_async(listener,(q,args.output,args.bam,args.what,args.verbose))
    jobs = [ pool.apply_async(args.func,
        args = (args.bam,args.window,sv,args.verbose,q))
        for sv in svlines ]
    output = [ p.get() for p in jobs ]
    # done
    q.put('kill')
    q.join()
    pool.close()
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)
