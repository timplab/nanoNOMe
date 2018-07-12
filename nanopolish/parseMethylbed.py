#! /usr/bin/env python

import os
import math
import bisect
import sys
import argparse
import gzip
from collections import namedtuple
from methylbed_utils import MethRead

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # get module
    module=sys.argv[1]
    # parse args
    parser = argparse.ArgumentParser( description='parse methylation bed files' )
    parser.add_argument('-i', '--input', type=str, required=False)
    parser.add_argument('-o', '--out', type=str, required=False,
            help="output file path - defaults out to stdout")
    parser.add_argument('-m','--mod',type=str,required=False,default="cpg")
    args = parser.parse_args(sys.argv[2:])
    args.srcdir=srcdir
    args.module=module
    if args.mod=="cpg":
        args.motif="CG"
    elif args.mod=="gpc":
        args.motif="GC"
    return args

class SiteStats:
    def __init__(self,methcall,chrom):
        self.rname=chrom
        self.pos=methcall.pos
        self.num_reads = 0
        self.num_methylated = 0
        self.seq=methcall.seq
    def update(self,methcall):
        if methcall.call==-1:
            return
        self.num_reads+=1
        self.num_methylated+=methcall.call
    def printFreq(self,context):
        motifidx=self.seq.index(context)
        motifcontext=self.seq[motifidx:motifidx+2]
        if context == "CG" :
            tricontext=self.seq[motifidx-1:motifidx+2]
        elif context == "GC" :
            tricontext=self.seq[motifidx:motifidx+3]
        print("\t".join([str(x) for x in [
            self.rname,
            self.pos+1,
            "*",
            self.num_methylated,
            self.num_reads-self.num_methylated,
            motifcontext,tricontext]]))
def printsite(stats):
    print("{}\t{}\t{}\t{}".format(stats.rname,stats.pos,
        stats.num_reads,stats.num_methylated))
def printcyto(stats,mod):
    if mod=="cpg":
        motif="CG"
    elif mod =="gpc" :
        motif="GC"
    conind=stats.seq.index(motif)
    context=stats.seq[conind:conind+4]
    print("{}\t{}\t*\t{}\t{}\tCG\t{}".format(stats.rname,stats.pos,
        stats.num_methylated,
        stats.num_reads-stats.num_methylated,
        context))

def getFreq(args,in_fh,out):
    sites=dict()
    for line in in_fh:
        try : 
            line = line.decode('ascii')
        except AttributeError :
            pass
        read=MethRead(line)
        sitekeys=sorted(sites.keys())
#        print(sitekeys)
        try : 
        # print everything in sites if chromosome is different
            if read.rname != sites[sitekeys[0]].rname :
                printind=len(sitekeys)
            else : 
                printind=bisect.bisect_left(sitekeys,read.keys[0])
        except IndexError : 
            printind=0
        if printind != 0 :
            for i in range(printind):
                key=sitekeys[i]
                sites[key].printFreq(args.motif)
                sites.pop(key)
        for key in read.keys:
            if  key not in sites.keys():
                sites[key] = SiteStats(read.calldict[key],read.rname)
            sites[key].update(read.calldict[key])
    for key in sorted(sites.keys()):
        sites[key].printFreq(args.motif)

def getReadlevel(args,in_fh,out):
    for line in in_fh:
        try : 
            line = line.decode('ascii')
        except AttributeError :
            pass
        read=MethRead(line)
        callkeys=sorted(read.calldict)
        for key in callkeys:
            methcall=read.calldict[key]
            outlist=[read.rname]+[str(x) for x in methcall]
            outlist.insert(3,read.qname)
            print("\t".join(outlist),file=out)

def getMethIntersect(args,in_fh,out) :
    for line in in_fh : 
        try : 
            line = line.decode('ascii')
        except AttributeError :
            pass
        intersect=methInt(line)
        read=intersect.methread
        region=intersect.region
        start=region.start
        end=region.end
        for key in read.keys :
            methcall=read.calldict[key]
            if key >= start and key < end :
                outlist=([region.rname]+
                        [str(x) for x in methcall]+
                        [str(x) for x in region.items]+
                        [str(methcall.pos-region.start)])
                outlist.insert(2,str(methcall.pos+1))
                outlist.insert(4,read.qname)
                print("\t".join(outlist),file=out)
       
def main() :
    args=parseArgs()
    # read in input
    if args.input:
        if args.input.split('.')[-1]=="gz":
            in_fh=gzip.open(args.input,'rt')
        else :
            in_fh = open(args.input,'r')
    else:
        in_fh = sys.stdin
    # read in output
    if args.out:
        out=open(args.out,'w')
    else :
        out=sys.stdout
    #
    if args.module == "frequency":
        getFreq(args,in_fh,out)
    elif args.module == "readlevel":
        getReadlevel(args,in_fh,out)
    elif args.module == "intersect" :
        getMethIntersect(args,in_fh,out)

    in_fh.close()

if __name__=="__main__":
    main()

