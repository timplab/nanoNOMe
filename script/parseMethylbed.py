#! /usr/bin/env python

import math
import bisect
import sys
import argparse
import gzip
from collections import namedtuple
from methylation_parsers import MethRead

class SiteStats:
    def __init__(self,methcall):
        self.rname=methcall.rname
        self.pos=methcall.pos
        self.num_reads = 0
        self.num_methylated = 0
        self.seq=methcall.seq
    def update(self,methcall):
        if methcall.call==-1:
            return
        self.num_reads+=1
        self.num_methylated+=methcall.call
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
def printBedgraph(stats):
    print("{}\t{}\t{}\t{}\t{}".format(
        stats.rname,
        stats.pos,
        stats.pos+1,
        stats.num_methylated,
        stats.num_reads-stats.num_methylated))

def getFreq(args,in_fh):
    sites=dict()
    for line in in_fh:
        read=MethRead(line)
        sitekeys=sorted(sites.keys())
#        print(sitekeys)
        try : 
            if read.rname != sites[sitekeys[0]].rname :
                printind=len(sitekeys)
            else : 
                printind=bisect.bisect_left(sitekeys,read.keys[0])
        except IndexError : 
            printind=0
        # print everything in sites if chromosome is different
        if printind != 0 :
            for i in range(printind):
                key=sitekeys[i]
                printBedgraph(sites[key])
                sites.pop(key)
        for key in read.keys:
            if  key not in sites.keys():
                sites[key] = SiteStats(read.calldict[key])
            sites[key].update(read.calldict[key])
    for key in sorted(sites.keys()):
        printBedgraph(sites[key])

def getReadlevel(args,in_fh):
    for line in in_fh:
        read=MethRead(line)
        callkeys=sorted(read.calldict)
        for key in callkeys:
            methcall=read.calldict[key]
            print("{}\t{}\t{}\t{}\t{}".format(
                methcall.rname,
                methcall.pos,
                methcall.call,
                read.qname,
                methcall.seq))
class Feature :
    def __init__(self,bedentry):
        self.fields=bedentry.strip().split("\t")
        self.rname=self.fields[0]
        self.start=int(self.fields[1])
        self.end=int(self.fields[2])
        self.strand=self.fields[-2]
    
def getMethIntersect(args,in_fh):
    for line in in_fh:
        fields=line.strip().split("\t")
        read=MethRead("\t".join(fields[0:6]))
        feature=Feature("\t".join(fields[6:]))
        for key in read.keys :
            call=read.calldict[key]
            if key >= feature.start and key < feature.end :
                print("{}\t{}\t{}\t{}\t{}".format(
                    call.rname,
                    call.pos,
                    call.call,
                    read.qname,
                    call.pos-feature.start))

        
def main() :
    # parse args
    parser = argparse.ArgumentParser( description='Calculate methylation frequency from sorted methylation bed file')
    parser.add_argument('-i', '--input', type=str, required=False)
    parser.add_argument('-t','--type',type=str,required=False,default="frequency",
            help="type of output: one of frequency,readlevel")
    parser.add_argument('-m','--mod',type=str,required=False,default="cpg")
    args = parser.parse_args()
    if args.input:
        if args.input.split('.')[-1]=="gz":
            in_fh=gzip.open(args.input,'rt')
        else :
            in_fh = open(args.input,'r')
    else:
        in_fh = sys.stdin
    if args.type == "frequency":
        getFreq(args,in_fh)
    elif args.type == "readlevel":
        getReadlevel(args,in_fh)
    elif args.type == "intersect":
        getMethIntersect(args,in_fh)
    in_fh.close()

if __name__=="__main__":
    main()

