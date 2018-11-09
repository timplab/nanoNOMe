#! /usr/bin/env python
import sys
import argparse

def parseArgs():
    parser = argparse.ArgumentParser( description='Generate methylation bedGraph file')
    parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.5,
            help="absolute value of threshold for methylation call (default : 2.5)")
    parser.add_argument('-i', '--input', type=str, required=False,help="input methylation tsv file (default stdin)")
    parser.add_argument('-m', '--mod',type=str,required=False,default='cpg',help="modification motif; one of cpg,gpc")
    parser.add_argument('-e', '--exclude',type=str,required=False,help="motif to exclude from reporting")
    parser.add_argument('-w', '--window',type=int,required=False,default=1,
            help="number of nucleotides to report on either side of called nucleotide")
    parser.add_argument('--offset',type=int,required=False,default=1,
            help="nanopolish coordinate offset (1-based)")
    args = parser.parse_args()
    assert(args.call_threshold is not None)
    return args

class MethQuery:
    def __init__(self,string,mod,offset):
        self.fields=string.strip().split("\t")
        self.rname=self.fields[0]
        self.mod=mod
        if mod=="cpg" :
            self.start=int(self.fields[2])+offset
            self.motif="CG"
        elif mod=="gpc" :
            self.start=int(self.fields[2])+1+offset # offset from methylation calls - 1-based
            self.motif="GC"
        self.strand=self.fields[1]
        self.qname=self.fields[4]
        self.ratio=float(self.fields[5])
        self.unmeth=float(self.fields[7])
        self.recsites=int(self.fields[9])
        self.seq=self.fields[10][3:]
    def methcall(self,thr):
        if abs(self.ratio) < thr :
            call="x"
        elif self.ratio > 0 :
            call="m"
        else : 
            call="u"
        return call
    def parse(self,motif,window):
        pos=[]
        seqs=[]
        newpos=self.seq.find(motif)
        while newpos != -1:
            seq=self.seq[newpos-window:newpos+window+2]
            pos.append(newpos)
            seqs.append(seq)
            newpos=self.seq.find(motif,newpos+1)
        coord=[x-y for x,y in zip(pos,[5]+pos)]
        return coord,seqs

class readQuery:
    def __init__(self,query,args):
        self.qname=query.qname
        self.rname=query.rname
        self.start=query.start
        self.end=query.start
        self.strand=query.strand
        self.dist=[]
        self.seq=[]
        self.ratio=[]
        self.unmeth=[]
        self.call=[]
        self.mod=args.mod
        self.motif=query.motif
        self.thr=args.call_threshold
        self.window=args.window
        self.update(query)

    def update(self,query):
        dists,seqs=query.parse(self.motif,self.window)
        dists[0]=query.start-self.end
        self.dist=self.dist+dists
        self.end=self.start+sum(self.dist)
        self.seq=self.seq+seqs
        call=query.methcall(self.thr)
        for i in range(query.recsites):
            self.ratio.append(query.ratio)
            self.unmeth.append(query.unmeth) 
            self.call.append(call)

    def summarizeCall(self):
        summary=""
        for ind,call in zip(self.dist,self.call):
            summary=summary+"{}{}".format(ind,call)
        return summary
    def concatenateCall(self):
        summary=""
        for ind,call in zip(self.dist,self.call) :
            summary=summary+(ind-1)*"."+call
        return summary

def catList(strlist,delim=""):
    return delim.join([str(x) for x in strlist])

def PrintRead(read):
    print("\t".join([read.rname,str(read.start-1),str(read.end),
        read.qname,
        read.summarizeCall(),
        catList(read.ratio,","),
        catList(read.seq,",")]))

def summarizeMeth(args):
    if args.input:
        in_fh = open(args.input)
    else:
        in_fh = sys.stdin
    for line in in_fh:
        try : 
            # skip headers (in concatenated tsvs) and otherwise faulty entries
            query=MethQuery(line,args.mod,args.offset)
        except :
            continue
        # skip queries that have a motif in the call group
        if args.exclude:
            if args.exclude in query.seq:
                continue
        # update query
        try :
            if ((query.qname != read.qname) or
                    (query.rname != read.rname) or
                    (query.start < read.end)):
                PrintRead(read)
                read=readQuery(query,args)
            else : 
                read.update(query)
        except :
            read=readQuery(query,args)
    try : 
        PrintRead(read)
    except :
        pass
    in_fh.close()

if __name__=="__main__":
    args=parseArgs()
    summarizeMeth(args)
