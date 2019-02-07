#!/home/isac/.conda/envs/isacenv/bin/python
# readlevel analysis of NOMe-seq
import scipy.stats as ss
import sys
import argparse
import gzip
from collections import namedtuple
from nomeseq_parsers import methInt
import numpy as np

def parseArgs():
    parser = argparse.ArgumentParser( description='detect NDRs by interrogating read-level')
    parser.add_argument('-c', '--cpg', type=str, required=True,
            help="output of bedtools intersect of region with read-level GpC methylation bed file")
    parser.add_argument('-o', '--out', type=str, required=False,
            help="output file path - defaults out to stdout")
    parser.add_argument('-w','--window',type=int,required=False,default=10,
            help="binning window for calculating methylation")
    parser.add_argument('-a','--alpha',type=float,required=False,default=-1,
            help="FDR q-value threshold")
    parser.add_argument('-t','--thr',type=int,required=False,default=3,
            help="coverage threshold for methylation calculation")
    args = parser.parse_args()
    return args

def readInput(inputpath) :
    if inputpath.split('.')[-1]=="gz":
        in_fh=gzip.open(inputpath,'rt')
    else :
        in_fh = open(inputpath,'r')
    return in_fh

testResult=namedtuple('testResult',['readsmeth','readsunmeth','region_freq','total_freq'])
class methCounts :
    def __init__(self):
        self.methcount,self.testcount,self.totmeth,self.totcount=[0,0,0,0]
class methregion:
# get number of methylated and total number of calls in that region
# get the average frequency in the entire region
    def __init__(self,intersect):
        self.region=intersect.region
        self.testreg=[self.region.start,self.region.end]#self.getTestRegion()
        self.metharrays=[]
        self.methreads=0
        self.totreads=0
        self.addread(intersect.methread)
    def getTestRegion(self): 
        if self.region.strand == "+":
            return [self.region.start-2000,self.region.end]
        elif self.region.strand == "-" :
            return [self.region.start,self.region.end+2000]
    def addread(self,read):
        self.metharrays.append(read.callarray)
        self.totreads+=1
    def getMeth(self,cov) :
        counts=methCounts()
        for metharray in self.metharrays:
            metharray=metharray[metharray[:,1]!=-1]
            pos=metharray[:,0]
            meth=metharray[:,1]
            testmeth=meth[(pos>=self.testreg[0])&
                    (pos<self.testreg[1])]
            counts.testcount+=len(testmeth)
            counts.methcount+=sum(testmeth)
            counts.totcount+=len(meth)
            counts.totmeth+=sum(meth)
            if len(testmeth) == 0 :
                self.totreads-=1
            else : 
                if (float(sum(testmeth))/len(testmeth)
                        > float(sum(meth))/len(meth)):
                    self.methreads+=1
        if counts.testcount < cov :
            totfreq=-1
            regfreq=-1
        else :
            totfreq=float(counts.totmeth)/counts.totcount
            regfreq=float(counts.methcount)/counts.testcount
        self.result=testResult(self.methreads,self.totreads-self.methreads,regfreq,totfreq)

def printResult(methreg,out):
    regionout="\t".join([str(x) for x in methreg.region])
    resultout="\t".join([str(x) for x in methreg.result])
    print("\t".join([regionout,resultout]),file=out)

def makekey(region):
    return ".".join([str(x) for x in region])

def readlevelCpGmeth(indata,out,cov=3,window=10,alpha=0.05) :
    regdict=dict()
    i=0
    for line in indata:
        intersect_read=methInt(line)
        key=makekey(intersect_read.region)
        if key not in regdict.keys():
            regdict[key]=methregion(intersect_read)
        else : 
            regdict[key].addread(intersect_read.methread)
        i+=1
        if i % 10000 == 0 :
            print("Read in {} reads".format(i),file=sys.stderr)
        # for debugging
#        if i == 10000 :
#            break
    print("Total {} reads read".format(i),file=sys.stderr)
    print("Calculating Methylation")
    i=0
    for key in list(regdict.keys()) :
        # NDR cnadidate is -100 to +50 bp of the TSS (stranded)
        regdict[key].getMeth(cov)
        if regdict[key].result.total_freq == -1 :
            regdict.pop(key)
        i+=1
        if i % 1000 == 0 :
            print("Processed {} regions".format(i),file=sys.stderr)
    print("Total {} regions processed".format(i),file=sys.stderr)
    print("Writing out the result",file=sys.stderr)
    [printResult(regdict[i],out) for i in regdict]
    

if __name__=="__main__":
    args=parseArgs()
    indata=readInput(args.cpg)
    if args.out:
        out=open(args.out,'w')
    else :
        out=sys.stdout
    readlevelCpGmeth(indata,out,args.thr,args.window,args.alpha)
    indata.close()
    if args.out:
        out.close()
