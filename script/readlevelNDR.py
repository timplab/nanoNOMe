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
    parser.add_argument('-g', '--gpc', type=str, required=True,
            help="read-level GpC methylation bed file")
    parser.add_argument('-o', '--out', type=str, required=False,
            help="output file path - defaults out to stdout")
    parser.add_argument('-w','--window',type=int,required=False,default=10,
            help="binning window for calculating methylation")
    parser.add_argument('-a','--alpha',type=float,required=False,default=0.05,
            help="FDR q-value threshold")
    parser.add_argument('-t','--thr',type=int,required=False,default=1,
            help="coverage threshold for NDR candidate")
    args = parser.parse_args()
    return args

def readInput(inputpath) :
    if inputpath.split('.')[-1]=="gz":
        in_fh=gzip.open(inputpath,'rt')
    else :
        in_fh = open(inputpath,'r')
    return in_fh

testResult=namedtuple('testResult',['rname','start','end','qname','nuc','NDRfreq'])
def printResult(testresult,out):
    print("\t".join([str(x) for x in testresult]),file=out)

class candidateInt :
# get number of methylated and total number of calls in that region
# get the average frequency in the entire region
# perform one-tailed binomial test
    def __init__(self,intersect):
        self.region=intersect.region
        self.methread=intersect.methread
        self.center=self.getCenter()
    def getCenter(self):
        start=self.region.start
        end=self.region.end
        return start+int((end-start)/2)
    def setTestRegion(self,upstream,downstream) :
        center=self.getCenter()
        if self.region.strand == "-":
            testreg=[center-downstream,center+upstream]
        else :
            testreg=[center-upstream,center+downstream]
        self.testRegion=testreg
    def NDRtest(self,cov) :
        metharray=self.methread.callarray
        metharray=metharray[metharray[:,1]!=-1]
        pos=metharray[:,0]
        meth=metharray[:,1]
        testmeth=meth[(pos>=self.testRegion[0])&
                (pos<self.testRegion[1])]
        testcnt=len(testmeth)
        methcnt=sum(testmeth)
        if testcnt==0 :
            methfreq=-1
            nuc=-1
        else : 
            methfreq=float(methcnt)/testcnt
            nuc=int(methfreq>0.5)
        self.result=testResult(self.methread.rname,
                self.methread.start,
                self.methread.end,
                self.methread.qname,
                nuc,methfreq)

def makekey(strlist):
    return ".".join(strlist)

def readlevelNDR(gpc,out,cov=3,window=10,alpha=0.05) :
    regdict=dict()
    i=0
    for line in gpc:
        i+=1
        cand=candidateInt(methInt(line))
        # NDR cnadidate is -100 to +50 bp of the TSS (stranded)
        cand.setTestRegion(100,50)
        cand.NDRtest(cov)
        key=makekey([cand.region.rname]+[str(x) for x in cand.testRegion])
        try : 
            regdict[key].append(cand.result)
        except KeyError :
            regdict[key]=[cand.result]
        # for debugging
#        if i==10000 : 
#            break
    reglist=[]
    keys=sorted(regdict.keys())
    for key in keys :
        calls=[x.nuc for x in regdict[key] if x.nuc != -1]
        if len(calls)==0:
            opencount=-1
            closecount=-1
            freq=-1
        else : 
            totcount=len(calls)
            opencount=sum(calls)
            closecount=totcount-opencount
            freq=float(opencount)/totcount
        reglist.append([opencount,closecount,freq])
    sortindex=np.argsort(np.array(reglist)[:,2])[::-1]
    [ print("\t".join(keys[i].split(".")+
        [str(x) for x in reglist[i]]),file=out) for i in sortindex if reglist[i][2] != -1]

if __name__=="__main__":
    args=parseArgs()
    gpc=readInput(args.gpc)
    if args.out:
        out=open(args.out,'w')
    else :
        out=sys.stdout
    readlevelNDR(gpc,out,args.thr,args.window,args.alpha)
    gpc.close()
    if args.out:
        out.close()
