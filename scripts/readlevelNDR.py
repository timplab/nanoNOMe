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
    parser.add_argument('-a','--alpha',type=float,required=False,default=-1,
            help="FDR q-value threshold")
    parser.add_argument('-t','--thr',type=int,required=False,default=3,
            help="coverage threshold for NDR candidate")
    args = parser.parse_args()
    return args

def readInput(inputpath) :
    if inputpath.split('.')[-1]=="gz":
        in_fh=gzip.open(inputpath,'rt')
    else :
        in_fh = open(inputpath,'r')
    return in_fh

testResult=namedtuple('testResult',['pval','readsopen','numreads','ndr_freq','region_freq'])
class methCounts :
    def __init__(self):
        self.methcount,self.testcount,self.totmeth,self.totcount=[0,0,0,0]
class candidateInt :
# get number of methylated and total number of calls in that region
# get the average frequency in the entire region
# perform one-tailed binomial test
    def __init__(self,intersect):
        self.region=intersect.region
        self.metharrays=[]
        self.openreads=0
        self.totreads=0
        self.center=self.getCenter()
        self.addread(intersect.methread)
        self.center=intersect.center
    def getCenter(self):
        start=self.region.start
        end=self.region.end
        return start+int((end-start)/2)
    def addread(self,read):
        self.metharrays.append(read.callarray)
        self.totreads+=1
    def isnew(self,newregion):
        if (self.region.rname == newregion.rname
                and self.region.start == newregion.start
                and self.region.end == newregion.end):
            return False
        else :
            return True
    def setTestRegion(self,upstream,downstream) :
        center=self.center#self.getCenter()
        if self.region.strand == "-":
            testreg=[center-downstream,center+upstream]
        else :
            testreg=[center-upstream,center+downstream]
        self.testRegion=testreg
    def NDRtest(self,cov) :
        counts=methCounts()
        for metharray in self.metharrays:
            metharray=metharray[metharray[:,1]!=-1]
            pos=metharray[:,0]
            meth=metharray[:,1]
            testmeth=meth[(pos>=self.testRegion[0])&
                    (pos<self.testRegion[1])]
            counts.testcount+=len(testmeth)
            counts.methcount+=sum(testmeth)
            counts.totcount+=len(meth)
            counts.totmeth+=sum(meth)
            if len(testmeth) == 0 :
                self.totreads-=1
            else : 
                if (float(sum(testmeth))/len(testmeth)
                        > float(sum(meth))/len(meth)):
                    self.openreads+=1
        if counts.testcount < cov :
            pval=-1
            totfreq=-1
            ndrfreq=-1
        else :
            pval=ss.binom_test(x=counts.methcount,
                    n=counts.testcount,
                    p=float(counts.totmeth)/counts.totcount)
            totfreq=float(counts.totmeth)/counts.totcount
            ndrfreq=float(counts.methcount)/counts.testcount
        self.result=testResult(pval,self.openreads,self.totreads,ndrfreq,totfreq)

fdrResult=namedtuple('fdrResult',['fdr','pval','readsopen','numreads','ndr_freq','region_freq'])
def getFDR(regdict,alpha):
    keys=list(sorted(regdict.keys()))
    pvals=np.asarray([regdict[i].result.pval for i in keys])
    sortind=np.argsort(pvals)
    pvals_sort=np.take(pvals,sortind)
    multiplier=len(pvals)/(np.array(range(len(pvals)))+1)
    fdr=pvals_sort*multiplier
    thr=np.argmax(fdr>alpha)
    if alpha == -1 :
        sigind=[keys[i] for i in sortind]
    else :
        sigind=[keys[i] for i in sortind[:thr]]
    regions_sorted=[regdict[i] for i in sigind]
    fdr_result=[(regions_sorted[i].region , fdrResult(fdr[i],
        regions_sorted[i].result.pval,
        regions_sorted[i].result.readsopen,
	regions_sorted[i].result.numreads,
        regions_sorted[i].result.ndr_freq,
        regions_sorted[i].result.region_freq)) for i in range(len(sigind))]
    return fdr_result

def printResult(testresult,out):
    print("\t".join([str(x) for x in testresult[0]]+[str(x) for x in testresult[1]]),file=out)

def makekey(region):
    return ".".join([str(x) for x in region])

def readlevelNDR(gpc,out,cov=3,window=10,alpha=0.05) :
    regdict=dict()
    i=0
    for line in gpc:
        intersect_read=methInt(line)
        key=makekey(intersect_read.region)
        if key not in regdict.keys():
            regdict[key]=candidateInt(intersect_read)
        else : 
            regdict[key].addread(intersect_read.methread)
        i+=1
        if i % 10000 == 0 :
            print("Read in {} reads".format(i),file=sys.stderr)
        # for debugging
#        if i == 10000 :
#            break
    print("Total {} reads read".format(i),file=sys.stderr)
    print("Performing binomial tests")
    i=0
    for key in list(regdict.keys()) :
        # NDR cnadidate is -100 to +50 bp of the TSS (stranded)
        regdict[key].setTestRegion(100,50)
        regdict[key].NDRtest(cov)
        if regdict[key].result.pval == -1 :
            regdict.pop(key)
        i+=1
        if i % 1000 == 0 :
            print("Processed {} regions".format(i),file=sys.stderr)
    print("Total {} regions processed".format(i),file=sys.stderr)
    print("Calculating FDR",file=sys.stderr)
    fdr_result=getFDR(regdict,alpha)
    print("Writing out the result",file=sys.stderr)
    [printResult(x,out) for x in fdr_result]

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
