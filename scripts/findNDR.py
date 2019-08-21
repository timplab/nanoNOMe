#!/home/isac/.conda/envs/isacenv/bin/python
# this is the script to detect NDRs from NOMe-seq data - it is a simple approach from the original NOMe-seq paper, and I probably want to implement a more robust approach such as BPRMeth from scNMT-seq paper in the future
import scipy.stats as ss
import sys
import argparse
import gzip
from collections import namedtuple
from nomeseq_parsers import distFreq
import numpy as np

def parseArgs():
    parser = argparse.ArgumentParser( description='detect nucleosome depleted regions')
    parser.add_argument('-g', '--gpc', type=str, required=True,
            help="parsed distance to motif file of gpc methylation frequency")
    parser.add_argument('-c', '--cpg', type=str, required=False,
            help="parsed distance to motif file of cpg methylation frequency")
    parser.add_argument('-o', '--out', type=str, required=False,
            help="output file path - defaults out to stdout")
    parser.add_argument('-a','--alpha',type=float,required=False,default=-1,
            help="alpha for thresholding FDR")
    parser.add_argument('-t','--thr',type=int,required=False,default=10,
            help="coverage threshold for NDR candidate")
    args = parser.parse_args()
    return args

def readInput(inputpath) :
    if inputpath.split('.')[-1]=="gz":
        in_fh=gzip.open(inputpath,'rt')
    else :
        in_fh = open(inputpath,'r')
    return in_fh

testResult=namedtuple('testResult',['pval','NDRfreq','totfreq'])

class NDRcandidate :
# for each TSS, perform the following
# NDR cnadidate is -100 to +50 bp of the TSS (stranded)
# get number of methylated and total number of calls in that region
# get the average frequency in the entire region
# perform one-tailed binomial test : scipy.stats.binom_test
    def __init__(self,site):
        self.region=site.region
        self.mfreq=site.mfreq
        self.freqdict=dict()
    def addline(self,site):
        self.freqdict[site.distance]=site
    def isnew(self,site):
        if (self.region.rname == site.region.rname 
                and self.region.start == site.region.start 
                and self.region.end == site.region.end ):
            return False
        else :
            return True
    def NDRtest(self,cov) :
        posdict=self.freqdict
        dists=posdict.keys()
        reg=[x for x in dists if x>=-100 and x<= 50]
        mfreq_region=[posdict[i].mfreq for i in reg]
        methcnt=sum([x.methylated for x in mfreq_region])
        unmethcnt=sum([x.unmethylated for x in mfreq_region])
        totcount=methcnt+unmethcnt
        if totcount<cov :
            pval=-1
            regfreq=-1
            totfreq=-1
        else : 
            self.ndrstart=min([x.start for x in mfreq_region])
            self.ndrend=max([x.end for x in mfreq_region])
            regfreq=float(methcnt)/(methcnt+unmethcnt)
            totmeth=sum([x.mfreq.methylated for x in posdict.values()])
            totunmeth=sum([x.mfreq.unmethylated for x in posdict.values()])
            totfreq=totmeth/(totmeth+totunmeth)
            pval=ss.binom_test(x=[methcnt,unmethcnt],p=totfreq,alternative='greater')
        self.result=testResult(pval,regfreq,totfreq)

fdrResult=namedtuple('fdrResult',['chrom','start','end','strand','feature','fdr','pval','NDRfreq','totfreq'])
def printResult(testresult,out):
    print("\t".join([str(x) for x in testresult]),file=out)
def getFDR(regions,alpha):
    pvals=np.asarray([x.result.pval for x in regions])
    sortind=np.argsort(pvals)
    pvals_sort=np.take(pvals,sortind)
    multiplier=len(pvals)/(np.array(range(len(pvals)))+1)
    fdr=pvals_sort*multiplier
    thr=np.argmax(fdr>alpha)
    if alpha == -1 :
        sigind=sortind
    else :
        sigind=sortind[:thr]
    regions_sorted=[regions[i] for i in sigind]
    fdr_result=[fdrResult(regions_sorted[i].region.rname,
        regions_sorted[i].ndrstart,
        regions_sorted[i].ndrend,
        regions_sorted[i].region.strand,
        regions_sorted[i].region.name,
        fdr[i],
        regions_sorted[i].result.pval,
        regions_sorted[i].result.NDRfreq,
        regions_sorted[i].result.totfreq) for i in range(len(sigind))]
    return fdr_result

def findNDR(gpc,out,cov,alpha) :
    reglist=[]
    i=0
    for line in gpc:
        i+=1
        sitefreq=distFreq(line)
        try :
            candidate
        except NameError :
            candidate=NDRcandidate(sitefreq)
        if candidate.isnew(sitefreq)==False :
            candidate.addline(sitefreq)
        else : 
            candidate.NDRtest(cov)
            if candidate.result.pval!=-1 :
                reglist.append(candidate)
            candidate=NDRcandidate(sitefreq)
        if i % 100000 == 0 : 
            print("Processed {} lines".format(i),file=sys.stderr)
    print("Calculating FDR",file=sys.stderr)
    fdr_result=getFDR(reglist,alpha)
    print("Writing out the result",file=sys.stderr)
    [printResult(x,out) for x in fdr_result]

if __name__=="__main__":
    args=parseArgs()
    gpc=readInput(args.gpc)
#    cpg=readInput(args.cpg)
    if args.out:
        out=open(args.out,'w')
    else :
        out=sys.stdout
    findNDR(gpc,out,args.thr,args.alpha)
    gpc.close()
    if args.out:
        out.close()
