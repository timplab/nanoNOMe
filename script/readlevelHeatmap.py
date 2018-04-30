#!/home/isac/.conda/envs/isacenv/bin/python
# readlevel analysis of NOMe-seq : plot heatmap of co-methylation
import scipy.stats as ss
import sys
import argparse
import gzip
from collections import namedtuple
from nomeseq_parsers import methInt,MethRead
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
blues=["#6BAED6","#4292C6","#2171B5","#08519C","#08306B","#041938"]


def parseArgs():
    parser = argparse.ArgumentParser( description='plot heatmap of methylation co-occurrence')
    parser.add_argument('-g', '--gpc', type=str, required=True,
            help="read-level GpC methylation bed file")
    parser.add_argument('-o', '--out', type=str, required=False,
            help="output file path - defaults out to stdout")
    parser.add_argument('-w','--window',type=int,required=False,default=10,
            help="binning window for calculating methylation")
    parser.add_argument('-t','--threshold',type=int,required=False,default=5,
            help="upper limit on number to be displayed as the darkest pixel")
    args = parser.parse_args()
    try : 
        args.out
    except NameError :
        args.out = "heatmap.pdf"
    return args

def readInput(inputpath) :
    if inputpath.split('.')[-2]=="gz":
        in_fh=gzip.open(inputpath,'rt')
    else :
        in_fh = open(inputpath,'r')
    return in_fh

class HeatmapRegion :
    def __init__(self,intersect):
        self.region=intersect.region
        self.metharrays=[]
        self.totreads=0
        self.addread(intersect.methread)
        self.center=intersect.center
    def initDict(self,win) :
        start=round(self.region.start/win)*win
        end=round(self.region.end/win)*win
        self.keys=list(range(start,end,win))
        self.matrix=np.zeros((len(self.keys),len(self.keys)))
        self.centerkey=round(self.center/win)*win
    def addread(self,read):
        self.metharrays.append(read.callarray)
        self.totreads+=1
    def getMatrix(self,win) :
        self.initDict(win)
        for metharray in self.metharrays:
            calls=[]
            loci=[]
            for pos,call in metharray :
                key=round(pos/win)*win
                loci.append(key)
                if call == 1 :
                    calls.append(key)
            calls=list({ self.keys.index(x) for x in calls if x in self.keys })
            loci=list({ self.keys.index(x) for x in loci if x in self.keys })
            for i in range(len(calls)) :
                for j in calls :
                    self.matrix[calls[i]][j]+=1
            # fill in the diagonal with call loci - for determining which loci had any call
            for i in range(len(loci)) :
                self.matrix[loci[i]][loci[i]]=self.totreads/4
        return self.matrix
    def plot(self,thr,out) :
        matrix_normalized=self.matrix/self.totreads
        df=pd.DataFrame(matrix_normalized,index=self.keys,columns=self.keys)
        plt.figure()
        tickind=round(len(self.keys)/4)
        ax=sns.heatmap(df,
                xticklabels=self.keys[0::tickind],
                yticklabels=self.keys[0::tickind],
                cmap=sns.color_palette(blues),
                vmin=0,
                mask=matrix_normalized<0.25)
        ax.set(xlabel='Coordinate on {}'.format(self.region.rname))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        centeridx=self.keys.index(self.centerkey)
        plt.plot([centeridx,centeridx],[0,len(self.keys)])
        plt.plot([0,len(self.keys)],[centeridx,centeridx])
        plt.show()
        try :
            name=self.region.id
        except NameError :
            name=self.region.name
        title="{} {}:{}-{}".format(self.region.id,
                self.region.rname,
                self.region.start,
                self.region.end)
        plt.title(title)
        out.savefig()
        plt.close()


def makekey(region):
    return ".".join([str(x) for x in region])

def readlevelHeatmap(gpc,out,thr=5,window=10) :
    regdict=dict()
    i=0
    for line in gpc:
        intersect_read=methInt(line)
        key=makekey(intersect_read.region)
        if key not in regdict.keys():
            regdict[key]=HeatmapRegion(intersect_read)
        else : 
            regdict[key].addread(intersect_read.methread)
        i+=1
        if i % 1000 == 0 :
            print("Read in {} reads".format(i),file=sys.stderr)
    print("Total {} reads read, plotting heatmaps".format(i),file=sys.stderr)
    i=0
    with PdfPages(out) as pdf :
        for key in list(regdict.keys()) :
            regdict[key].getMatrix(window)
            regdict[key].plot(thr,pdf)
            i+=1
            print("Plotted {}".format(regdict[key].region.name),
                    file=sys.stderr)
            # for debugging
            if i == 10:
                break
    print("Total {} regions processed".format(i),file=sys.stderr)

if __name__=="__main__":
    args=parseArgs()
    gpc=readInput(args.gpc)
    readlevelHeatmap(gpc,args.out,args.threshold,args.window)
    gpc.close()
