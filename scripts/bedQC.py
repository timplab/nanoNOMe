#! /usr/bin/env python

import math
import sys
import csv
import argparse
from collections import namedtuple
import numpy as np

def parseArgs() :
    parser = argparse.ArgumentParser( description='yield and read length qc from bed file')
    parser.add_argument('input', nargs='?',type=argparse.FileType('r'),
            default=sys.stdin,help="input bed")
    parser.add_argument('-l','--length',metavar='N',type=int,required=False,
            default=200,help='filter reads that are shorter than N')
    args = parser.parse_args()
    return args

def getN50(lengths_list) :
    lengths_sorted=np.sort(lengths_list)
    lengths_cumsum=np.cumsum(lengths_sorted)
    n50_index=np.argmax(lengths_cumsum>max(lengths_cumsum)/2)
    n50_length=lengths_sorted[n50_index]
    return n50_length

QCnumbers = namedtuple('QCnumbers', ['read_num','total_bp','N50_length','Q1_length','mean_length','Q3_length','max_length'])
class bedSummary :
    def __init__(self,cutoff_length=200) :
        self.thr=cutoff_length
        self.read_lengths = []
        self.read_num=0
    def update(self,fields) :
        rlen = int(fields[2])-int(fields[1])
        if rlen < self.thr :
            return
        self.read_num+=1
        self.read_lengths.append(rlen)
    def getnumbers(self) :
        mean_length=int(np.mean(self.read_lengths))
        total=sum(self.read_lengths)
        q1=int(np.percentile(self.read_lengths,25))
        q3=int(np.percentile(self.read_lengths,75))
        max_length=max(self.read_lengths)
        n50=getN50(self.read_lengths)
        qcnumbers=QCnumbers(self.read_num,total,
                n50,q1,mean_length,q3,max_length)
        return qcnumbers

def bedQC(in_fh,cutoff_length) : 
    summary=bedSummary(cutoff_length)
    for line in in_fh :
        fields = line.strip().split('\t')
        summary.update(fields)
    qcnumbers=summary.getnumbers()
    print('\t'.join([str(x) for x in qcnumbers]))


if __name__=="__main__":
    args=parseArgs()
    bedQC(args.input,args.length)
