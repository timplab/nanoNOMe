#! /usr/bin/env python
# script for indexing reads to paths in tarballs - currently needs the fast5list.txt to work
import sys
import os
import argparse
import tarfile
from collections import namedtuple

def parseArgs():
    parser = argparse.ArgumentParser( description='dealing with tar index for fast5s')
    parser.add_argument('-v', '--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i', '--input', nargs='?', type=argparse.FileType('r'),
            required=False, default=sys.stdin,help="input vcf file path (or stdin)")
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('w'),
            required=False, default=sys.stdout,help="output file path (or stdout)")
    parser.add_argument('what',nargs='+',help="element(s) to extract")
    args = parser.parse_args()
    return args

def parsewhats(whats) :
    outwhats=list()
    if "coords" in whats :
        outwhats.append("chrom")
        outwhats.append("start")
        outwhats.append("end")
    if "svtype" in whats :
        outwhats.append("svtype")
    return outwhats

def findelement(s,what,delimiter) :
    if delimiter == ";" :
        firststring=what.upper()+"=" 
        endstring=";"
    try:
        start = s.rindex( firststring ) + len( firststring )
        end = s.index( endstring, start )
        return s[start:end]
    except ValueError:
        return ""


def getwhats(line,whats):
    out=list()
    fields=line.strip().split("\t")
    if "chrom" in whats :
        if "chr" in fields[0] :
            chrom=fields[0]
        else : chrom="chr"+fields[0]
        out.append(chrom)
    if "start" in whats :
        out.append(fields[1])
    if "end" in whats :
        out.append(findelement(fields[7],"end",";"))
    if "svtype" in whats :
        out.append(fields[4])
    return out

def printBed(outlist,out):
    outlist.insert(4,".")
    outlist.insert(5,".")
    print("\t".join(outlist),file=out)

if __name__=="__main__":
    args=parseArgs()
#    print(args)
#    sys.exit()
    whats=parsewhats(args.what)
    for line in args.input :
        for what in whats :
            outlist=getwhats(line,whats)
            printBed(outlist,out=args.output)
                
