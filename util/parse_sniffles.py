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
    parser.add_argument('what',help="element to extract")
    args = parser.parse_args()
    return args

class SnifflesEntry :
    def __init__(self,line) :
        self.line=line.strip()
        self.fields=self.line.split("\t")
        (self.chrom,self.pos,self.id,self.ref, 
                self.type,self.qual,self.filter,self.infostring,
                self.format,self.genotype) = self.fields
        self.pos = int(self.pos)
        self.parseinfo()
        self.parsegenotype()
    def parseinfo(self) :
        self.infofields = [ x.split("=") for x in self.infostring.strip().split(";")]
        self.info = dict()
        self.info["CONFIDENCE"] = self.infofields[0][0]
        for entry in self.infofields[1:] :
            self.info[entry[0]] = entry[1]
        self.info["END"] = int(self.info["END"])
    def parsegenotype(self) :
        self.allele = self.genotype.split(":")[0]
        self.num_against = int(self.genotype.split(":")[1])
        self.num_for = int(self.genotype.split(":")[2])
        self.coverage = self.num_against+self.num_for
    def test(self,what):
        if what == "multiregion" :
            if ( self.type == "<TRA>" and
                    self.coverage > 10 and
                    self.coverage < 100 and
                    self.allele == "0/1" ) :
                return True
    def printBedpe(self,name,out) :
        n=500
        print("\t".join([str(x) for x in [self.chrom,self.pos-n,self.pos+n,
            self.info["CHR2"],self.info["END"]-n,self.info["END"]+n,
            name,".",".",self.info["RNAMES"]]]),file=out)
    def printBed(self,name,out) :
        n=500
        print("\t".join([
            str(x) for x in [
                self.chrom,self.pos-n,self.pos+n,name,
                ".",".",self.info["RNAMES"]]
            ]),file=out)

        print("\t".join([
            str(x) for x in [
                self.info["CHR2"],self.info["END"]-n,self.info["END"]+n,
                name,".",".",self.info["RNAMES"]]
            ]),file=out)

if __name__=="__main__":
    args=parseArgs()
#    print(args)
#    sys.exit()
    # removing commented lines
    cmtflag = 1
    while cmtflag == 1 :
        line = args.input.readline()
        if line[0] != "#" : cmtflag = 0
    i = 1
    for line in args.input :
        entry = SnifflesEntry(line)
        if entry.test(args.what) : 
            entry.printBed(i,args.output)
            i+=1
