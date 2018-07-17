#! /usr/bin/env python
import sys
import os
import argparse
import gzip
import math

def parseArgs():
    parser = argparse.ArgumentParser( description='dealing with bed file')
    subparsers = parser.add_subparsers()
    # parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-b','--bed',type=os.path.abspath,
            help="input bed file (default stdin)")
    parent_parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),required=False,
            default=sys.stdout,help="output bed path (default stdout)")
    parent_parser.add_argument('-v', '--verbose', action='store_true',default=False,
            help="verbose output")
    # parser for center
    parser_center = subparsers.add_parser('getcenter',parents=[parent_parser], 
            help = 'get center')
    parser_center.set_defaults(func=get_center)
    # parser for start
    parser_start = subparsers.add_parser('getstart',parents=[parent_parser],
            help = 'get start in strand-specific manner')
    parser_start.set_defaults(func=get_start)
    # parser for region
    parser_region = subparsers.add_parser('region',parents=[parent_parser],
            help = 'get region')
    parser_region.add_argument('-u','--upstream',type=int,required=False,
            default=1000,help="bp upstream")
    parser_region.add_argument('-d','--downstream',type=int,required=False,
            default=1000,help="bp downstream")
    parser_region.add_argument('--ignore_strands',action='store_true',required=False,
            default=False,help="ignore strands in determining upstream vs downstream")
    parser_region.set_defaults(func=get_region)
    args = parser.parse_args()
    return args

class BedQuery :
    def __init__(self,bedline) :
        self.fields=bedline.strip().split("\t")
        self.chrom=self.fields[0]
        self.start=int(self.fields[1])
        self.end=int(self.fields[2])
        self.name=self.fields[3]
        self.score=self.fields[4]
        self.strand=self.fields[5]
    def update_fields(self) :
        self.fields=[self.chrom,self.start,self.end,
                self.name,self.score,self.strand]+self.fields[6:]
    def printbed(self,out) :
        print("\t".join([str(x) for x in self.fields]),file=out)

def get_center(bed,args) :
    center=math.floor((bed.start+bed.end-1)/2)
    bed.start=center
    bed.end=center+1
    bed.update_fields()
    return bed

def get_start(bed,args) :
    if bed.strand == "-" :
        bed.start=bed.end-1
    else :
        bed.end=bed.start+1
    bed.update_fields()
    return bed

def get_region(bed,args) :
    bed.fields[6:]=[bed.start,bed.end]
    start_original=bed.start 
    bed.start=start_original-args.upstream-1
    bed.end=start_original+args.downstream+1
    bed.update_fields()
    return bed

if __name__=="__main__":
    args=parseArgs()
    if args.bed :
        if args.bed.split(".")[-1] == "gz" :
            bedfh = gzip.open(args.bed,'rt')
        else : bedfh = open(args.bed,'r')
    else : bedfh=sys.stdin
    for line in bedfh :
        bedline=BedQuery(line)
        newbed=args.func(bedline,args)
        if newbed.start >= 0 :
            newbed.printbed(args.out)

