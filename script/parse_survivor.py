#! /usr/bin/env python
import sys
import os
import argparse
import gzip
import pysam
import numpy as np
from collections import namedtuple
from methylbed_utils import SnifflesEntry,make_coord,read_bam
import multiprocessing as mp
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='dealing with SURVIVOR vcf')
    subparsers = parser.add_subparsers()
    # parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-v', '--verbose', action='store_true',default=False,
            help="verbose output")
    parent_parser.add_argument('-s', '--sniffles', type=os.path.abspath,
            required=False, help="input vcf file path")
    parent_parser.add_argument('-o', '--output', type=argparse.FileType('w'),
            required=False, default=sys.stdout,help="output file path (or stdout)")
    parent_parser.add_argument('-f', '--format', type=str,
            required=False, default="bed",help="output file format (default : bed)")
    parent_parser.add_argument('-r','--readsupport',type=int,required=False,
            default=5,help="minimum read suppport (default 5)")
    parent_parser.add_argument('--one',type=str,required=True,
            help="Name of sample one in comparison")
    parent_parser.add_argument('--two',type=str,required=True,
            help="Name of sample two in comparison")
    parent_parser.add_argument('--heterozygous', action='store_true',
            required=False, default=False,help="report het svs too (default False)")
    # parser for del
    parser_del = subparsers.add_parser('del',parents=[parent_parser],
            help = "deletions")
    parser_del.add_argument('-l','--length',type=int,required=False,
            default=500,help="minimum length (default 500)")
#    parser_del.add_argument('-w','--window',type=int,required=False,
#            default=200,help="window for flanking regions")
    parser_del.set_defaults(what="del",func=parseDeletions)
    # insertion
    parser_ins = subparsers.add_parser('ins',parents=[parent_parser],
            help = "insertions")
    parser_ins.add_argument('-l','--length',type=int,required=False,
            default=200,help="minimum length (default 500)")
    parser_ins.set_defaults(what="ins",func=parseInsertions)
    # tra
    parser_ins = subparsers.add_parser('tra',parents=[parent_parser],
            help = "translocations")
    parser_ins.set_defaults(what="tra",func=parseTranslocations)

    args = parser.parse_args()
#    print(args)
    return args

def svcoord_to_bedcoord(in_coord) :
    coord_split = in_coord.split("-")
    chroms = [ x.split("_")[0] for x in coord_split ]
    pos = [ int(x.split("_")[1]) for x in coord_split ]
    return [chroms[0],pos[0]-1,pos[1]]

def parseDeletions(fields,one_ind,two_ind,zyg_test,sv_id,args) :
    samp_fields = [ fields[i] for i in [one_ind,two_ind] ]
    if "DEL" not in ",".join(samp_fields) : return sv_id # test for SV type
    info = [ x.split(":") for x in samp_fields ]
    # test for zygosity
    zygosity = [ x[0] for x in info ]
    sv_ind = zyg_test(zygosity)
    if sv_ind == -1 : return sv_id
    # test for read support
    sv_info = info[sv_ind]
    re = int(sv_info[3].split(",")[-1])
    if re < args.readsupport : return sv_id
    # test for length
    sv_sample = [args.one,args.two][sv_ind]
    coords_list = sv_info[-1].split(",")
    bedcoords_list = [svcoord_to_bedcoord(x) for x in coords_list ]
    lengths = [x[2]-x[1] for x in bedcoords_list ]
    if not any( x > args.length for x in lengths ) : return sv_id
    bed_filt = [ x for i,x in enumerate(bedcoords_list) 
            if lengths[i] > args.length ] 
    for x in bed_filt :
        sv_id += 1
        out_str = '\t'.join([str(y) for y in 
            x + [str(sv_id)+"_del",sv_sample]])
        print(out_str,file=args.output) 
    return sv_id

def parseInsertions(fields,one_ind,two_ind,zyg_test,sv_id,args) :
    samp_fields = [ fields[i] for i in [one_ind,two_ind] ]
    if "INS" not in ",".join(samp_fields) : return sv_id # test for SV type
    info = [ x.split(":") for x in samp_fields ]
    # test for zygosity
    zygosity = [ x[0] for x in info ]
    sv_ind = zyg_test(zygosity)
    if sv_ind == -1 : return sv_id
    # test for read support
    sv_info = info[sv_ind]
    re = int(sv_info[3].split(",")[-1])
    if re < args.readsupport : return sv_id
    # test for length
    svlen = int(sv_info[2])
    if not  svlen > args.length : return sv_id
    sv_sample = [args.one,args.two][sv_ind]
    sv_id += 1
    outfields = [fields[0],int(fields[1])-1,fields[1],str(sv_id)+"_ins",sv_sample]
    bedcoord = '\t'.join([str(x) for x in outfields])
    print(bedcoord,file=args.output) 
    return sv_id

def parseTranslocations(fields,one_ind,two_ind,zyg_test,sv_id,args) :
    def coord_to_bed(coord) :
        chrom,end = coord.split("_")
        start = str(int(end)-1)
        return [chrom,start,end]
    op_strand = {"+":"-","-":"+"}
    samp_fields = [ fields[i] for i in [one_ind,two_ind] ]
    if "TRA" not in ",".join(samp_fields) : return sv_id # test for SV type
    info = [ x.split(":") for x in samp_fields ]
    # test for zygosity
    zygosity = [ x[0] for x in info ]
    sv_ind = zyg_test(zygosity)
    if sv_ind == -1 : return sv_id
    # test for read support
    sv_info = info[sv_ind]
    re = int(sv_info[3].split(",")[-1])
    if re < args.readsupport : return sv_id
    sv_sample = [args.one,args.two][sv_ind]
    sv_id += 1
    coords_list = sv_info[-1].split(",")
    strands = sv_info[4]
    if len(coords_list) == 1 :
        coords = coords_list[0].split("-")
        bp1 = "\t".join(coord_to_bed(coords[0])+
                [str(sv_id)+"_tra",sv_sample,strands[0]])
        print(bp1,file=args.output)
        bp2 = "\t".join(coord_to_bed(coords[1])+
                [str(sv_id)+"_tra",sv_sample,op_strand[strands[1]]])
        print(bp2,file=args.output)
        # need to work on cases when there are multiple coords
    return sv_id

if __name__=="__main__":
    args=parseArgs()
    in_fh = gzip.open(args.sniffles,'rt')
    # different test depending on zygosity 
    if not args.heterozygous :
        def zyg_test(zyg_list) :
            zyg_str = ",".join(zyg_list)
            if (any(tag in zyg_str for tag in [ "0/0","./." ]) and
                    "1/1" in zyg_str) :
                return zyg_list.index("1/1")
            else :
                return -1
    # remove commented lines and figure out what fields to use
    for line in in_fh :
        if line[0:2] == "##" : continue
        if line[1:6] == "CHROM" :
            fields = line.split("\t")
            for i,s in enumerate(fields) :
                if args.one in s : one_ind = i
                elif args.two in s : two_ind = i
            break
    # parsing
    i = 0
    for line in in_fh :
        fields = line.strip().split("\t")
        i = args.func(fields,one_ind,two_ind,zyg_test,i,args)
    in_fh.close()
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)
    sys.exit()
