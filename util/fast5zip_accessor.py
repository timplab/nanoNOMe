#! /usr/bin/env python
# script for indexing reads to paths in tarballs - currently needs the fast5list.txt to work
import sys
import os
import gzip
import argparse
import tarfile
import datetime
from collections import namedtuple
import zipfile
import pysam

def parseArgs():
    parser = argparse.ArgumentParser( description='dealing with zip index for fast5s')
    subparsers = parser.add_subparsers()
    # parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-v', '--verbose', action='store_true',default=False,
            help="verbose output")
    parent_parser.add_argument('-z', '--zipfile', type=os.path.abspath,required=False,
            help="fast5 zip file") 
    # parser for zipping
    parser_zip = subparsers.add_parser('zip', parents = [parent_parser],
            help = 'zip fast5s' )
    parser_zip.add_argument('-d','--dir',type=os.path.abspath,required=True,
            help = 'parent directory for fast5 files')
    parser_zip.set_defaults(func=zip_main)
    # parser for indexing
    parser_index = subparsers.add_parser('index',parents=[parent_parser], 
            help = 'index the fast5s in a zip')
    parser_index.add_argument('-l', '--list', type=os.path.abspath,required=False,
            help="log or list of fast5 zip file") 
    parser_index.add_argument('-r', '--readdb', type=argparse.FileType('r'),required=True,
            help="readdb output from `nanpolish index`")
    parser_index.add_argument('-b', '--bed', type=str,required=False,
            help="(optional) bed file of reads (bamtobed); if not supplied just outputs readname : zip path")
    parser_index.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),required=False,
            default=sys.stdout,help="output index file path (optional)")
    parser_index.set_defaults(func=index_main)
    # parser for accessing
    parser_access = subparsers.add_parser('access',parents=[parent_parser],
            help = 'access fast5(s) in a zip')
    parser_access.add_argument('-r','--reads',type=str,required=False,
            help="reads to extract in the form of entries from index file; use 'stdin' for piped input")
    parser_access.add_argument('-i','--index',type=os.path.abspath,required=False,
            help="index bed file")
    parser_access.add_argument('-w','--window',type=str,required=False,
            help="window from index file to extract [chrom:start-end]")
    parser_access.add_argument('-o','--outdir',type=os.path.abspath,required=True,
            help='output file dir')
    parser_access.set_defaults(func=access_main)
    args = parser.parse_args()
    return args

def openFile(path) :
    if path.split(".")[-1] == "gz" :
        out = gzip.open(path,'rt')
    else :
        out = open(path,'r')
    return out

def zip_main(args) :
    t0 = datetime.datetime.now()
    if args.verbose : print("zipping {} into {}".format(args.dir,args.zipfile),file=sys.stderr)
    with zipfile.ZipFile(args.zipfile,'w',compression=zipfile.ZIP_DEFLATED) as zip_fh :
        i = 0
        for root,dirs,files in os.walk(args.dir) :
            for fn in files :
                i += 1
                inpath = os.path.join(root,fn)
                newpath = inpath[len(args.dir):].strip("/")
                zip_fh.write(inpath,arcname=newpath)
                if args.verbose : 
                    print(newpath,file=sys.stdout)
                    if i%10000 == 0 :
                        print("archived {} files".format(i),file=sys.stderr)
    t1 = datetime.datetime.now()
    print("finished archiving {} files in {} seconds".format(i,(t1-t0).total_seconds()),
            file=sys.stderr)

def index_main(args) : 
    assert ((args.zipfile is not None) or
            (args.list is not None)), "either zip file or list must be provided"
    # f5list from zip
    ts0 = datetime.datetime.now()
    if args.list :
        if args.verbose : print("reading in the zip list",file=sys.stderr)
        f5list = list()
        with openFile(args.list) as fh :
            i = 0
            for line in fh :
                i+=1
                fields = line.split()
                index = [i for i,s in enumerate(fields) if '.fast5' in s ]
                if index : f5list.append(fields[index[0]])
#                if i == 10000 : break
    else : 
        if args.verbose : print("getting list from zip",file=sys.stderr)
        print("not yet supported",file=sys.stderr)
        sys.exit()
    # generate dict from the list
    f5dict = dict(zip([ x.split("/")[-1] for x in f5list ],f5list))
    del f5list
    ts1 = datetime.datetime.now() 
    if args.verbose : print("{} seconds taken to read in zip file list".format(
        (ts1-ts0).total_seconds()),file=sys.stderr)
    # readdb
    if args.verbose : print("readding in readdb",file=sys.stderr)
    readdb = dict()
    i = 0
    for line in args.readdb :
        i+=1
        fields = line.strip().split("\t")
        readdb[fields[1].split("/")[-1]] = fields[0]
#        if i == 10000 : break
    ts2 = datetime.datetime.now()
    if args.verbose : print("{} seconds taken to read in readdb".format(
        (ts2-ts1).total_seconds()),file=sys.stderr)
    # get readname : zip path pairs
    if args.verbose : print("associating read names with path within zip",file=sys.stderr)
    indexdict = dict()
    for key in readdb.keys() :
        try :
            indexdict[readdb[key]] = f5dict[key]
        except KeyError : continue
    for x in [f5dict,readdb] : del x
    ts3 = datetime.datetime.now()
    if args.verbose : print("{} seconds taken to read in readdb".format(
        (ts3-ts2).total_seconds()),file=sys.stderr)
    # bed association
    if args.bed : 
        if args.verbose : print("Outputting bed with the zip path",file=sys.stderr)
        if args.bed == "stdin" :
            in_fh = sys.stdin
        else : 
            in_fh = openFile(args.bed)
        i = 0
        for line in in_fh :
            key = line.split("\t")[3]
            try :
                print("\t".join([line.strip(),indexdict[key]]),file=args.out)
                i+=1
            except KeyError : continue
#            if i == 10 : break
        in_fh.close()
    else : 
        if args.verbose : print("Outputting read name and zip path associations",file=sys.stderr)
        for key,val in indexdict.items() :
            print("\t".join([key,val]),file=args.out)
    ts4 = datetime.datetime.now()
    if args.verbose : 
        print("{} seconds taken to ouput bed coordinates with zip path".format(
            (ts4-ts3).total_seconds()),file=sys.stderr)
        print("{} seconds taken total".format(
            (ts4-ts0).total_seconds()),file=sys.stderr)


def access_main(args) :
    assert ((args.reads is not None) or
            (args.index is not None)), "either index or entries from it must be given"
    t0 = datetime.datetime.now()
    if args.reads:
        if args.verbose : print("extracting reads from list given",file=sys.stderr)
        if args.reads == "stdin" :
            fh = sys.stdin
        else :
            fh = openFile(args.reads)
        entries = [x for x in fh]
        fh.close()
    elif args.window :
        if args.verbose : print("extracting reads in region given",file=sys.stderr)
        with pysam.TabixFile(args.index) as tabix :
            entries = [x for x in tabix.fetch(args.window)]
    else :
        if args.verbose : print("extracting all reads in index file",file=sys.stderr)
        with pysam.TabixFile(args.index) as tabix :
            entries = [x for x in tabix.fetch()]
    paths=set()
    for read in entries :
        paths.add(read.strip().split("\t")[-1])
    if args.verbose : print("extracting {} files into {}".format(len(paths),args.outdir))
    i = 0
    with zipfile.ZipFile(args.zipfile,'r') as zip_fh :
        for read in paths :
            i += 1 
            zip_fh.extract(read,path=args.outdir)
            if args.verbose and i%1000 == 0:
                print("extracted {} files".format(i),file=sys.stderr)
    t1 = datetime.datetime.now()
    if args.verbose : 
        print("finished extracting in {} seconds".format((t1-t0).total_seconds()),
                file=sys.stderr)
        
    
if __name__=="__main__":
    args=parseArgs()
#    print(args)
#    sys.exit()
    args.func(args)
