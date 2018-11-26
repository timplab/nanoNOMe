#! /usr/bin/env python
import sys
import argparse
import os
import boto3
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='copy in files from aws s3')
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i','--input',type=argparse.FileType('r'),required=False,
            default=sys.stdin,help="input file names")
    parser.add_argument('-b','--bucket',type=str,required=False,
            help="bucket")
    parser.add_argument('-d','--dir',type=str,required=False,
            default="",help="directory within bucket")
    parser.add_argument('-s','--s3index',type=argparse.FileType('r'),required=False,
            help="s3 index, with field 1 as file name and field 2 as path within s3")
    parser.add_argument('-o','--outdir',type=str,required=False,
            default=".",help="directory to output files")
    args = parser.parse_args()
    assert ( args.bucket is not None or args.s3index is not None ),"either bucket or s3 index must be provided"
    return args

if __name__=="__main__":
    args=parseArgs()
    fnames = [ x.strip() for x in args.input ]
    key_dict = dict()
    if args.verbose : print("making a dict of the index",file=sys.stderr)
    if args.s3index is not None :
        for line in args.s3index :
            fields = line.strip().split()
            path=fields[1].split(":")
            key_dict[fields[0]] = path[1]
            args.bucket=path[0]
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(args.bucket)
    if args.s3index is None :
        if args.verbose : print("indexing s3 dir",file=sys.stderr)
        i = 0
        for obj in bucket.objects.filter(Prefix=args.dir):
            key_dict[obj.key.split('/')[-1]] = obj.key
            i += 1
            if i%10000 == 0 :
                if args.verbose : print(i,file=sys.stderr)
    if args.verbose : print("done making a dict of the index",file=sys.stderr)
    if args.verbose : print("downloading files",file=sys.stderr)
    i = 0
    for fn in fnames :
        outpath=os.path.join(args.outdir,fn)
        try :
            bucket.download_file(key_dict[fn],outpath)
            i +=1
        except KeyError :
            print("{} does not exist in the specified location".format(fn),file=sys.stderr)
        if args.verbose : 
            if i%10000 == 0 : print(i,file=sys.stderr)
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)


