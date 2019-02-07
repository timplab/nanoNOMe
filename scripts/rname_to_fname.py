#! /usr/bin/env python
import sys
import argparse
import os
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='get fname from rnames')
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i','--index',type=argparse.FileType('r'),required=True,
            help="index file")
    parser.add_argument('-r','--readnames',type=argparse.FileType('r'),required=False,
            default=sys.stdin,help="read names")
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args = parseArgs()
    rnames = [ x.strip() for x in args.readnames ]
    if args.verbose : 
        print("extracting",file=sys.stderr)
    for line in args.index :
        fields = line.strip().split()
        if fields[0] in rnames :
            print(fields[-1],file=sys.stdout)
    if args.verbose : 
        print("elapsed time : {} seconds".format(
            time-time()-start_time),file=sys.stderr)
