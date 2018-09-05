#! /usr/bin/env python
import sys
import os
import argparse
import pysam
from collections import namedtuple
import multiprocessing as mp
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='dealing with sniffles vcf ')
    subparsers = parser.add_subparsers()
    # parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-v', '--verbose', action='store_true',default=False,
            help="verbose output")
    parent_parser.add_argument('-t','--threads',type=int,default=1,
            help="threads (default 1)")
    parent_parser.add_argument('-i', '--input', nargs='?', type=argparse.FileType('r'),
            required=False, default=sys.stdin,help="input vcf file path (or stdin)")
    parent_parser.add_argument('-o', '--output', nargs='?', type=str,
            required=False, default="stdout",help="output file path (or stdout)")
    # parser to print predicted bed coordinate ranges of SVs
    parser_bed = subparsers.add_parser('bed',parents=[parent_parser])
    parser_bed.add_argument('what',help="element to extract")
    parser_bed.set_defaults(func=snifflesBed)
    # parser to subset bam entries 
    parser_bam = subparsers.add_parser('bam',parents=[parent_parser])
    parser_bam.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help = "bam file to subset")
    parser_bam.set_defaults(what="bam",func=snifflesBam)
    args = parser.parse_args()
#    print(args)
    return args

def make_coord(chrom,start,end) :
    return chrom+":"+str(start)+"-"+str(end)

# https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def listener(q,out,inbam=None,what="bam",verbose=False) :
    '''listens for messages on the q, writes to file. '''
    if verbose : print("writing output to {}".format(out),file=sys.stderr)
    in_fh = pysam.AlignmentFile(inbam,'rb')
    if out == "stdout" :
        out_fh = sys.stdout
        print(in_fh.header,file=out_fh)
        def printline(m,out_fh) :
            print(m,file=out_fh) 
    else :
        out_fh = pysam.AlignmentFile(out,'wb',template=in_fh)
        def printline(m,out_fh) :
            read = pysam.AlignedSegment.fromstring(m,out_fh.header)
            out_fh.write(read)
    in_fh.close()
    key_list = list()
    while 1:
        m = q.get()
        if m == 'kill':
            break
        key,line = m
        if key not in key_list:
            key_list.append(key)
            printline(line,out_fh)
        else :
            if verbose: print("{} already written - skipping".format(key),file=sys.stderr)
        q.task_done()
    if verbose : print("total {} lines written".format(len(key_list)),file=sys.stderr)
    out_fh.close()
    q.task_done()


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
        self.rnamelist = self.info["RNAMES"].split(",")
    def parsegenotype(self) :
        self.allele = self.genotype.split(":")[0]
        self.num_against = int(self.genotype.split(":")[1])
        self.num_for = int(self.genotype.split(":")[2])
        self.coverage = self.num_against+self.num_for
    def test(self,what):
        if what == "multiregion" :
#            if ( self.type == "<TRA>" and
#                    self.coverage > 10 and
#                    self.coverage < 100 and
#                    allele == "0/1" ) :
            if ( self.type == "<TRA>" and
                    self.coverage > 10 and
                    self.coverage < 100 ) :
                return True
    def printBedpe(self,name,out) :
        n=500
        print("\t".join([str(x) for x in [self.chrom,self.pos-n,self.pos+n,
            self.info["CHR2"],self.info["END"]-n,self.info["END"]+n,
            name,".",".",self.info["RNAMES"]]]),file=out)
    def printBed(self,name,out) :
        n=500
        if self.pos < n : self.pos = 500
        if self.info["END"] < n : self.info["END"] = 500
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

def snifflesBed(args) :
    print("not yet implemented")

def snifflesBam(bamfn,svline,verbose,q) :
    win = 500
    sv = SnifflesEntry(svline)
    if sv.allele == "0/0" : return
    start = sv.pos-win
    end = sv.pos+win
    coord = make_coord(sv.chrom,start,end)
    with pysam.AlignmentFile(bamfn,"rb") as bam :
        bam_entries = [ x for x in bam.fetch(region=coord) ]
    for bam in bam_entries : 
        if bam.query_name in sv.rnamelist : 
            key = '.'.join([bam.query_name,bam.reference_name,str(bam.reference_start)])
            q.put((key,bam.to_string()))



if __name__=="__main__":
    args=parseArgs()
    # removing commented lines
    cmtflag = 1
    while cmtflag == 1 :
        line = args.input.readline()
        if line[0] != "#" : cmtflag = 0
    inlines = args.input.readlines()
    # initialze mp
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    if args.verbose : print("using {} parallel processes".format(args.threads),file=sys.stderr)
    # watcher for output
    watcher = pool.apply_async(listener,(q,args.output,args.bam,args.what,args.verbose))
    time.sleep(0.05)
    jobs = [ pool.apply_async(args.func,
        args = (args.bam,line,args.verbose,q))
        for line in inlines ]
    output = [ p.get() for p in jobs ]
    # done
    q.put('kill')
    q.join()
    pool.close()
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)



