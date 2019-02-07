#! /usr/bin/env python
import os
import sys
import gzip
import pysam
import multiprocessing as mp
import time
start_time = time.time()

# https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def listener(q,out,verbose=False) :
    '''listens for messages on the q, writes to file. '''
    if verbose : print("writing output to {}".format(out),file=sys.stderr)
    if out == "stdout" :
        out_fh = sys.stdout
    else :
        out_fh = open(out,'w')
    reg_list= list()
    while 1:
        m = q.get()
        if m == 'kill':
            break
        key,line = m
        if key not in reg_list:
            reg_list.append(key)
            print(line,file=out_fh)
        else :
            if verbose: print("{} already written - skipping".format(key),file=sys.stderr)
        q.task_done()
    if verbose : print("total {} lines written".format(len(reg_list)),file=sys.stderr)
    out_fh.close()
    q.task_done()

def chunkify(fname,size=1024*1024):
    """
    partition file into small chunks for each process
    """
    fileEnd = os.path.getsize(fname)
    with open(fname,'rb') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                f.close()
                break

def init_mp(threads) :
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(threads)
    return manager,q,pool

def close_mp(q,pool) :
    q.put('kill')
    q.join()
    pool.close()
