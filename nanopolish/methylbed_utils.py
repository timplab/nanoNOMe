import math
import os
import csv
from collections import namedtuple
import re
import numpy as np

# for read from a methylation bed file
methylCall = namedtuple('methylCall', ['pos','call','ratio','seq'])
class MethRead :
    def __init__(self,string):
        self.fields=string.strip().split("\t")
        self.rname=self.fields[0]
        self.start=int(self.fields[1])
        self.end=int(self.fields[2])
        self.rlen=self.end-self.start
        self.qname=self.fields[3]
        self.methstring=self.fields[4]
        self.ratios=self.fields[5].strip().split(",")
        self.seqs=self.fields[6].strip().split(",")
        self.calldict=self.parseMeth()
        self.keys=sorted(self.calldict.keys())
        self.callarray=self.getArray(self.calldict)
    def make_key(self,pos):
        return int(pos)
    def parseMeth(self):
        calldict=dict()
        pos=int(self.start)
        calls=re.findall('(\d+)([umx])?',self.methstring)
        for i,(dist,call) in enumerate(calls):
            pos+=int(dist)
            if call=="x" :
                is_meth=-1
            else :
                is_meth=int(call=="m")
            calldict[self.make_key(pos)]=methylCall(
                    pos,
                    is_meth,
                    float(self.ratios[i]),
                    self.seqs[i])
        return calldict
    def getArray(self,calldict) :
        callarray=np.array([(x,calldict[x].call) for x in sorted(calldict.keys())])
        return callarray

# sniffles entry
class SnifflesEntry :
    def __init__(self,line) :
        self.line=line.strip()
        self.fields=self.line.split("\t")
        (self.chrom,self.pos,self.id,self.ref,
                self.type,self.qual,self.filter,self.infostring,
                self.format,self.genotype) = self.fields
        self.pos = int(self.pos)
        self.type = self.type.strip("<").strip(">")
    def activate(self) :
        self.parseinfo()
        self.parsegenotype()
    def parseinfo(self) :
        self.infofields = [ x.split("=") for x in self.infostring.strip().split(";")]
        self.info = dict()
        self.info["CONFIDENCE"] = self.infofields[0][0]
        for entry in self.infofields[1:] :
            self.info[entry[0]] = entry[1]
        self.info["END"] = int(self.info["END"])
        self.rnames = self.info["RNAMES"].split(',')
    def parsegenotype(self) :
        self.allele = self.genotype.split(":")[0]
        self.num_against = int(self.genotype.split(":")[1])
        self.num_for = int(self.genotype.split(":")[2])
        self.coverage = self.num_against+self.num_for

