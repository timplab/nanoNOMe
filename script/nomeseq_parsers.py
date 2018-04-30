import math
import os
import csv
from collections import namedtuple
import re
import numpy as np
import pandas as pd

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
                    self.ratios[i],
                    self.seqs[i])
        return calldict
    def getArray(self,calldict) :
        callarray=np.array([(x,calldict[x].call) for x in sorted(calldict.keys())])
        return callarray
    # this fxn is not complete
    def setRollmean(self,win) :
        metharray=np.array(self.callarray)
        start=metharray[0,0]
        end=metharray[-1,0]
        # initialize array
        raw_array=np.zeros(end-start+1)
        raw_array[:]=NA
        # replace values
        ind=metharray[:,0]-start
        meth=np.array(metharray[:,1])
#        meth[metharray[:,1]==0]=-1
#        meth[metharray[:,1]==-1]=0
        np.put(raw_array,ind,meth)
        smooth_array=pd.rolling_mean(raw_array,win,min_periods=1)
#        smooth_array=pd.Series(self.callarray[:,1],index=self.callarray[:,0])

        self.smoothmeth=smooth_array
        return

# for a line in intersect result
class Region :
    def __init__(self,l) :
        self.items=l
        self.rname=l[0]
        self.start=int(l[1])
        self.end=int(l[2])
        self.name=l[3]
        self.score=l[4]
        self.strand=l[5]
        self.id=l[6]
        self.fxn=l[7]
        try : 
            self.center=int(l[8])
        except IndexError :
            pass

class methInt :
    def __init__(self,string):
        f=string.strip().split("\t")
        self.fields=f
        mread="\t".join(self.fields[:7])
        reg=self.fields[7:]
        self.methread=MethRead(mread)
        self.region=Region(reg[:8])
        try : 
            self.center=int(reg[8])
        except IndexError :
            pass

# parse frequency line
class methFreq : 
    def __init__(self,string):
        f=string.strip().split("\t")
        self.fields=f
        self.rname=f[0]
        self.start=int(f[1])
        self.end=int(f[2])
        self.methylated=int(f[3])
        self.unmethylated=int(f[4])
# parse distance to frequency line
class distFreq : 
    def __init__(self,string):
        f=string.strip().split("\t")
        self.fields=f
        self.mfreq=methFreq("\t".join(f[:5]))
        self.region=region(f[5],int(f[6]),int(f[7]),f[10],f[8])
        self.distance=int(f[-1])

