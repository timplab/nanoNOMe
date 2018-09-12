#!/usr/bin/python
import os
import sys
root = "/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/sv"
cells = [ "GM12878","MCF10A","MCF7","MDAMB231" ]
vcfs = [ root+"/"+x+".sniffles.vcf" for x in cells ]
out = os.path.join(root,"SVcomparison.vcf")
n = 3 # rounding digits

def read_vcf(fp) :
    outdict = dict()
    i = 0
    with open(fp,'r') as fh :
        for line in fh :
            if line[0] == "#" : 
                continue
            i+=1
            fields = line.strip().split("\t")
            start = str(round(int(fields[1])/10**n))
            infolist = fields[7].split(';')
            key = '.'.join([fields[0],start,fields[4],
                infolist[2],infolist[3][:len(infolist[3])-n]])
            outdict[key] = line.strip()
    return outdict

print("reading data",file=sys.stderr)
data = [ read_vcf(x) for x in vcfs ]

cseq = list(range(len(cells)))
keylist = [ key for d in data for key in d.keys() ]
keys = set(keylist)
print("{} unique SVs out of total {}".format(len(keys),len(keylist)),file=sys.stderr)

out_fh = open(out,'w')
for key in keys :
    tag = ["x"]*len(cells)
    for i in cseq :
        try : 
            line = data[i][key]
            tag[i] = "o"
        except : pass
    print(line+"\t"+''.join(tag),file=out_fh)
out_fh.close()
