#! /usr/bin/env python
import sys
import os
import math

if sys.argv[1] :
    in_fh = open(sys.argv[1])
else :
    in_fh = sys.stdin

for line in in_fh :
    fields = line.split("\t")
    center = math.floor((int(fields[1])+int(fields[2])-1)/2)
    fields[1] = str(center)
    fields[2] = str(center+1)
    print("\t".join(fields),file=sys.stdout)

in_fh.close()
