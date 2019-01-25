#!/usr/bin/env python
# remove alternative contigs
import sys
import os

flag=0
for line in sys.stdin:
    if line[0] == ">": 
        if any(x in line for x in ["_","."]) : flag=0
        else : flag=1
    if flag==1:
        sys.stdout.write(line)
