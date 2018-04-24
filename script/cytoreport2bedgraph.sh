#!/bin/bash

infile=$1

# for mat is chr,start,end,strand,meth,unmeth,context,actualcontext
gunzip -c $infile |\
  awk '{ OFS="\t" }{ print $1,$2-1,$2,$3,$4,$5,$6,$7 }' |\
  sort -k1,1 -k2,2n -T /dilithium/Data/tmp/
