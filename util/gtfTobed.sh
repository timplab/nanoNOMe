#!/bin/bash

# changing to 0-based start, 1-based end

infile=$1
ann=$2
ann=${ann:=gene}

gunzip -c $1 | grep -v "#" | awk -v ann=$ann \
  '{ OFS="\t" }{ if($3==ann) print "chr"$1,$4-1,$5,$10,".",$7,$14,$18 }' \
  | tr -d '";' | sort -k1,1 -k2,2n
