#!/bin/bash
export LC_ALL=C
inbed="$1"
gs="$2"
prefix=${inbed%.bed}

for width in 40 2000 5000 10000; do
  side=$(($width/2))
  outbed=$prefix.${width}bp.bed
  bedtools slop -b $side -i $inbed -g $gs |\
    sort -k1,1 -k2,2n > $outbed
done
