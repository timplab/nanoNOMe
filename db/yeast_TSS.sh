#!/bin/bash
dir=/home/isac/Dropbox/Data/nome-seq/db/yeast
tss=$dir/SacCer3_TSS.txt
bed=$dir/SacCer3_TSS.bed

awk 'BEGIN{ OFS="\t" }NR>1{ if($5==1){ print $1,$2-1,$2,$6,".","+" }else { print $1,$3-1,$3,$6,".","-" }}' $tss |\
  sort -k1,1 -k2,2n > $bed



