#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/ctcf
ctcf=$root/GM12878_CTCF.hg19.bed
center=$root/GM12878_CTCF.hg19.center.bed
awk '{ $2=($2+$3-1)/2 ;$3=$2+1 ; printf "%s\t%i\t%i\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6 }' $ctcf |\
 sort -k1,1 -k2,2n  > $center
