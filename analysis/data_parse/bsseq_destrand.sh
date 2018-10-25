#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
mfreqs=$(find $root/pooled/methylation/mfreq_all -name "*BSseq_1.*gz")
parser=../../script/methylFreq_destrand.py

for m in $mfreqs;do
  echo $m
  base=$(basename "$m")
  dir=$(dirname "$m")
  base=${base%%.*}
  mod=${m#*$base.}
  mod=${mod%%.*}
  out=$dir/${base}_destrand.$mod.methfreq.txt.gz
  echo $out
  gunzip -c $m | python $parser -v | bgzip > $out
  tabix -b 2 -e 2 $out
done
