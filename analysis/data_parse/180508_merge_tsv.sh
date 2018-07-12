#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
pooldir=$root/pooled/methylation/meth
catscr=/home/isac/Code/ilee/util/catFiles.py

for cell in GM12878;do
  dir=${root}/${cell,,}/ngmlr
  for mod in cpg gpc; do
    mcalldir=$dir/mcall-$mod
    pooltsv=$pooldir/$cell.$mod.pooled.meth.bed.gz
    echo $pooltsv
    tsvs=`find $mcalldir -name "*$cell*$mod*.meth.tsv"`
    awk 'FNR>1{ print }' $tsvs |\
      sort -T $pooldir -k1,1 -k2,2n |\
      bgzip > $pooltsv
    tabix -p bed $pooltsv
  done
done

