#!/bin/bash

root=/dilithium/Data/Nanopore/Analysis/171205_nomeseq
wdir=$root/$mod
samp=dam

if [ 1 -eq 1 ];then
  for mod in dam; do
    echo "calculate methylation frequency for $mod"
    methtsv=`find $wdir/$samp/methcall-b37 -name "*$mod.meth.tsv"`
    base=$(basename $methtsv)
    pre=${base%.meth.tsv}
    echo $pre
    freqout=$wdir/$samp/methcall-b37/$pre.meth.freq.tsv
    python ~/Code/ilee/nanopolish/calculate_methylation_frequency.py \
      -i $methtsv -m $mod > $freqout
  done
fi
