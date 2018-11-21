#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
outdir=$root/pooled/methylation/methbyread_all
d1=$root/pooled/methylation/methbyread_old
d2=$root/gm12878/ngmlr/mbed
samp=GM12878
for mod in cpg gpc;do
  echo $mod
  f1=$(find $d1 -name "*$samp*$mod*gz")
  f2=$(find $d2 -name "181102*$samp*$mod*gz")
  out=$outdir/$samp.$mod.pooled.meth.bed.gz
  echo $out
  gunzip -c $f1 $f2 |\
    sort -k1,1 -k2,2n |\
    bgzip > $out
  tabix -p bed $out
done
