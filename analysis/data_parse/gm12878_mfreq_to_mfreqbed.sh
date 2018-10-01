#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation
mfreqdir=$root/mfreq_all
outdir=$root/mfreq_bed
script=../../script/methylFreq_destrand.py
[ -e $outdir ]||mkdir $outdir
mods="cpg gpc"
bases="GM12878 GM12878_BSseq_1"
bases="GM12878_BSseq_1"


for base in $bases;do
  for mod in cpg gpc;do
    mfreq=$(find $mfreqdir -name "$base.$mod.methfreq.txt.gz")
    if [ "$base" == "GM12878_BSseq_1" ];then
      subcom="gunzip -c $mfreq | python $script"
    else
      subcom="gunzip -c $mfreq | grep -v GCG"
    fi
    out=$outdir/$base.$mod.methfreq.bedGraph
    eval $subcom | awk 'OFS="\t"{ print $1,$2-1,$2,$4,$5 }' > $out
  done
done
  

