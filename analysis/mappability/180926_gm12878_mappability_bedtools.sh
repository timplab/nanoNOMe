#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/
bismap=$root/database/hg38/bismap/hg38/k100.Bismap.MultiTrackMappability.bedGraph
indir=$root/pooled/methylation/mfreq_bed
outdir=$root/mappability_bismap
mods="cpg gpc"
mods="cpg"
samps="nanonome bsseq"
samps="bsseq"
for samp in $samps;do
  if [ "$samp" == "nanonome" ];then
    base=GM12878
  elif [ "$samp" == "bsseq" ];then
    base=GM12878_BSseq_1
  fi
  infp=$indir/$base.{}.methfreq.bedGraph
  out=$outdir/$base.{}.bismap.bedGraph
  parallel "bedtools intersect -sorted -loj -a $infp -b $bismap > $out" ::: $mods
done

