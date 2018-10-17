#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/
#db=$root/database/hg38/hg38_repeats_SINE.bed
#dbname=SINE
db=$root/database/hg38/hg38_repeats_LINE.bed
dbname=LINE
indir=$root/pooled/methylation/mfreq_bed
outdir=$root/intersect
mods="cpg gpc"
samps="nanonome bsseq"
for samp in $samps;do
  if [ "$samp" == "nanonome" ];then
    base=GM12878
  elif [ "$samp" == "bsseq" ];then
    base=GM12878_BSseq_1
  fi
  for mod in $mods;do
    pre="$pre $base.$mod"
  done
done
infp=$indir/{}.methfreq.bedGraph
out=$outdir/{}.$dbname.methfreq.bedGraph
if [ "$1" == "bam" ];then
  pre=""
  for samp in $samps;do
    if [ "$samp" == "nanonome" ];then
      samp=GM12878
      base=$root/pooled/bed/GM12878.pooled
    elif [ "$samp" == "bsseq" ];then
      samp=GM12878_BSseq_1
      base=/dilithium/Data/NGS/projects/gm12878/bsseq/bed/GM12878_BSseq_ENCLB794YYH.nodup
    fi
    pre="$pre $base"
  done
  infp={}.bed
  out={}.$dbname.bedGraph
fi
echo $pre

com="bedtools intersect -sorted -wb -a $infp -b $db > $out"
echo $com
parallel "$com"  ::: $pre
