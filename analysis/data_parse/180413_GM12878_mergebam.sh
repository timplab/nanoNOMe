#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/gm12878/ngmlr
bamdir=$root/bam
poolbam=$bamdir/GM12878.pooled.bam

if [ "$1" == "within" ];then
  for dir in `find $bamdir/* -maxdepth 0 -type d`;do
    base=$(basename "$dir")
    sampbam=$bamdir/$base.bam
    bams=`find $dir -name "*sorted.bam"`
    echo $dir
    com="sambamba merge -p -t 10 $sampbam $bams"
  done
fi
if [ "$1" == "between" ];then
  bams=`find $bamdir/* -maxdepth 0 -type f -name "1803*bam"`
  com="sambamba merge -p -t 10 $poolbam $bams"
fi

echo $com
#$com
