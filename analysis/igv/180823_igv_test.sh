#!/bin/bash
cells="GM12878"
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
bed=$root/database/hg38/hg38_repeats_LINE_long.bed
reg=LINE_long
pooldir=$root/pooled
script=../../nanopolish/convertBam.py
for cell in $cells;do
  bam=$pooldir/bam/$cell.pooled.bam
  mbeddir=$pooldir/methylation/methbyread_all
  cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz

  if [ "$1" == "meth" ];then
    log=$pooldir/igv/$cell.$reg.log
    out=$pooldir/igv/$cell.$reg.sorted.bam
    com="python -u $script -v -t 8 \
      -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
      samtools sort -o $out"
  fi
  if [ "$1" == "original" ];then
    out=$pooldir/igv/$cell.$reg.original.sorted.bam
    com="samtools view -@ 10 -hb -L $bed $bam |\
      samtools sort -@ 10 -o $out"
  fi
  eval $com
  samtools index $out
done
