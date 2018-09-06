#!/bin/bash
cells="MCF10A MCF7 MDAMB231"
cells="GM12878"
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
bed=$root/database/hg38/hg38_repeats_LINE_long.bed
reg=LINE_long
pooldir=$root/pooled
for cell in $cells;do
  bam=$pooldir/bam/$cell.pooled.bam
  mbeddir=$pooldir/methylation/methbyread_all
  cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
  script=../../nanopolish/convertBam.py
  log=$pooldir/igv/$cell.$reg.log
  outsort=$pooldir/igv/$cell.$reg.sorted.bam

  com="python -u $script -v -t 10 \
    -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
    samtools sort -o $outsort"
  echo $com
  eval $com
  samtools index $outsort
done
