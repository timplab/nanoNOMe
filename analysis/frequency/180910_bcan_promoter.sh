#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
plotdir=$root/plots/promoters
bamdir=$root/pooled/bam
beddir=$root/pooled/methylation/methbyread_all
script=../../nanopolish/convertBam.py
regname=mdaVS10a_promoters
reg=$plotdir/$regname.bed
cells="MCF10A MCF7 MDAMB231"

# Rscript 180910_bcan_promoter.R

for cell in $cells; do
  bam=$bamdir/$cell.pooled.bam
  cpg=$beddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$beddir/$cell.gpc.pooled.meth.bed.gz
  log=$plotdir/$cell.$regname.log
  out=$plotdir/$cell.$regname.sorted.bam
  com="python -u $script -v -t 10 \
    -b $bam -c $cpg -g $gpc -r $reg 2> $log |\
    samtools sort -o $out"
  echo $com
  eval $com
  samtools index $out
done
