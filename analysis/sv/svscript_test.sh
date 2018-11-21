#!/bin/bash
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
methdir="$root/pooled/methylation/methbyread_all"
bamdir=$root/pooled/bam
script=../../script/SVmethylation.py
cells=GM12878
cells="GM12878 MCF10A MCF7 MDAMB231"
for cell in $cells; do
  sv=$(find $root/pooled/sv/ -name "$cell.pooled*vcf.gz")
  bam=$(find $bamdir -name "$cell*pooled.bam")
  cpg=$(find $methdir -name "$cell.*cpg*pooled*bed.gz")
  gpc=$(find $methdir -name "$cell.*gpc*pooled*bed.gz")

  gunzip -c $sv | $script -t 8 -v -b $bam -c $cpg -g $gpc
done
