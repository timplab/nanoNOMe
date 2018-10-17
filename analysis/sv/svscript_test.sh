#!/bin/bash
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
methdir="$root/pooled/methylation/methbyread_all"
bamdir=$root/pooled/bam
sv=$root/pooled/sv/SVcomparison.vcf
script=../../script/SVmethylation.py
cell=MDAMB231
bam=$(find $bamdir -name "$cell*pooled.bam")
cpg=$(find $methdir -name "$cell*cpg*pooled*bed.gz")
gpc=$(find $methdir -name "$cell*gpc*pooled*bed.gz")

grep -E 'xxxo|xxoo' $sv | $script -t 8 -v -b $bam -c $cpg -g $gpc
