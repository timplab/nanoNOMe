#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
cell=GM12878
bamdir=$root/pooled/bam
anndir=$root/database
hapdir=$root/pooled/hap
ref=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa
hcdir=~/Code/HapCUT2/build

bam=$(find $bamdir -name "*$cell*subset.bam")
var=$anndir/gm12878/variants/NA12878.vcf

hairout=$hapdir/$cell.hairs
if [ "$1" == "hair" ];then
  log=$hairout.log
  $hcdir/extractHAIRS --ONT 1 --ref $ref \
    --bam $bam --VCF $var --out $hairout &> $log
fi

hapout=$hapdir/$cell.haplotypes
if [ "$1" == "hapcut" ];then
  log=$hapout.log
  $hcdir/HAPCUT2 --ea 1 --fragments $hairout --VCF $var --output $hapout &> $log
fi

if [ "$1" == "prune" ];then
  python $hcdir/../utilities/prune_haplotype.py -i $hapout -o $hapout.pruned \
    --min_mismatch_qual 30 --min_switch_qual 30
fi
