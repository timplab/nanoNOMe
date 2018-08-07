#!/bin/bash
sniffles=~/Code/Sniffles*/bin/sniffles*/sniffles
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled
bases="GM12878 MCF10A MCF7 MDAMB231"
#base="{}"
bam=$root/bam/$base.pooled.bam
outdir=$root/sv
[ -e $outdir ]||mkdir $outdir
if [ "$1" == "vcf" ];then
  echo "generating vcf"
  for base in $bases; do
    bam=$root/bam/$base.pooled.bam
    svout=$outdir/$base.sniffles.vcf
    $sniffles -m $bam -s 3 -n 30 \
      -v $svout --genotype --cluster --tmp_file $outdir/$base.tmp &> $svout.log
  done
fi
if [ "$1" == "bed" ];then
  echo "generating bedpe"
  for base in $bases; do
    echo $base
    bam=$root/bam/$base.pooled.bam
    svout=$outdir/$base.sniffles.bed
    $sniffles -m $bam -s 3 -n 30 \
      -b $svout --tmp_file $svout.tmp &> $svout.log
  done
fi
