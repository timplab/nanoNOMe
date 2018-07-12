#!/bin/bash
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis/validation/scNOMe"

if [ "$1" == "fq" ];then
  yield=$root/fastqyield.txt
  zcat $root/fastq/GM12878_sample_rep*gz |\
    grep length |\
    awk 'FS="="{ sum+=$2 }END{ print NR,sum }' > $yield
fi
bamdir=$root/bam
beddir=$root/bed
if [ "$1" == "bam" ];then
  bamyield=$root/qc/bamyield.csv
  [ -e $bamyield ]&&rm $bamyield
  for bam in $(find $bamdir -name "*pooled*bam");do
    base=$(basename "$bam")
    samp=${base%%.*}
    tmpdir=$root/tmp/bed/$samp
    [ -e $tmpdir ]||mkdir $tmpdir -p
    echo $samp
    bed=$beddir/$samp.pooled.bed.gz
    [ -e $bed ]||\
      bedtools bamtobed -i $bam |\
      sort -T $tmpdir -k1,1 -k2,2n |\
      bgzip > $bed
    qc=$(gunzip -c $bed |\
      awk 'OFS=","{ sum+=$3-$2 }END{ print NR,sum }')
    echo $samp,$qc >> $bamyield
  done
fi

