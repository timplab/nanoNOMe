#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
if [ ! -z $1 ];then
  root=$1
fi
dir=$root/database/gm12878/mnase
samp=GM12878
bigwig=$dir/*bigWig
converter=../../util/bigWigToBedGraph

bed=$dir/${samp}_mnase.bedGraph
if [ ! -e $bed.gz ];then
  $converter $bigwig $bed
  bgzip $bed ; tabix -p bed $bed.gz
fi
bedhg38=$dir/${samp}_mnase.hg38.bedGraph

if [ ! -e $bedhg38.gz ];then
  chain=/mithril/Data/NGS/Reference/human/liftover/hg19ToHg38.over.chain
  echo "liftover"
  liftOver $bed.gz $chain $bedhg38 $bedhg38.unmapped.bed
  echo "sort"
  sort -k1,1 -2,2n -T $dir $bedhg38 | bgzip 
  tabix -p bed $bedhg38.gz
fi


