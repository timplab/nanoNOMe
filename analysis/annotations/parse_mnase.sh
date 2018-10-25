#!/bin/bash
dir="$1"
prefix="$dir/GM12878_MNase"
bigwig=$prefix.bigWig

bed=$prefix.bedGraph
if [ ! -e $bed.gz ];then
  bigWigToBedGraph $bigwig $bed
  bgzip $bed ; tabix -p bed $bed.gz
fi
bedhg38=$dir/${samp}_mnase.hg38.bedGraph

if [ ! -e $bedhg38.gz ];then
  chain=/mithril/Data/NGS/Reference/human/liftover/hg19ToHg38.over.chain
  echo "liftover"
#  ~/liftOver $bed.gz $chain $bedhg38 $bedhg38.unmapped.bed
  echo "sort"
#  sort -k1,1 -k2,2n -T $dir $bedhg38 > $bedhg38.sorted
#  bgzip $bedhg38.sorted 
  mv $bedhg38.sorted.gz $bedhg38.gz
  tabix -p bed $bedhg38.gz
fi


