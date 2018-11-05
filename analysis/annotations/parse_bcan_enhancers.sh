#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
echo "parsing MCF10A enhancers"
root="$1"
PRE=$root/MCF10A_enhancer

# first convert to bed
echo "covert to bed"
awk 'OFS="\t"{ print $1,$2,$3,$7,".","." }' ${PRE}_hg19.txt \
  > ${PRE}_hg19.bed

# liftover 
echo "liftover to hg38"
chain=$root/../hg38/hg19ToHg38.over.chain
bed=${PRE}_hg38.bed
liftOver ${PRE}_hg19.bed $chain $bed.unsorted $bed.unmapped
sort -k1,1 -k2,2n $bed.unsorted > $bed
rm $bed.unsorted $bed.unmapped

# corresponding transcripts
trans=$root/bcan_transcripts.bed
Rscript $srcdir/match_enhancer_promoter.R $root 2> /dev/null


