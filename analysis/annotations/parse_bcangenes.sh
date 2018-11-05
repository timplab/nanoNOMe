#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root="$1"
dbdir="$root/annotations/breastcancer"
prefix="$dbdir/bcan_genes"
gs="$root/annotations/hg38/hg38_genomesize.txt"

# parse genes
echo "subsetting bcan genes"
genes=$root/annotations/hg38/hg38_genes.bed
names=$(awk 'NR>1{ print $1 }' \
  $srcdir/../../annotations/bcan_genes_CancerGeneCensus.txt)
bed=$prefix.bed
grep -F "$names" $genes > $bed
# get TSS
tss=$prefix.TSS.bed
com="$srcdir/../../util/getTSS.sh $bed > $tss"
echo $com
eval $com
# get regions around TSS
for width in 400 2000 5000 10000;do
  side=$(($width/2))
  region=$prefix.TSS.${width}bp.bed
  com="bedtools slop -b $side -i $tss -g $gs|\
    sort -k1,1 -k2,2n > $region"
  echo $com
  eval $com
done
# repeats
echo "repetitive elements"
rep=$root/annotations/hg38/hg38_repeats.bed
proxrep=$dbdir/repeats_bcan_proximity.bed
bedtools closest -D b -a $rep -b $bed |\
  awk '{ if(sqrt($NF*$NF) <= 5000) print }' > $proxrep
bodyrep=$dbdir/repeats_bcan_body.bed
bedtools intersect -wo -a $rep -b $bed > $bodyrep

# for transcripts
echo "transcripts"
base=${prefix%%_*}
trans=$root/annotations/hg38/hg38_transcripts.bed
prefix=${base}_transcripts
bed=$prefix.bed
grep -F "$names" $trans > $bed

# get TSS
tss=$prefix.TSS.bed
com="$srcdir/../../util/getTSS.sh $bed > $tss"
echo $com
eval $com

# get regions around TSS
for width in 400 2000 5000 10000;do
  side=$(($width/2))
  region=$prefix.TSS.${width}bp.bed
  com="bedtools slop -b $side -i $tss -g $gs|\
    sort -k1,1 -k2,2n > $region"
  echo $com
  eval $com
done
