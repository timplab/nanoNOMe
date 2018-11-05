#!/bin/bash
dir=$(dirname $(readlink -f $0))
# parse genes
echo "parsing genes"
db="$1"
prefix=${db%%.*}
dbdir=$(dirname "$db")
gs=$dbdir/../hg38/hg38_genomesize.txt

# convert to bed
bed=$prefix.bed
com="$dir/../../util/gtfTobed.sh $db > $bed"
echo $com
eval $com

# get TSS
tss=$prefix.TSS.bed
com="$dir/../../util/getTSS.sh $bed > $tss"
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


# for transcripts
echo "transcripts"
base=${prefix%%_*}
prefix=${base}_transcripts
bed=$prefix.bed
$dir/../../util/gtfTobed.sh $db transcript > $bed
# get TSS
tss=$prefix.TSS.bed
python $dir/../../util/getTss.sh $bed > $tss
# get regions around TSS
for width in 400 2000 5000 10000;do
  side=$(($width/2))
  region=$prefix.TSS.${width}bp.bed
  com="bedtools slop -b $side -i $tss -g $gs|\
    sort -k1,1 -k2,2n > $region"
  echo $com
  eval $com
done
