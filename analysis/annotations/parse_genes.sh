#!/bin/bash
dir=$(dirname $(readlink -f $0))
# parse genes
echo "parsing genes"
db="$1"
prefix=${db%%.*}

# convert to bed
bed=$prefix.bed
com="$dir/../../util/gtfTobed.sh $db > $bed"
echo $com
eval $com

# get TSS
tss=$prefix.TSS.bed
com="python $dir/../../util/bed_parser.py getstart -b $bed -o $tss"
echo $com
eval $com

# get regions around TSS
for width in 400 2000;do
  side=$(($width/2))
  region=$prefix.TSS.${width}bp.bed
  com="python $dir/../../util/bed_parser.py region -u $side -d $side \
    -b $tss | sort -k1,1 -k2,2n > $region"
  echo $com
  eval $com
done


# for transcripts
base=${prefix%%_*}
prefix=${base}_transcripts
bed=$prefix.bed
[ -e $bed ]||\
  $dir/../../util/gtfTobed.sh $db transcript > $bed
# get TSS
tss=$prefix.TSS.bed
[ -e $tss ]||\
  python $dir/../../util/bed_parser.py getstart -b $bed -o $tss
# get regions around TSS
for width in 400 2000;do
  side=$(($width/2))
  region=$prefix.TSS.${width}bp.bed
  [ -e $region ]||\
    python $dir/../../util/bed_parser.py region -u $side -d $side \
    -b $tss | sort -k1,1 -k2,2n > $region
done
