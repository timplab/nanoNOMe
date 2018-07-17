#!/bin/bash
# parse genes
# path to the ensembl gtf.gz provided as an arugment
echo "parsing genes"
db="$1"
prefix=${db%%.*}

# convert to bed
bed=$prefix.bed
[ -e $bed ]||\
  ../../util/gtfTobed.sh $db > $bed

# get TSS
tss=$prefix.TSS.bed
[ -e $tss ]||\
  python ../../util/bed_parser.py getstart -b $bed -o $tss

# get regions around TSS
for width in 400 2000;do
  side=$(($width/2))
  region=$prefix.TSS.${width}bp.bed
  [ -e $region ]||\
    python ../../util/bed_parser.py region -u $side -d $side \
    -b $tss | sort -k1,1 -k2,2n > $region
done

