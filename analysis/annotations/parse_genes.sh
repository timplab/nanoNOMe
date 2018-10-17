#!/bin/bash
# parse genes
echo "parsing genes"
db=$(readlink -f ../../annotations/hg38/hg38_genes.gtf.gz)
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


# for transcripts
base=${prefix%%_*}
prefix=${base}_transcripts
bed=$prefix.bed
[ -e $bed ]||\
  ../../util/gtfTobed.sh $db transcript > $bed
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
