#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
dir="$1"
ctcf="$2"
prefix=$dir/GM12878_CTCF_hg38

gmctcf=$prefix.bed
# first get gm12878-specific ctcf binding sites
com="bedtools intersect -a $ctcf -b $dir/GM12878_CTCF_ChIP.bed -u |\
  sort -k1,1 -k2,2n > $gmctcf"
echo $com
eval $com

# get center
center=$prefix.center.bed
com="python $srcdir/../../util/bed_parser.py getcenter -v -b $gmctcf |\
  sort -k1,1 -k2,2n > $center"
echo $com
eval $com

# windows
for width in 400 2000; do
  side=$(($width/2))
  region=$prefix.center.${width}bp.bed
  com="python $srcdir/../../util/bed_parser.py region -u $side -d $side \
    -b $center | sort -k1,1 -k2,2n > $region"
  echo $com
  eval $com
done

