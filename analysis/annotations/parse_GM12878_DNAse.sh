#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
dir="$1"
bed="$2"
prefix=$dir/GM12878_DNAse_top

# get top 10k peaks
sort -k7,7nr $bed | head -n10000 | sort -k1,1 -k2,2n > $prefix.bed

# get center
center=$prefix.center.bed
com="python $srcdir/../../util/bed_parser.py getcenter -v -b $bed |\
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

