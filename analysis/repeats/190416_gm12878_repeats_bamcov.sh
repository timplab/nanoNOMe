#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data

nomepre=$root/nanonome/pooled/bed/GM12878_nanoNOMe.pooled
bspre=$root/bsseq/GM12878_BSseq.nodup
prefixes="$bspre $nomepre"
outdir=$root/subset/repeats/bed
[ -e $outdir ]||mkdir $outdir

bed=$root/hg38/hg38_repeats.bed

patterns="LINE SINE LTR DNA"
patterns="DNA"

if [ "$1" == "intersect" ];then
  for pre in $prefixes; do
    base=$(basename "$pre")
#    for pattern in $patterns; do 
#      reg=$root/hg38/hg38_repeats.$pattern.bed
#      out=$outdir/$base.$pattern.cov.bed
#      args="$args $reg $pre.bed $out"
#    done
    args="$args $bed $pre.bed $outdir/$base.repeats.bed"
  done
  db="{1}"
  infp="{2}"
  out="{3}"
  com="bedtools intersect -sorted -wb -a $infp -b $db |\
    awk 'OFS=\"\t\"{ if(\$5>=20){print}}' > $out" # filter by mapq >=20 b/c nanopore reads are
#  com="bedtools coverage -hist -sorted -a $db -b $infp > $out"
#  com="echo $db $infp $out"
  parallel -N 3 "$com" ::: $args
fi

if [ "$1" == "cov" ];then
  for pre in $prefixes; do
    base=$(basename "$pre")
    samp=$outdir/$base.repeats.bed
    out=$outdir/$base.repeats.cov.bed
    args="$args $bed $pre.bed $out"
  done
  db="{1}"
  infp="{2}"
  out="{3}"
  com="bedtools coverage -hist -sorted -a $db -b $infp > $out"
  parallel -N 3 "$com" ::: $args
fi

if [ "$1" == "mapq" ];then
  for pre in $prefixes; do
    base=$(basename "$pre")
    echo $base
    samp=$outdir/$base.repeats.bed
    awk '{ print $NF,$5 }' | sort -k1,1 > $out
  done
fi


