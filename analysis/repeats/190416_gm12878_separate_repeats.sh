#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data

nomepre=$root/nanonome/pooled/mfreq/GM12878_nanoNOMe.pooled
bspre=$root/bsseq/mfreq/GM12878_BSseq
prefixes="$bspre $nomepre"
outdir=$root/subset/repeats/mfreq
[ -e $outdir ]||mkdir $outdir

bed=$root/hg38/hg38_repeats.bed

# first subset bed
patterns="LINE SINE LTR DNA"
if [ "$1" == "reg" ];then
  for pattern in $patterns; do 
    reg=$root/hg38/hg38_repeats.$pattern.bed
    if [ ! -e $reg ];then
      awk -v pat=$pattern 'OFS="\t"{ if ($7==pat) print }' $bed |\
        sort -k1,1 -k2,2n > $reg
    fi
  done
fi

if [ "$1" == "data" ];then
  scr="/home/isac/Code/nanoNOMe/scripts/methylFreq_destrand.py"
  mods="gpc cpg"
#  mods="gpc"
  patterns="DNA"
  for pattern in $patterns; do
    reg=$root/hg38/hg38_repeats.$pattern.bed
    for mod in $mods; do
      for pre in $prefixes; do
        fp=$pre.$mod.mfreq.txt.gz
        base=$(basename "${fp%%.mfreq*}")
        echo $base.$pattern
        out=$outdir/$base.$pattern.mfreq.txt.gz
        args="$args $fp $reg $mod $out"
      done
    done
  done
  com="tabix {1} -T {2} |\
    $scr -m {3} | bgzip > {4}"
  parallel -N 4 $com ::: $args
fi

