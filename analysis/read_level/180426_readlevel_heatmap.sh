#!/bin/bash
codedir=`readlink -f ../../`
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
scriptdir=~/Code/ilee/nanopolish/script
outdir=$root/intersect
[ -e $outdir ]||mkdir -p $outdir
db=/home/isac/Dropbox/Data/nome-seq/db/breastcancer/bcanpoorprog.4kb.bed
dbname=$(basename "$db")
dbname=${dbname%.bed}

for cell in MCF10A MCF7 MDAMB231;do
  echo $cell
  intersect=$outdir/$cell.gpc.$dbname.bedGraph
  if [ "$1" == "intersect" ];then
    beddir=`readlink -f $root/${cell,,}/ngmlr/pooled`
    bed=`find $beddir -name "$cell.gpc.pooled.meth.bedGraph.gz"`
    gunzip -c $bed |\
      bedtools intersect -b $db -a stdin -F 1 -sorted -wo \
      > $intersect
  fi
  plotpath=~/Dropbox/Data/nome-seq/plots/heatmap/$cell.$dbname.heatmap.pdf
  if [ "$1" == "heatmap" ];then
    python $codedir/../../../nomeseq/script/readlevelHeatmap.py -w 20 \
      -g $intersect -o $plotpath
  fi
done

