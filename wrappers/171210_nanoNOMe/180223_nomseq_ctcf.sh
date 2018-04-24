#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171210_nomeseq
beddir=$root/bed
ctcfdir=$root/ctcf
[ -e $ctcfdir ]||mkdir $ctcfdir
ctcf="/home/isac/Dropbox/Data/nome-seq/db/ctcf/ctcf_hg38.sorted.bed"
plotdir="/home/isac/Dropbox/Data/Nanopore/171210_nomeseq/plots"
plotscript=/home/isac/Code/ilee/methylation/scripts/methylation_plots.R

for bedgraph in `find $beddir -name "*gpc.*freq*sorted*bedGraph.gz"`;do
  base=$(basename $bedgraph)
  base=${base%%.*}
  echo $base
  if [ "$1" == "getdistance" ];then
    gunzip -c $bedgraph |\
      bedtools closest -D b -b $ctcf -a stdin |\
      awk '{ OFS="t" }{ if(($NF<5000)&&($NF>-5000)) print }' >  $ctcfdir/$base.gpc.dist2CTCF.bed
  fi
  if [ "$1" == "parsedistance" ];then
    awk '{ OFS="\t" }{ print $NF,$4,$5 }' $ctcfdir/$base.gpc.dist2CTCF.bed > $ctcfdir/$base.gpc.dist2CTCF.txt
  fi 
  if [ "$1" == "plotmeth" ];then
    [ -e $plotdir ]||mkdir $plotdir
    Rscript $plotscript methByDistance -i $ctcfdir/$base.gpc.dist2CTCF.txt -o $plotdir/$base.ctcfmeth.pdf
  fi
done

