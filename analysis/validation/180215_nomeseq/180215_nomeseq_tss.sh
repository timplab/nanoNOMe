#!/bin/bash
root=/dilithium/Data/NGS/Aligned/180215_nomeseq
beddir=$root/bed
tssdir=$root/tss
[ -e $tssdir ]||mkdir $tssdir
tss="/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.bed"
plotdir="/home/isac/Dropbox/Data/Genetics/MethSeq/180215_nomeseq/plots"
plotscript=/home/isac/Code/ilee/methylation/scripts/methylation_plots.R

for bedgraph in `find $beddir -name "*bedGraph.gz"`;do
  base=$(basename $bedgraph)
  base=${base%%.*}
  echo $base
  if [ "$1" == "getdistance" ];then
    gunzip -c $bedgraph |\
      bedtools closest -D b -b $tss -a stdin |\
      awk '{ OFS="\t" }{ if(($NF<5000)&&($NF>-5000)) print }' >  $tssdir/$base.gpc.distance2tss.bed
  fi
  if [ "$1" == "parsedistance" ];then
    awk '{ OFS="\t" }{ print $NF,$5,$6 }' $tssdir/$base.gpc.distance2tss.bed > $tssdir/$base.gpc.dist2tss.txt
  fi 
  if [ "$1" == "plotmeth" ];then
    [ -e $plotdir ]||mkdir $plotdir
    Rscript $plotscript methByDistance -i $tssdir/$base.gpc.dist2tss.txt -o $plotdir/$base.tssmeth.pdf
  fi
done

