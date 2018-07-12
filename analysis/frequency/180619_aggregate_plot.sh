#!/bin/bash
cell=GM12878
cell=GM12878_sample
#cell=GM12878_control
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
if [ "$cell" == "GM12878" ];then
  freqdir=$root/pooled/methylation/mfreq
  samp=gm12878
elif [ "$cell" == "GM12878_sample" ]||
  [ "$cell" == "GM12878_control" ];then #scNOMe sample
  freqdir=$root/validation/scNOMe/methfreq
  samp=gm12878
fi
echo $cell
cpg=$freqdir/$cell.cpg.methfreq.txt.gz
gpc=$freqdir/$cell.gpc.methfreq.txt.gz
plotdir=$root/plots/aggregate
[ -e $plotdir ]||mkdir -p $plotdir
dbroot=$root/database
if [ "$1" == "ctcf" ];then
  #dbpath="/home/isac/Dropbox/Data/nome-seq/db/ctcf/gm12878/GM12878_CTCF.ENCFF536RGD.2kb.bed"
  dbpath="$dbroot/$samp/ctcf/GM12878_ctcf.2000bp.bed"
elif [ "$1" == "tss" ];then
  dbpath=/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.2kb.bed
elif [ "$1" == "cgitss" ];then
  dbpath=/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.2kb.CGI.bed
elif [ "$1" == "dnase" ];then
  dbpath="/home/isac/Dropbox/Data/nome-seq/db/gm12878/dnase/GM12878_dnase.2kbregion.bed"
elif [ "$1" == "atac" ];then
  dbpath="/home/isac/Dropbox/Data/nome-seq/db/gm12878/atac/GM12878_atac.2kbregion.bed"
fi

plotpath=$plotdir/$cell.$1.aggregate.pdf
heatpath=$plotdir/$cell.$1.heatmap.pdf
plotter=~/Code/ilee/plot/nanonome_plots.R

Rscript $plotter aggregateByDistance -c $cpg -g $gpc -r $dbpath -o $plotpath 
#Rscript $plotter heatmapByDistance -c $cpg -g $gpc -r $dbpath -o $heatpath 
