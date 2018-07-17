#!/bin/bash
cell=GM12878
cell=GM12878_illumina
cell=GM12878_illumina_negcontrol
cell=GM12878_BSseq_ENCLB794YYH
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
freqdir=$root/pooled/methylation/mfreq_all
plotdir=$root/plots/aggregate
samp=gm12878
echo $cell
cpg=$freqdir/$cell.cpg.methfreq.txt.gz
gpc=$freqdir/$cell.gpc.methfreq.txt.gz
if [ "$1" == "ctcf" ];then
  #dbpath="/home/isac/Dropbox/Data/nome-seq/db/ctcf/gm12878/GM12878_CTCF.ENCFF536RGD.2kb.bed"
  dbpath="/home/isac/Dropbox/Data/nome-seq/db/$samp/ctcf/GM12878_CTCF.2kb.bed"
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
