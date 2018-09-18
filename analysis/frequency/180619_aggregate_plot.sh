#!/bin/bash
cells="GM12878 GM12878_scNOMe GM12878_scNOMe_negcontrol GM12878_BSseq_1 GM12878_wgs"
#cells="GM12878 GM12878_scNOMe GM12878_scNOMe_negcontrol GM12878_BSseq_1"
#cells="GM12878_wgs"
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
freqdir=$root/pooled/methylation/mfreq_all
plotdir=$root/plots/aggregate
plotter=../../script/nanonome_plots.R
if [ "$1" == "ctcf" ];then
  #dbpath="/home/isac/Dropbox/Data/nome-seq/db/ctcf/gm12878/GM12878_CTCF.ENCFF536RGD.2kb.bed"
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/ctcf/GM12878_CTCF.2kb.bed"
elif [ "$1" == "tss" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38/hg38_genes.TSS.2000bp.bed"
elif [ "$1" == "cgitss" ];then
  dbpath=/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.2kb.CGI.bed
elif [ "$1" == "dnase" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/dnase/GM12878_dnase.2kbregion.bed"
elif [ "$1" == "atac" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/atac/GM12878_atac.2kbregion.bed"
fi


for cell in $cells;do
  cpg=$freqdir/$cell.cpg.methfreq.txt.gz
  gpc=$freqdir/$cell.gpc.methfreq.txt.gz
  plotpath=$plotdir/$cell.$1.aggregate.pdf
  heatpath=$plotdir/$cell.$1.heatmap.pdf

  Rscript $plotter aggregateByDistance -c $cpg -g $gpc -r $dbpath -o $plotpath 
  #Rscript $plotter heatmapByDistance -c $cpg -g $gpc -r $dbpath -o $heatpath 
done
