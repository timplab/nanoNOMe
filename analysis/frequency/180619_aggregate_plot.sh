#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
cells="GM12878 GM12878_scNOMe GM12878_scNOMe_negcontrol GM12878_BSseq_1 GM12878_wgs"
#cells="GM12878 GM12878_scNOMe GM12878_scNOMe_negcontrol GM12878_BSseq_1"
#cells="GM12878_wgs"
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -n $1 ]||root="$1"
freqdir=$root/pooled/methylation/mfreq_all
plotdir=$root/plots/aggregate
plotter=$srcdir/../../script/nanonome_plots.R
dbpath="$root/database/gm12878/ctcf/GM12878_CTCF.2kb.bed" #default CTCF
if [ "$2" == "ctcf" ];then
  echo "CTCF"
elif [ "$2" == "tss" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38/hg38_genes.TSS.2000bp.bed"
elif [ "$2" == "cgitss" ];then
  dbpath=/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.2kb.CGI.bed
elif [ "$2" == "dnase" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/dnase/GM12878_dnase.2kbregion.bed"
fi


for cell in $cells;do
  cpg=$freqdir/$cell.cpg.methfreq.txt.gz
  gpc=$freqdir/$cell.gpc.methfreq.txt.gz
  plotpath=$plotdir/$cell.$1.aggregate.pdf

  Rscript $plotter aggregateByDistance -c $cpg -g $gpc -r $dbpath -o $plotpath 
done
