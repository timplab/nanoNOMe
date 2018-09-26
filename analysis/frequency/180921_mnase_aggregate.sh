#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
dat=$root/database/gm12878/mnase/GM12878_mnase.hg38.bedGraph.gz
outdir=$root/distancebed
[ -e $outdir ]||mkdir $outdir
plotdir=$root/plots/aggregate
if [ "$1" == "ctcf" ];then
  #dbpath="/home/isac/Dropbox/Data/nome-seq/db/ctcf/gm12878/GM12878_CTCF.ENCFF536RGD.2kb.bed"
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/ctcf/GM12878_ctcf.center.bed"
elif [ "$1" == "tss" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38/hg38_genes.TSS.2000bp.bed"
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38/hg38_genes.TSS.sorted.bed"
elif [ "$1" == "cgitss" ];then
  dbpath=/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.2kb.CGI.bed
elif [ "$1" == "dnase" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/dnase/GM12878_dnase.2kbregion.bed"
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/dnase/GM12878_dnase.center.bed"
elif [ "$1" == "atac" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/atac/GM12878_atac.2kbregion.bed"
elif [ "$1" == "shuffle" ];then
  dbpath="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/ctcf/GM12878_CTCF.center.shuffle.bed"
fi

dist=$outdir/GM12878_mnase_distance_$1.bedGraph
if [ ! -e $dist ];then
  echo "$dist"
  gunzip -c $dat | bedtools closest -D b -b $dbpath -a stdin |\
    awk '{ if(sqrt($NF*$NF) <= 1000) print }' > $dist
fi

#script=../../script/parseDNAse.py
#python $script by-distance -v -i $dist

plotpath=$plotdir/GM12878_mnase_aggregate_$1.pdf
script=../../script/dnase_plot.R
Rscript $script aggregateByDistance -i $dist -o $plotpath
