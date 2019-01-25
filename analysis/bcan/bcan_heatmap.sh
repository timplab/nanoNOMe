#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
methdir=$root/pooled/methylation/methbyread_all
outdir=$root/plots/readlevel
[ -e $outdir ]||mkdir $outdir
annodir=$root/annotations/breastcancer
reg="bcanselect"
[ -z $2 ]||reg="$2"
if [ "$reg" == "10AvsMCF7" ];then
  bed=$annodir/MCF10A_vs_MCF7_top_epigenetic_state_genes.TSS.2000bp.bed
elif [ "$reg" == "10AvsMDAMB231" ];then
  bed=$annodir/MCF10A_vs_MDAMB231_top_epigenetic_state_genes.TSS.2000bp.bed
elif [ "$reg" == "bcangenes" ];then
  bed=$annodir/breastcancer_genes.TSS.4kb.bed
elif [ "$reg" == "bcanbothcgi" ];then
  bed=$annodir/bcan_diffexp_bothcancer_cgi.TSS.2kb.bed
elif [ "$reg" == "bcanselect" ];then
  bed=$annodir/bcan_select.2kb.bed
elif [ "$reg" == "bcandel" ];then
  bed=$root/sv/bcan_hom_DEL.flank.2000bp.bed
fi

plotter="$srcdir/../../script/readlevelHeatmap_diagonal.py"
cells="MCF10A $reg"
cells="MCF10A MCF7 MDAMB231"

for cell in $cells;do
  echo $cell
  for mod in cpg gpc;do
    prefix=$cell.$mod
    pres="$pres $prefix"
  done
done
echo $pres

meth=$methdir/{}.pooled.meth.bed.gz
plotpath=$outdir/{}_readlevel_heatmap_${reg}.pdf
log=$outdir/{}_readlevel_heatmap_${reg}.log
com="python -u $plotter -v \
  -i $meth -r $bed -o $plotpath 2> $log"
echo $com
parallel "$com" ::: $pres
