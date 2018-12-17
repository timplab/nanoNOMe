#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
methdir=$root/pooled/methylation/methbyread_all
outdir=$root/plots/readlevel
[ -e $outdir ]||mkdir $outdir
annodir=$root/annotations/breastcancer
reg="bcangenes"
[ -z $2 ]||reg="$2"
if [ "$reg" == "10AvsMCF7" ];then
  bed=$annodir/MCF10A_vs_MCF7_top_epigenetic_state_genes.TSS.2000bp.bed
elif [ "$reg" == "10AvsMDAMB231" ];then
  bed=$annodir/MCF10A_vs_MDAMB231_top_epigenetic_state_genes.TSS.2000bp.bed
elif [ "$reg" == "bcangenes" ];then
  bed=$annodir/breastcancer_genes.TSS.4kb.bed
fi

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
plotpath=$outdir/{}_readlevel_heatmap_expression_${reg}.pdf
log=$outdir/{}_readlevel_heatmap_expression_${reg}.log
com="python -u $srcdir/../../script/readlevelHeatmap.py -v \
  -i $meth -r $bed -o $plotpath 2> $log"
echo $com
parallel "$com" ::: $pres
