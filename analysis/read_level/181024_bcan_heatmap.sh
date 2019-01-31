#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/kyber/Data/Nanopore/projects/nanonome/analysis
plotdir=$root/plots/heatmap
[ -e $plotdir ]||mkdir $plotdir
methdir=$root/data/nanonome/pooled/mbed
reg=bcan_promoters
[ -z $1 ]||reg="$1"
annodir=$root/data/bcan
if [ "$reg" == "bcan_promoters" ];then
  bed=$annodir/bcan_select.2kb.bed
elif [ "$reg" == "comparison" ];then
  bed=$annodir/bcan_10a_vs_231_promoters.bed
elif [ "$reg" == "rna" ];then
  bed=$annodir/bcan_rnaseq_methylation.CGI.TSS.5000bp.bed
fi

for cell in MCF10A MCF7 MDAMB231;do
  echo $cell
  for mod in cpg gpc;do
    prefix=${cell}_nanoNOMe.pooled.$mod
    pres="$pres $prefix"
  done
done

meth=$methdir/{}.meth.bed.gz
plotpath=$plotdir/{}.$reg.heatmap.pdf
log=$plotdir/{}.$reg.heatmap.log
com="python -u $srcdir/../../script/readlevelHeatmap_diagonal.py -v \
  -i $meth -r $bed -o $plotpath 2> $log"

parallel "$com" ::: $pres
