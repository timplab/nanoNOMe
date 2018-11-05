#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
plotdir=$root/plots/readlevel
methdir=$root/pooled/methylation/methbyread_all
reg=bcan_promoters
[ -z $2 ]||reg="$2"
annodir=$root/annotations/breastcancer
if [ "$reg" == "bcan_promoters" ];then
  bed=$annodir/bcan_genes.TSS.5000bp.bed
elif [ "$reg" == "comparison" ];then
  bed=$annodir/bcan_10a_vs_231_promoters.bed
elif [ "$reg" == "rna" ];then
  bed=$annodir/bcan_rnaseq_methylation.CGI.TSS.5000bp.bed
fi

for cell in MCF10A MCF7 MDAMB231;do
  echo $cell
  for mod in cpg gpc;do
    base=$cell.$mod
    meth=$methdir/$base.pooled.meth.bed.gz
    plotpath=$plotdir/$base.$reg.heatmap.pdf
    com="python -u $srcdir/../../script/readlevelHeatmap.py -v -w 20 -c 5 \
      -i $meth -r $bed -o $plotpath"
    echo $com
    eval $com
  done
done

