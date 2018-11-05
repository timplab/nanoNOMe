#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
plotdir=$root/plots/readlevel
methdir=$root/pooled/methylation/methbyread_all
reg=ctcf
[ -z $2 ]||reg="$2"
annodir=$root/annotations/gm12878
if [ "$reg" == "ctcf" ];then
  bed=$annodir/GM12878_CTCF_hg38.cgi.center.2000bp.bed
fi

for cell in GM12878;do
  echo $cell
  for mod in cpg gpc;do
    prefix=$cell.$mod
    pres="$pres $prefix"
  done
done
echo $pres

meth=$methdir/{}.pooled.meth.bed.gz
plotpath=$plotdir/{}.$reg.heatmap.pdf
com="python -u $srcdir/../../script/readlevelHeatmap.py -v -w 20 -c 5 \
  -i $meth -r $bed -o $plotpath"
echo $com
parallel "$com" ::: $pres
