#!/bin/bash
echo "MNase-seq metaplots"
srcdir=$(dirname $(readlink -f $0))
root="$1"
dbdir=$root/data/gm12878
dat=$dbdir/GM12878_MNase_hg38.bedGraph.gz
plotdir=$root/plots/metaplots

for reg in CTCF DNAse;do
  echo $reg
  if [ "$reg" == "CTCF" ];then
    db=$root/annotations/gm12878/GM12878_CTCF_hg38.center.bed
  elif [ "$reg" == "DNAse" ];then
    db=$root/annotations/gm12878/GM12878_DNAse_top.center.bed
  fi
  dist=$dbdir/GM12878_MNase_hg38.$reg.bedGraph
  if [ ! -e "$dist" ];then
    echo "dist"
    gunzip -c $dat | bedtools closest -D b -b $db -a stdin |\
      awk '{ if(sqrt($NF*$NF) <= 1000) print }' > $dist
  fi
  echo "plot"
  plotpath=$plotdir/GM12878_MNase.$reg.metaplot.pdf
  script=$srcdir/../../script/mnase_plot.R
  Rscript $script metaplotByDistance -i $dist -o $plotpath 2> /dev/null
done
