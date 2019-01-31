#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/kyber/Data/Nanopore/projects/nanonome/analysis
outdir=$root/igv
[ -e $outdir ]||mkdir $outdir
annodir=$root/data/bcan
reg="bcanpromoter"
[ -z $1 ]||reg="$1"
if [ "$reg" == "bcanpromoter" ];then
  bed=$annodir/bcan_diffexp_bothcancer_cgi.TSS.20kb.bed
elif [ "$reg" == "bcansv" ];then
  bed=$annodir/sv_bcangenes.bed
fi

cells="MCF10A MCF7 MDAMB231"
script=$srcdir/../../script/convertBam.py
mbeddir=$root/data/nanonome/pooled/mbed

for cell in $cells;do
  echo $cell
  base=${cell}_nanoNOMe.pooled
  bam=$root/data/nanonome/pooled/bam/$base.bam
  cpg=$mbeddir/$base.cpg.meth.bed.gz
  gpc=$mbeddir/$base.gpc.meth.bed.gz
  log=$outdir/$base.$reg.log
  out=$outdir/$base.$reg.bam

  com="python -u $script -v -t 10 \
    -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
    samtools sort -o $out"
  echo $com
  eval $com
  samtools index $out
done


