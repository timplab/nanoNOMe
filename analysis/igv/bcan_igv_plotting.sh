#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
outdir=$root/igv
[ -e $outdir ]||mkdir $outdir
annodir=$root/annotations/breastcancer
reg="bcanpromoter"
[ -z $2 ]||reg="$2"
if [ "$reg" == "bcanpromoter" ];then
  bed=$annodir/bcan_10a_vs_231_promoters.bed
elif [ "$reg" == "bcansv" ];then
  bed=$annodir/sv_bcangenes.bed
fi

cells="MCF10A MCF7 MDAMB231"

for cell in $cells;do
  echo $cell
  bam=$root/pooled/bam/$cell.pooled.bam
  mbeddir=$root/pooled/methylation/methbyread_all
  cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
  script=$srcdir/../../script/convertBam.py
  log=$outdir/${cell}_$reg.log
  out=$outdir/${cell}_$reg.bam

  com="python -u $script -v -t 10 \
    -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
    samtools sort -o $out"
  echo $com
  eval $com
  samtools index $out
done


