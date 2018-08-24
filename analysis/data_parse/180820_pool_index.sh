#!/bin/bash
root=/shared/data/analysis
outdir=$root/pooled/fastq
zipdir=/shared/data/zip
[ -e $zipdir ]||mkdir $zipdir
rawroot=/shared/data/raw
np=/shared/Code/nanopolish/nanopolish
batchcom=../../util/batchcommand.scr
cells="MCF10A"

for cell in $cells;do
  sum=$outdir/$cell.pooled.summary.txt
  celldir=$(echo $cell | tr "[A-Z]" "[a-z]")
  echo $cell
  if [ "$1" == "summary" ];then
    dir=$root/$celldir
    sums=$(find $dir -name "*summary.txt")
    awk 'FNR==1 && NR!=1 {next}{print}' $sums > $sum
  fi
  rawdir=$rawroot/$celldir
  fq=$outdir/$cell.pooled.fastq.gz
  if [ "$1" == "index" ];then
    com="$np index -v -d $rawdir -s $sum $fq"
    echo $com
    log=$outdir/$cell.index.log
    sbatch -e $log -o $log -c 1 -J "${cell}-index" $batchcom $com
  fi
  outzip=$zipdir/${cell}_nanoNOMe.zip
  if [ "$1" == "zip" ];then
    log=$zipdir/$cell.zip.log
    com="python -u ../../util/fast5zip_accessor.py zip -v -z $outzip -d $rawdir"
    echo $com
    sbatch -e $log -o $log -J "${cell}-zip" -t 200:0:0 $batchcom $com
  fi
  if [ "$1" == "zipindex" ];then
    readdb=$fq.index.readdb
    ziplog=$zipdir/$cell.zip.log
    bed=$root/pooled/bed/$cell.pooled.bed
    log=$zipdir/$cell.zipindex.log
    index=$zipdir/${cell}_nanoNOMe.zipindex.bed.gz
    com="python -u ../../util/fast5zip_accessor.py index -v -l $ziplog -r $readdb -b $bed | bgzip > $index; tabix -p bed $index"
    echo $com
    sbatch -e $log -o $log -J "${cell}-index" -t 200:0:0 $batchcom $com
  fi
done
