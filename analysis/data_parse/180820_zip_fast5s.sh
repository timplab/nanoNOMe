#!/bin/bash
cell=GM12878
cellname=$(echo $cell | tr "[A-Z]" "[a-z]")
root=/shared/data
rawdir=$root/raw/$cellname/fast5
outdir=$root/tar
outzip=$outdir/${cell}_nanoNOMe.zip
batch=../../script/batchcommand.scr

log=$outzip.log
com="cd $rawdir && zip -vr $outzip ."
echo $com
sbatch -e $log -o $log -t 48:0:0 $batch "$com"
