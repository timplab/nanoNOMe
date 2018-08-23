#!/bin/bash
cell=GM12878
cellname=$(echo $cell | tr "[A-Z]" "[a-z]")
root=/shared/data
rawdir=$root/raw/$cellname/fast5
outdir=$root/tar
outzip=$outdir/${cell}_nanoNOMe.zip
batch=../../util/batchcommand.scr

log=$outzip.log
com="python -u ../../util/fast5zip_accessor.py zip -v -z $outzip -d $rawdir > $outzip.filelist.txt"
echo $com
sbatch -e $log -o $log -t 200:0:0 $batch "$com"
