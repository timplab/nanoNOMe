#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/atac
listpath=./atacseq_sra_runtable.txt

logdir=$root/log
fqdir=$root/fastq
[ -e $fqdir ]||mkdir $fqdir
log=$logdir/gm12878_atacseq_ascp.log
[ -e $log ]&&rm $log

sras=$(awk '{ if (NR>1) print $5 }' $listpath | tr "\n" " " )
echo $sras
cd $fqdir 
parallel fastq-dump --split-3 --gzip {} ::: $sras
