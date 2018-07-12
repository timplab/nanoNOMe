#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
pooldir=$root/pooled/fastq
for cell in MDAMB231;do
  echo $cell
  fqdir=$root/${cell,,}
  fqs=`find $fqdir -name "*$cell*fastq.gz"`
  echo $fqs > $root/log/$cell.fastqmerge.log
  fq=$pooldir/$cell.pooled.fastq.gz
  cat $fqs > $fq

done
