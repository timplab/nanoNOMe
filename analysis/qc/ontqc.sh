#!/bin/bash
sumdir=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/fastq
sums=$(find $sumdir -name "*summary.txt")
for sum in $sums; do
  base=${sum%%.*}
  out=$base.ontqc.txt
  python ../../util/ontQC.py $sum > $out
done
