#!/bin/bash
cells="MCF10A MCF7 MDAMB231"
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
bed=$root/plots/readlevel/enhancerwindow.bed
reg=testset
pooldir=$root/pooled
for cell in $cells;do
  bam=$pooldir/bam/$cell.pooled.bam
  mbeddir=$pooldir/methylation/methbyread_all
  cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
  script=../../nanopolish/convertBam.py
  log=$pooldir/igv/$cell.$reg.log
  out=$pooldir/igv/$cell.$reg.bam

  python -u $script -v -t 10 -b $bam -c $cpg -g $gpc -r $bed -o $out 2> $log
  outsort=$pooldir/igv/$cell.$reg.sorted.bam
  samtools sort -o $outsort $out
  samtools index $outsort
done
