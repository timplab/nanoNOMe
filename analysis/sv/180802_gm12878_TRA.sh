#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
bamdir=$root/pooled/bam
svdir=$root/pooled/sv
cell=GM12878
parser=../../util/parse_sniffles.py

#bam=$bamdir/$cell.pooled.bam
sv=$svdir/$cell.sniffles.vcf
out=$svdir/$cell.multiregion.bed

python $parser -i $sv multiregion |\
  sort -k1,1 -k2,2n > $out
