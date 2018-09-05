#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
bamdir=$root/pooled/bam
mbeddir=$root/pooled/methylation/methbyread_all
svdir=$root/pooled/sv
parser=../../nanopolish/SVmethylation.py
bammer=../../script/parse_sniffles.py
cells="MCF10A MCF7 MDAMB231"
cells="MCF10A"

for cell in $cells;do
  echo $cell
  bam=$bamdir/$cell.pooled.bam
  cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
  sv=$svdir/$cell.sniffles.vcf
  out=$svdir/$cell.TRAmultiregion.bed
  outbam=$svdir/$cell.TRAmultiregion.bam

#  tail $sv -n 10 | python $bammer bam -v -t 8 -b $bam |\
#    awk '{if ($1 != "@PG" || n < 1) print; if ($1 == "@PG") n +=1 }' |\
#    samtools sort -o $outbam
#  samtools index $outbam
  tail $sv -n 10 | python $parser -v -t 8 -b $bam -c $cpg -g $gpc 
done
