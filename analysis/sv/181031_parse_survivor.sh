#!/bin/bash
dir=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/sv
vcf=$dir/merged_SURVIVOR_1kbp_typesafe_10_31_2018.sort.anno.vcf.gz
out=$dir/exclusive_sv.vcf
header=$dir/sniffles_header.txt

gunzip -c $vcf | grep "#" > $header
cat $header > $out
gunzip -c $vcf | grep -E "001000|101000|010111|110111" >> $out

for svtype in DEL DUP TRA INS;do
  sub=$dir/exclusive_$svtype.vcf
  cat $header > $sub
  grep "<$svtype>" $out >> $sub
done
 
sub=$dir/exclusive_hom.vcf
cat $header > $sub
grep "1/1" $out >> $sub

