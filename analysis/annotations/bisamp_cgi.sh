#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38
bismap=$root/bismap/bismapM50_chr22_part.bed
cgi=$root/hg38_cgi.txt.gz
cgibed=$root/hg38_cgi.bed
[ -e $cgibed ]||\
  gunzip -c $cgi |\
  awk 'OFS="\t"{ print $2,$3,$4,".",".","." }' |\
  grep -v _ \
  > $cgibed

bedtools intersect -wao -b $bismap -a $cgibed |\
  awk '{ print $0,$3-$2 }' |\
  grep chr22 |\
  awk '{ if ($2<10933654) print }'
