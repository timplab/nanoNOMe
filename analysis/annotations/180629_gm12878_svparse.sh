#!/bin/bash
svdir="/dilithium/Data/Nanopore/projects/nomeseq/analysis/gm12878/ngmlr/sniffles"
vcf="$svdir/GM12878.sniffles.vcf"
parser=/home/isac/Code/ilee/sv/parse_sniffles.py

if [ "$1" == "getcoords" ];then
  svtype=$(echo "$2" | tr a-z A-Z)
  vcfpre=${vcf%.vcf}
  subsetout=$vcfpre.$2.vcf
  outbed=$vcfpre.$2.bed
  grep "<$svtype>" $vcf |\
    python $parser svtype coords |\
    uniq > $outbed
fi

if [ "$1" == "getregion" ];then
  parser="/home/isac/Code/ilee/annotation/script/bed_parser.py"
  vcfpre=${vcf%.vcf}
  outbed=$vcfpre.$2.bed
  regbed=$vcfpre.$2.400b.region.bed
  python $parser region -u 200 -d 200 -b $outbed -o $regbed

#  awk 'OFS="\t"{ $7=$2;$8=$3; if ($2<1000) $2=0; else $2=$2-1000;$3=$3+1000; print }' $outbed > $regbed
fi
#cut -f1,2,5,10,11 $vcfout
