#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
annodir=$root/annotations/CTCF
[ -e $annodir ]||mkdir $annodir

all=$annodir/CTCFBSDB_allcomp.txt.gz
if [ "$1" == "dl" ];then
  wget http://insulatordb.uthsc.edu/download/allcomp.txt.gz \
    -O $all
fi
bed18=$annodir/CTCFBSDB_allcomp_hg18.bed
if [ "$1" == "bed" ];then
  gunzip -c $all | awk '{OFS="\t"}$2=="Human"{ $2=$2-1; print $1,$2,$3,".",".","." }' | cut -f3 |\
    tr ":" "\t" | tr "-" "\t" |\
    awk '{OFS="\t"}NR>1{ $2=$2-1; print $1,$2,$3,".",".","." }' > $bed18
fi
bed38=$annodir/CTCFBSDB_allcomp_hg38.bed
center=$annodir/CTCFBSDB_allcomp_hg38.center.bed
nogene=$annodir/CTCFBSDB_allcomp_hg38.center.noTSS.bed
nogene2kb=$annodir/CTCFBSDB_allcomp_hg38.center.noTSS.2000bp.bed
genes5kb=$root/annotations/hg38/hg38_genes.TSS.5000bp.bed
gs=$root/annotations/hg38/hg38_genomesize.txt
if [ "$1" == "liftover" ];then
  chain18=$annodir/../hg38/hg18ToHg38.over.chain
  liftOver $bed18 $chain18 $bed38 $bed38.unmapped
  python ../../util/bed_parser.py getcenter -v -b $bed38 |\
    sort -k1,1 -k2,2n > $center
  bedtools intersect -v -a $center -b $genes5kb > $nogene
  bedtools slop -b 1000 -i $nogene -g $gs |\
    sort -k1,1 -k2,2n > $nogene2kb
fi

gmbed=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_allcomp.bed
gmcenter=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_allcomp.center.bed
gmnogene=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_allcomp.center.noTSS.bed
gm2kb=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_allcomp.center.noTSS.2000bp.bed
if [ "$1" == "gm" ];then
  chip=$root/annotations/gm12878/GM12878_CTCF_ChIP.bed
  bedtools intersect -a $bed38 -b $chip -u |\
    sort -k1,1 -k2,2n > $gmbed
  python ../../util/bed_parser.py getcenter -v -b $gmbed |\
    sort -k1,1 -k2,2n > $gmcenter
  bedtools intersect -v -a $gmcenter -b $genes5kb > $gmnogene
  bedtools slop -b 1000 -i $gmnogene -g $gs |\
    sort -k1,1 -k2,2n > $gm2kb
fi
