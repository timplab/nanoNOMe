#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
annodir=$root/annotations/CTCF
[ -e $annodir ]||mkdir $annodir
genes5kb=$root/annotations/hg38/hg38_genes.TSS.5000bp.bed
gs=$root/annotations/hg38/hg38_genomesize.txt

all=$annodir/CTCFBSDB_all.txt.gz
if [ "$1" == "dl" ];then
  wget http://insulatordb.uthsc.edu/download/CTCFBSDB_all_exp_sites_Sept12_2012.txt.gz \
    -O $all
fi
hg19=$annodir/CTCFBSDB_hg19.txt
hg18=$annodir/CTCFBSDB_hg18.txt
if [ "$1" == "separate" ];then
  gunzip -c $all | awk '$3=="hg19"{print}' > $hg19
  gunzip -c $all | awk '$3=="hg18"{print}' > $hg18
fi
bed19=$annodir/CTCFBSDB_hg19.bed
bed18=$annodir/CTCFBSDB_hg18.bed
if [ "$1" == "bed" ];then
  cut -f5,6,7 $hg19 > $bed19
  cut -f5,6,7 $hg18 > $bed18
fi
bed1938=$annodir/CTCFBSDB_hg19_to_hg38.bed
bed1838=$annodir/CTCFBSDB_hg18_to_hg38.bed
if [ "$1" == "liftover" ];then
  chain19=$annodir/../hg38/hg19ToHg38.over.chain
  chain18=$annodir/../hg38/hg18ToHg38.over.chain
  liftOver $bed19 $chain19 $bed1938 $bed1938.unmapped
  liftOver $bed18 $chain18 $bed1838 $bed1838.unmapped
fi
outbed=$annodir/CTCFBSDB_hg38.bed
if [ "$1" == "merge" ];then
  intersect=$annodir/CTCFBSDB_hg38_intersect.bed
  justhg19=$annodir/CTCFBSDB_hg38_from_hg19.bed
  justhg18=$annodir/CTCFBSDB_hg38_from_hg18.bed
  bedtools intersect -v -a $bed1938 -b $bed1838 > $justhg19
  bedtools intersect -v -b $bed1938 -a $bed1838 > $justhg18
  bedtools intersect -wa -a $bed1938 -b $bed1838 | uniq > $intersect
  cat $intersect $justhg19 $justhg38 |\
    sort -k1,1 -k2,2n | uniq > $outbed
fi
sortbed=$annodir/CTCFBSDB_hg19_to_hg38.sorted.bed
if [ "$1" == "sort" ];then
  sort -k1,1 -k2,2n $bed1938 | uniq > $sortbed
fi
cleanbed=$annodir/CTCFBSDB_hg19_to_hg38.clean.bed
if [ "$1" == "cleanup" ];then
  python ./cleanup_bed.py -v -i $sortbed > $cleanbed
fi
notss=$annodir/CTCFBSDB_hg19_to_hg38.clean.noTSS.bed
if [ "$1" == "noTSS" ];then
  bedtools intersect -v -a $cleanbed -b $genes5kb > $notss
fi

gmbed=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_hg19_to_hg38.noTSS.bed
gmcenter=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_hg19_to_hg38.noTSS.center.bed
gm2kb=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_hg19_to_hg38.noTSS.center.2000bp.bed
if [ "$1" == "gm" ];then
  chip=$root/annotations/gm12878/GM12878_CTCF_ChIP.bed
  bedtools intersect -a $notss -b $chip -u > $gmbed
  python ../../util/bed_parser.py getcenter -v -b $gmbed |\
    sort -k1,1 -k2,2n > $gmcenter
  bedtools slop -b 1000 -i $gmcenter -g $gs |\
    sort -k1,1 -k2,2n > $gm2kb
fi
