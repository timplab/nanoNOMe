#!/bin/bash
echo "parsing enhancers"
root="$1"

# hg38 annotation
dir=$root/hg38
prefix=$dir/hg19_enhancers
txt=$prefix.txt.gz
# convert to bed
bed=$prefix.bed
[ -e $bed ]||\
  gunzip -c $txt | awk 'OFS="\t"{ print $2,$3,$4,$5,".","." }' > $bed
# liftover
lift=$root/liftOver
if [ ! -e $lift ];then
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver \
  -O $lift
  chmod a+x $lift
fi
chain=$root/hg19ToHg38.over.chain
if [ ! -e $chain ];then
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
  -O $chain.gz
  gunzip $chain.gz
fi
new=$dir/hg38_enhancers.bed
if [ ! -e $new ];then
  $lift $bed $chain $new.unsorted $new.unmapped
  sort -k1,1 -k2,2n $new.unsorted > $new
  rm $new.unsorted
fi
# combine with promoters
combine=$dir/hg38_enhancers_promoters.bed
[ -e $combine ]||\
  cat $new $dir/hg38_genes.TSS.400bp.bed |\
  sort -k1,1 -k2,2n > $combine

# for gm12878
pre=$root/gm12878/enhancer/GM12878_enhancer_promoter
txt=${pre}_hg19.txt
bed=${pre}_hg19.bed
[ -e $bed ]||\
  awk 'OFS="\t"{ print $1,$2,$3,$7,".","." }' $txt > $bed
# liftover
new=${pre}_hg38.bed
if [ ! -e $new ];then
  $lift $bed $chain $new.unsorted $new.unmapped
  sort -k1,1 -k2,2n $new.unsorted > $new
  rm $new.unsorted
fi

