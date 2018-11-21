#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
dir=$root/subset/ctcf
f5dir=$dir/fast5
[ -e $f5dir ]||mkdir -p $f5dir
bamdir=$root/pooled/bam
fqdir=$root/pooled/fastq
methdir=$root/pooled/methylation/methbyread_all
cell=GM12878
ref="/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa"
s3idx="/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/index/181110_GM12878_s3index.txt"
readdb=
regbed="/dilithium/Data/Nanopore/projects/nomeseq/analysis/annotations/gm12878/GM12878_CTCF_hg38.center.2000bp.bed"
regname="ctcf_2kb"

np=~/Code/nanopolish/nanopolish
methmarker=../../script/convertBam.py

# subset data in this region
prefix=$cell.$regname
bam=$bamdir/$cell.pooled.bam
regbam=$dir/$prefix.bam
fq=$fqdir/$cell.pooled.fastq.gz
regfq=$dir/$prefix.fastq
names=$dir/$prefix.names.txt
fnames=$dir/$prefix.f5names.txt
snp=$dir/$prefix.variants.vcf
phasebam=$dir/$prefix.phased.bam
methbam=$dir/$prefix.methylation.bam
phasemethbam=$dir/$prefix.phasemethylation.bam
if [ "$1" == "subset" ];then
  if [ "$2" == "bam" ];then
    samtools view -hb $bam -L $regbed > $regbam
    samtools index $regbam
  fi
  cut -f1 $regbam > $names
  if [ "$2" == "fastq" ];then
    python ../../util/fqsubset.py -i $fq -n $names -o $regfq
  fi
  if [ "$2" == "fast5" ];then
    grep -f "$names" $readdb | awk 'BEGIN{FS="/"}{ print $NF }'> $fnames
    python ../../util/s3_cp.py -v -b "timp.nanonome" -d "fast5/gm12878" -s $s3idx -i $fnames -o $f5dir
  fi
  if [ "$2" == "meth" ];then
    for mod in cpg gpc;do
      meth=$methdir/$cell.$mod.pooled.meth.bed.gz
      regmeth=$dir/$prefix.$mod.meth.bed.gz
      tabix $meth $reg | bgzip > $regmeth
      tabix -p bed $regmeth
    done
  fi
fi

if [ "$1" == "index" ];then
  $np index -v -d $f5dir $regfq
fi

