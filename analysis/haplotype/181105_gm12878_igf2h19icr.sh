#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
dir=$root/subset
[ -e $dir ]||mkdir $dir
f5dir=$dir/fast5/gm12878
[ -e $f5dir ]||mkdir -p $f5dir
bamdir=$root/pooled/bam
fqdir=$root/pooled/fastq
methdir=$root/pooled/methylation/methbyread_all
cell=GM12878
ref="/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa"
s3idx="/dilithium/Data/Nanopore/projects/nomeseq/raw/index/GM12878_s3index.txt"
reg="chr11:1998745-2003509"
regname="igf2h19icr"

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
    samtools view -hb $bam $reg > $regbam
    samtools index $regbam
  fi
  cut -f1 $regbam > $names
  if [ "$2" == "fastq" ];then
    python ../../util/fqsubset.py -i $fq -n $names -o $regfq
  fi
  if [ "$2" == "fast5" ];then
    readdb=$fqdir/$cell.pooled.fastq.gz.index.readdb
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

if [ "$1" == "snp" ];then
  $np variants -r $regfq -b $regbam -g $ref \
    -p 2 -w $reg > $snp
fi

if [ "$1" == "phase" ];then
  log=$dir/$prefix.phase.log
  $np phase-reads --progress -v -t 10 -r $regfq -b $regbam -g $ref \
    -w $reg $snp 2> $log | samtools sort > $phasebam 
  samtools index $phasebam
fi

if [ "$1" == "markmeth" ];then
  cpg=$dir/$prefix.cpg.meth.bed.gz
  gpc=$dir/$prefix.gpc.meth.bed.gz
  python $methmarker -v -b $regbam -c $cpg -g $gpc -w "$reg" |\
    samtools sort > $methbam
  python $methmarker -v -b $phasebam -c $cpg -g $gpc -w $reg |\
    samtools sort > $phasemethbam
  samtools index $methbam
  samtools index $phasemethbam
fi
