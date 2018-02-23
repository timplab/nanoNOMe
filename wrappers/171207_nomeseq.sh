#!/bin/bash

root=/dilithium/Data/Nanopore/Analysis/171205_nomeseq
mod=gpc
mod=dam
wdir=$root/$mod
tar=`find /dilithium/Data/Nanopore/oxford/171205_nomeseq -name "*$mod*"`
base=$(basename "$tar")
tag=${base%%.*}
nanopolish=~/Code/nanopolish-$mod/nanopolish
#nanopolish=~/Code/nanopolish/nanopolish
echo $tag

refind=/mithril/Data/NGS/Reference/human38/GRCH38.mmi
ref=/mithril/Data/NGS/Reference/human38/GRCH38.fa
modeldir=$root/model-fofn

if [ 0 -eq 1 ];then
  echo "untar"
  rawdir=$wdir/raw
  [ -e $rawdir ]||mkdir $rawdir
  tar -xvzf $tar -C $rawdir
fi 

if [ 1 -eq 1 ];then
  echo "basecall"
  rawroot=$wdir/raw
  rawdirs=`find $rawroot/* -maxdepth 0 -type d`
  calldir=$wdir/called
  [ -e $calldir ]||mkdir $calldir
  for dir in $rawdirs;do
    ind=$(basename "$dir")
    outdir=$calldir/$ind
    [ -e $outdir ]||mkdir $outdir
    read_fast5_basecaller.py -i $dir -t 12 -s $outdir --output_format fastq --flowcell FLO-MIN106 --kit SQK-LSK108
  done
  fqdir=$wdir/fastq
  [ -e $fqdir ]||mkdir $fqdir
  find $calldir -name "*fastq" -exec cat {} \; > $fqdir/$tag.fq
fi

if [ 0 -eq 1 ];then
  echo "index"
  fqdir=$wdir/fastq
  fq=`find $fqdir -name "*$tag*.fq"`
  logout=$wdir/log/index.log
  $nanopolish index -v -d $wdir/raw $fq 2> $logout
fi

if [ 0 -eq 1 ];then
  echo "align"
  fq=`find $wdir/fastq -name "*$tag*.fq"`
  bamdir=$wdir/bam
  [ -e $bamdir ]||mkdir $bamdir
  bamout=$bamdir/$tag.sorted.bam
  /home/isac/Code/minimap2/minimap2 -ax map-ont $refind $fq | \
  samtools view -q 20 -b - | \
  samtools sort -o $bamout
  samtools index $bamout
fi

if [ 0 -eq 1 ];then
  echo "eventalign"
  fq=`find $wdir/fastq -name "*$tag*.fq"`
  bam=`find $wdir/bam -name "*$tag*.sorted.bam"`
  fofn=`find $modeldir -name "*nuc*"`
  eventdir=$wdir/eventalign
  event=$eventdir/$tag.eventalign
  $nanopolish eventalign --summary $event.summary \
    --scale-events --print-read-names \
    -b $bam \
    -r $fq \
    -g $ref\
    -t 10 > $event
fi

if [ 0 -eq 1 ];then
  echo "call $mod methylation"
  fq=`find $wdir/fastq -name "*$tag*.fq"`
  bam=`find $wdir/bam -name "*$tag*.sorted.bam"`
  fofn=`find $modeldir -name "*$mod*"`
  methdir=$wdir/methcall
  [ -e $methdir ]||mkdir $methdir
  out=$methdir/$tag.$mod.meth.tsv
  $nanopolish call-methylation -t 10 -r $fq -g $ref -b $bam -m $fofn > $out
fi

if [ 0 -eq 1 ];then
  echo "call CpG methylation"
  fq=`find $wdir/fastq -name "*$tag*.fq"`
  bam=`find $wdir/bam -name "*$tag*.sorted.bam"`
  methdir=$wdir/methcall
  out=$methdir/$tag.cpg.meth.tsv
  ~/Code/nanopolish/nanopolish call-methylation -t 10 -r $fq -g $ref -b $bam > $out
fi

if [ 0 -eq 1 ];then
  echo "calculate methylation frequency for $mod"
  methtsv=`find $wdir/methcall -name "*$mod.meth.tsv"`
  freqout=$wdir/methcall/$tag.$mod.meth.freq.tsv
  python ~/Code/ilee/nanopolish/calculate_methylation_frequency.py \
    -i $methtsv -m $mod > $freqout
fi

if [ 0 -eq 1 ];then
  echo "calculate methylation frequency for cpg"
  methtsv=`find $wdir/methcall -name "*cpg.meth.tsv"`
  freqout=$wdir/methcall/$tag.cpg.meth.freq.tsv
  python ~/Code/ilee/nanopolish/calculate_methylation_frequency.py \
    -i $methtsv -m cpg > $freqout
fi

