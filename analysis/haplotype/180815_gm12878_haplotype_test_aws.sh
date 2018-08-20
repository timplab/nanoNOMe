#!/bin/bash
root=/shared/data/analysis/gm12878
pooldir=$root/pooled
rawroot=/shared/data/raw/gm12878/fast5
cell=GM12878

npdir=/shared/Code/nanopolish
np=$npdir/nanopolish
wrapperdir=/shared/Code/nanoNOMe/nanopolish/slurm
ref=/shared/data/Reference/hg38_noalt.fa
batchcom=/shared/Code/nanoNOMe/script/batchcommand.scr

fq=$pooldir/$cell.pooled.fastq.gz
bam=$pooldir/$cell.pooled.bam
sum=$pooldir/$cell.pooled.summary.txt
# prep data
if [ "$1" == "prep" ];then
  if [ ! -e $sum ];then
    sums=$(find $root/fastq -name "*summary.txt")
    awk 'FNR==1 && NR!=1 {next}{print}' $sums > $sum
  fi
  index=$pooldir/$cell.pooled.fastq.gz.index
  if [ ! -e $index ];then
    inds=$(find $root/fastq -name "*index")
    cat $inds > $index
  fi
  fai=$pooldir/$cell.pooled.fastq.gz.index.fai
  if [ ! -e $fai ];then
    fais=$(find $root/fastq -name "*fai")
    cat $fais > $fai
  fi
  gzi=$pooldir/$cell.pooled.fastq.gz.index.gzi
  if [ ! -e $gzi ];then
    gzis=$(find $root/fastq -name "*gzi")
    cat $gzis > $gzi
  fi
  readdb=$pooldir/$cell.pooled.fastq.gz.index.readdb
  if [ ! -e $readdb ];then
    rdbs=$(find $root/fastq -name "*readdb")
    cat $rdbs > $readdb
  fi
fi
if [ "$1" == "readdb" ];then
  com="$np index -v -d $rawroot -s $sum $fq"
  echo $com
  log=$pooldir/$cell.index.log
  sbatch -e $log -o $log -c 1 -J "index" $batchcom $com
fi

if [ "$1" == "variants" ];then
  snpdir=$root/snps
  [ -e $snpdir ]||mkdir $snpdir
  segs=$snpdir/segments.txt
  [ -e $segs ]||\
    python $npdir/scripts/nanopolish_makerange.py --segment-length 300000 $ref > $segs
  numsegs=$(wc -l "$segs" | awk '{ print $1 }')
  array="1-$numsegs"
  array="4-5"
  echo $array
  log=$snpdir/$cell.%a.variants.log
  com="head -n %a $segs | tail -n1 |\
    xargs -Ireg $np variants -v -o $snpdir/$cell.%a.vcf -w reg -p 2 -r $fq -b $bam -g $ref -t 36 --min-candidate-frequency 0.1" 
  sbatch -e $log -o $log -c 36 -J %a-variants -a $array $batchcom \
    $com
fi

if [ "$1" == "variants-guided" ];then
  snpdir=$root/snps
  segs=$snpdir/segments.txt
  array="4-5"
  echo $array
  log=$snpdir/$cell.%a.variants.guided.log
  vcf=/shared/data/annotation/NA12878.vcf
  com="head -n %a $segs | tail -n1 |\
    xargs -Ireg $np variants -v -o $snpdir/$cell.%a.guided.vcf -w reg -p 2 -r $fq -b $bam -g $ref -t 36 --min-candidate-frequency 0.1 -c $vcf" 
  sbatch -e $log -o $log -c 36 -J %a-variants-guided -a $array $batchcom \
    $com
fi
  
# number of possible settings, np-np,guidednp-np,vcf-np,vcf-hapcut2,np-hapcut2
# I'm going to try : np-np,vcf-hapcut2 to check - benchmarking another time
if [ "$1" == "haplotype" ];then
  snpdir=$root/snps
  segs=$snpdir/segments.txt
  array="4-5"
  echo $array
  # np-np
  var=$snpdir/$cell.%a.vcf
  log=$snpdir/$cell.%a.phase.log
  com="head -n %a $segs | tail -n1 |\
    xargs -Ireg $np phase-reads -v --progress -r $fq -b $bam -g $ref -t 36 -w reg $var > $snpdir/$cell.%a.phase.txt"
  sbatch -e $log -o $log -c 36 -J phase -a $array $batchcom $com
fi


