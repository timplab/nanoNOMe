#!/bin/bash
t=92
root=/data
samp=181102_GM12878_NOMe
outdir=$root/analysis/$samp
[ -e $outdir ]||mkdir -p $outdir
logroot=$outdir/log
[ -e $logroot ]||mkdir $logroot
bamdir=$outdir/bam
[ -e $bamdir ]||mkdir $bamdir
rawroot=$root/sorted/$samp
ref=/home/prom/Data/Reference/hg38_noalt/hg38_noalt.fa
fqs=$(find $rawroot -maxdepth 2 -type f -name "*fastq")
for fq in $fqs;do
  fqdir=$(dirname "$fq")
  idx=${fq##*_}
  idx=${idx%%.*}
  inds="$inds $idx"
done
#inds=$(echo $inds | tr " " "\n" | head -n30)
echo $inds

if [ "$1" == "index" ];then
  logdir=$logroot/index
  [ -e $logdir ]||mkdir $logdir
  log="$logdir/$samp.{}.index.log"
  fq="$fqdir/fastq_{}.fastq"
  sum="$fqdir/sequencing_summary_{}.txt"
  f5dir="$rawroot/fast5/{}"
  com="nanopolish index -v -d $f5dir -s $sum $fq &> $log"
  parallel "$com" ::: $inds
fi

if [ "$1" == "align" ];then
  logdir=$logroot/align
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/$samp.{}.align.log
  fq=$fqdir/fastq_{}.fastq
  bam=$bamdir/$samp.{}.bam
  # parallelize with 30 jobs to minimize RAM overusage
  com="ngmlr -t 3 -r $ref -q $fq -x ont 2> $log |\
    samtools view -q 20 -b - |\
    samtools sort -@ 3 -T $bam.sorting -o $bam &&\
    samtools index $bam"
  parallel -j 30 "$com" ::: $inds
fi

if [ "$1" == "mcall" ];then
  for mod in cpg gpc;do
    mdir=$outdir/$mod
    [ -e $mdir ]||mkdir $mdir
    logdir=$logroot/$mod
    [ -e $logdir ]||mkdir $logdir
    log=$logdir/$samp.{}.$mod.log
    fq=$fqdir/fastq_{}.fastq
    bam=$bamdir/$samp.{}.bam
    out=$mdir/$samp.{}.$mod.tsv
    com="nanopolish call-methylation -t 3 --progress\
      -r $fq -b $bam -g $ref -q $mod > $out 2> $log"
    echo $com
    parallel -j 30 "$com" ::: $inds
  done
fi

if [ "$1" == "mbed" ];then
  converter="../../nanopolish/mtsv2bedGraph.py"
  for mod in gpc;do
    bdir=$outdir/mbed_$mod
    [ -e $bdir ]||mkdir $bdir
    logdir=$logroot/mbed_$mod
    [ -e $logdir ]||mkdir $logdir
    log=$logdir/$samp.{}.mbed.$mod.log
    tsv=$outdir/$mod/$samp.{}.$mod.tsv
    bed=$bdir/$samp.{}.meth.bed
    com="python $converter -i $tsv -m $mod > $bed 2> $log"
    parallel "$com" ::: $inds
  done
fi

if [ "$1" == "merge" ];then
  for mod in cpg gpc;do
#    dir=$outdir/$mod
#    out=$outdir/$samp.$mod.meth.tsv.gz
#    find $dir -name "*$mod.*tsv" -exec cat {} \; |\
#      awk 'NR == 1 || $1 != "chromosome" { print }' |\
#      gzip > $out
    # mbed
    dir=$outdir/mbed_$mod
    out=$outdir/$samp.$mod.meth.bed.gz
    find $dir -name "$samp*.*bed" -exec cat {} \; |\
      sort -k1,1 -k2,2n | bgzip > $out
    tabix -p bed $out
  done
#  echo bam
#  outbam=$outdir/$samp.bam
#  samtools merge -f $outbam $outdir/bam/*bam
#  echo fq
#  outfq=$outdir/$samp.fastq.gz
#  find $fqdir -type f -name "*fastq" -exec cat {} \; |\
#    gzip > $outfq
fi
