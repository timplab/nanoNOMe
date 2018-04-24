#!/bin/bash
root=/shared/Data
tardir=$root/tar
srcroot=/shared/Code
ontpath=$srcroot/ilee/oxford
ontslurm=$ontpath/slurm
logroot=$root/log

nslurm=$srcroot/ilee/dnamods/slurm
nanopolish=$srcroot/nanopolish/nanopolish

target=hg38
if [ "target" == "b37" ];then
  ref=$root/ref/human_g1k_v37_decoy.fasta
else
  ref=$root/ref/GRCH38.fa
fi

for raw in `find $tardir -maxdepth 1 -name "**_2.*tgz" `;do
  base=$(basename "$raw")
  base=${base%%.*}
  echo $base
  wdir=$root/$base
  [ -e $wdir ]||mkdir -p $wdir
  if [ 0 -eq 1 ];then
    logdir=$logroot/untar
    [ -e $logdir ]||mkdir -p $logdir
    echo "untar"
    rawdir=$root/$base/raw
    if [ -e $rawdir ];then 
      echo "already doing/done this, moving on"
      continue
    else
      mkdir $rawdir
    fi
    log=$logdir/$base.untar.log
    sbatch -o $log -e $log -D $rawdir \
      $ontslurm/untar.scr $raw 
  fi
  if [ 0 -eq 1 ];then
    echo "basecall"
    logdir=$logroot/basecall
    [ -e $logdir ]||mkdir $logdir
    rawdir=$wdir/raw
    calldir=$wdir/called
    [ -e $calldir ]||mkdir $calldir
    # identify fast5 folder path and perform basecalling
    $ontslurm/call_wrapper.sh -i $rawdir -o $calldir \
      -f FLO-MIN106 -s $ontslurm --aws &> $logdir/$base.callwrapper.log
  fi
  if [ 0 -eq 1 ];then
    echo "cat fastq"
    calldir=$wdir/called
    fqdir=$wdir/fastq
    [ -e $fqdir ]||mkdir $fqdir
    logdir=$logroot/fqcat
    for dir in `find $calldir/* -maxdepth 0 -type d`;do
      ind=$(basename "$dir")
      fqpre=$fqdir/$base.$ind
      log=$logdir/$base.$ind.cat.log
      echo $fqpre
      sbatch -e $log -o $log \
        $ontslurm/cat_bcall.scr -i $dir -p $fqpre
    done
  fi
  if [ 0 -eq 1 ];then
    echo "index"
    rawdir=$wdir/raw
    logdir=$logroot/index
    fqdir=$wdir/fastq
    [ -e $logdir ]||mkdir $logdir
    for dir in `find $rawdir/* -maxdepth 0 -type d`;do
      ind=$(basename "$dir")
      echo $ind
      log=$logdir/$base.$ind.index.log
      fq=`find $fqdir -name "*$base*.$ind.*fq.gz"`
#      echo $fq
      sbatch -e $log -o $log $nslurm/fast5index.scr \
        -s $nanopolish -i $fq -d $dir
    done
  fi
  if [ 0 -eq 1 ];then
    echo "align"
    fqlist=`find $wdir/fastq -name "*fq.gz"`
    logdir=$logroot/align-$target
    [ -e $logdir ]||mkdir $logdir
    bamdir=$wdir/bam-$target
    [ -e $bamdir ]||mkdir $bamdir
    for fq in $fqlist;do
      b=$(basename "$fq")
      t=${b%.fq.gz}
      echo $t
      log=$logdir/$t.align.log
      pre=$bamdir/$t
      sbatch -e $log -o $log -c 36 $ontslurm/oxford_align.scr \
        -r $ref -i $fq -a minimap2 -b $pre --aws --samarg "-q 20"
    done
  fi
  if [ 1 -eq 1 ];then
    mod=${base%_nome*}
    mod=${mod#*_}
#    mod=cpg
    echo $mod
    methdir=$wdir/methcall/${mod}-$target
    [ -e $methdir ]||mkdir -p $methdir
    nroot=$srcroot/nanopolish
    echo "call $mod methylation"
    fqlist=`find $wdir/fastq -name "**fq.gz"`
    logdir=$logroot/methcall-$target
    [ -e $logdir ]||mkdir $logdir
    bamdir=$wdir/bam-$target
    for fq in $fqlist;do
      b=$(basename "$fq")
      t=${b%.fq.gz}
      echo $t
      bam=`find $bamdir -name "$t.sorted.bam"`
      log=$logdir/$t.$mod.methcall.log
      methout=$methdir/$t.$mod.meth.tsv
      sbatch -o $log -e $log -c 36 $nslurm/callmeth.scr \
        -s $nroot -r $fq -b $bam -g $ref -m $mod -o $methout --aws
    done
  fi
done
