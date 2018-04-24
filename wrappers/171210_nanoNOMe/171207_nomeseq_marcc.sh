#!/bin/bash
tardir=/scratch/users/ilee29@jhu.edu/tar
tag=dam
root=/scratch/users/ilee29@jhu.edu/171207_nomeseq
srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford
slurmpath=$srcpath/slurm
logroot=$root/log

nslurm=/home-2/ilee29@jhu.edu/Code/ilee/dnamods/slurm
nanopolish=/home-2/ilee29@jhu.edu/Code/nanopolish/nanopolish

ref=/scratch/groups/wtimp1/Reference/human/hg38/GRCH38.fa
ref=/scratch/groups/wtimp1/Reference/human/b37/human_g1k_v37_decoy.fasta

for raw in `find $tardir -name "*_$tag*" `;do
  base=$(basename "$raw")
  base=${base%%.*}
  echo $base
  wdir=$root/$base
  [ -e $wdir ]||mkdir -p $wdir
  if [ 0 -eq 1 ];then
    logdir=$logroot/untar
    [ -e $logdir ]||mkdir $logdir
    echo "untar"
    rawdir=$root/$base/raw
    [ -e $rawdir ]||mkdir $rawdir
    log=$logdir/$base.untar.log
    sbatch -o $log -e $log -D $rawdir \
      $slurmpath/untar.scr $raw 
  fi
  if [ 0 -eq 1 ];then
    echo "basecall"
    logdir=$logroot/basecall
    [ -e $logdir ]||mkdir $logdir
    rawdir=$root/$base/raw
    calldir=$wdir/called
    [ -e $calldir ]||mkdir $calldir
    # identify fast5 folder path and perform basecalling
    $slurmpath/call_wrapper.sh -i $rawdir -o $calldir \
      -f FLO-MIN106 -s $slurmpath --marcc &> $logdir/callwrapper.log
  fi
  if [ 0 -eq 1 ];then
    echo "cat fastq"
    calldir=$wdir/called
    fqdir=$wdir/fastq
    [ -e $fqdir ]||mkdir $fqdir
    for dir in `find $calldir/* -maxdepth 0 -type d`;do
      ind=$(basename "$dir")
      fqpre=$fqdir/$base.$ind
      log=$logroot/$base.$ind.cat.log
      echo $fqpre
      sbatch -e $log -o $log --partition=shared \
        $slurmpath/cat_bcall.scr -i $dir -p $fqpre
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
      echo $fq
      sbatch --partition=shared -e $log -o $log $nslurm/fast5index.scr \
        -s $nanopolish -i $fq -d $dir
    done
  fi
  if [ 0 -eq 1 ];then
    echo "align"
    fqlist=`find $wdir/fastq -name "*fq.gz"`
    logdir=$logroot/align-b37
    [ -e $logdir ]||mkdir $logdir
    bamdir=$wdir/bam-b37
    [ -e $bamdir ]||mkdir $bamdir
    for fq in $fqlist;do
      b=$(basename "$fq")
      t=${b%.fq.gz}
      echo $t
      log=$logdir/$t.align.log
      pre=$bamdir/$t
      sbatch --partition=shared -e $log -o $log $slurmpath/oxford_align.scr \
        -r $ref -i $fq -a minimap2 -b $pre --marcc --samarg "-q 20"
    done
  fi
  if [ 1 -eq 1 ];then
    mod=dam
    methdir=$wdir/methcall/$mod-b37
    [ -e $methdir ]||mkdir -p $methdir
    nanopolish=/home-2/ilee29@jhu.edu/Code/nanopolish-dam/nanopolish
    echo "call methylation"
    fqlist=`find $wdir/fastq -name "*fq.gz"`
    logdir=$logroot/methcall
    [ -e $logdir ]||mkdir $logdir
    bamdir=$wdir/bam-b37
    for fq in $fqlist;do
      b=$(basename "$fq")
      t=${b%.fq.gz}
      echo $t
      bam=`find $bamdir -name "$t.sorted.bam"`
      log=$logdir/$t.$mod.methcall.log
      methout=$methdir/$t.$mod.meth.tsv
      sbatch -o $log -e $log $nslurm/callmeth.scr \
        -s $nanopolish -r $fq -b $bam -g $ref -o $methout
    done
  fi
done
