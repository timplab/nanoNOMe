#!/bin/bash
rawdir=/dilithium/Data/NGS/Raw/180215_nomeseq
outroot=/dilithium/Data/NGS/Aligned/180215_nomeseq/subset
refdir=/mithril/Data/NGS/Reference/hg38_noalt
bismarkdir=/home/isac/Code/Bismark_v0.19.0
picard=~/Code/picard.jar

meta=`readlink -f  $outroot/*human.txt`
trimdir=$outroot/fqtrim
bamdir=$outroot/bam
logroot=$outroot/log
tmpdir=$outroot/tmp
reportdir=$outroot/report
beddir=$outroot/bed

for r1 in `find $rawdir -name "*human*50*R1*"`;do
  fqpre=${r1%R1*}
  r2=`readlink -f $fqpre*R2*fastq.gz`
  base=$(basename "$r1")
  base=${base%%_*}
  echo $base
  bampre=$bamdir/$base
  if [ "$1" == "subsample" ];then
    [ -e $bamdir ]||mkdir $bamdir
    samtools view -hb $outroot/../bam/$base.nodup.bam -s 1.10 > $bampre.subsample.bam
  fi

  if [ "$1" == "extract" ];then
    logdir=$logroot/extract
    [ -e $logdir ]||mkdir $logdir
    [ -e $reportdir ]||mkdir $reportdir
    samtools sort -n -T $tmpdir/$base.qsort -o $bampre.qsort.bam $bampre.subsample.bam
    $bismarkdir/bismark_methylation_extractor -p --multicore 8 --gzip \
      --genome_folder $refdir --multicore 3 \
      --ignore_r2 1 --ignore_3prime_r2 1 \
      $bampre.qsort.bam -o $reportdir --no_header &> $logdir/$base.extract.log
  fi
  if [ "$1" == "report" ];then
    logdir=$logroot/report
    [ -e $logdir ]||mkdir $logdir
    [ -e $reportdir ]||mkdir $reportdir
    $bismarkdir/bismark2bedGraph --counts --CX \
      --dir $reportdir -o $base.bedGraph ${reportdir}/*$base.*txt.gz &> $logdir/$base.bedgraph.log
  fi
  if [ "$1" == "splitreport" ];then
    ~/Code/ilee/util/splitByChromosome.py ${reportdir}/$base.bismark.cov.gz
  fi
  if [ "$1" == "nome" ];then
    logdir=$logroot/nome
    [ -e $logdir ]||mkdir $logdir
    find $reportdir/$base -name *.bismark.cov.gz |\
      parallel --jobs 4 \
      $bismarkdir/coverage2cytosine --NOMe-seq --gzip \
        --dir $reportdir/$base --genome_folder $refdir \
        -o {#} {} 

#    for chr in do
#      b=$(basename "$chr")
#      b=${b%.bismark.cov.gz}
#      echo $b
#      log=$logdir/$b.nome.log
#      $bismarkdir/coverage2cytosine --NOMe-seq --gzip \
#        --dir $reportdir/$base --genome_folder $refdir \
#        -o $b $chr &> $log
#    done
  fi
  if [ "$1" == "cyto2bed" ];then
    [ -e $beddir ]||mkdir $beddir
    ~/Code/ilee/bsseq/scripts/cytoreport2bedgraph.sh \
      $reportdir/$base.NOMe.GpC_report.txt.gz | \
      gzip > $beddir/$base.cytoreport.gpc.bedGraph.gz
  fi

done 
