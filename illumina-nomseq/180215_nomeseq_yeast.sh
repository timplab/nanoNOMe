#!/bin/bash
rawdir=/dilithium/Data/NGS/Raw/180215_nomeseq
outroot=/dilithium/Data/NGS/Aligned/180215_nomeseq/yeast
refdir=/mithril/Data/NGS/Reference/yeast
bismarkdir=/home/isac/Code/Bismark_v0.19.0
picard=~/Code/picard.jar

trimdir=$outroot/fqtrim
bamdir=$outroot/bam
logroot=$outroot/log
tmpdir=$outroot/tmp
reportdir=$outroot/report
beddir=$outroot/bed

for r1 in `find $rawdir -name "*yeast*R1*"`;do
  fqpre=${r1%R1*}
  r2=`readlink -f $fqpre*R2*fastq.gz`
  base=$(basename "$r1")
  base=${base%%_*}
  echo $base
  if [ "$1" == "trim" ];then
    logdir=$logroot/fqtrim
    [ -e $logdir ]||mkdir -p $logdir
    [ -e $trimdir ]||mkdir $trimdir
    trim_galore -q 28 --paired \
      $r1 $r2 -o $trimdir &> $logdir/$base.fqtrim.log
  fi
  trim1=`find $trimdir -name "${base}_*val_1.fq.gz"`
  trim2=`find $trimdir -name "${base}_*val_2.fq.gz"`
  if [ "$1" == "align" ];then
    logdir=$logroot/align
    [ -e $logdir ]||mkdir -p $logdir
    [ -e $bamdir ]||mkdir -p $bamdir
    $bismarkdir/bismark --bam --non_direction --bowtie2 \
      --un -p 4 --genome $refdir \
      -1 $trim1 -2 $trim2 -B $base \
      --temp_dir $tmpdir --output_dir $bamdir &> $logdir/$base.align.log
  fi
  bampre=$bamdir/${base}
  if [ "$1" == "removedup" ];then
    logdir=$logroot/removedup
    [ -e $logdir ]||mkdir $logdir
    samtools sort -T $tmpdir/$base.sorted -o $bampre.sorted.bam ${bampre}_pe.bam
    echo "done sorting"
    java -Xmx10G -jar $picard MarkDuplicates \
      I=$bampre.sorted.bam \
      O=$bampre.nodup.bam \
      M=$bampre.dupmetrics.txt \
      REMOVE_DUPLICATES=true &> $logdir/$base.removedup.log
    echo "done removing duplicates"
  fi
  if [ "$1" == "extract" ];then
    logdir=$logroot/extract
    [ -e $logdir ]||mkdir $logdir
    [ -e $reportdir ]||mkdir $reportdir
    samtools sort -n -T $tmpdir/$base.qsort -o $bampre.qsort.bam $bampre.nodup.bam
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
    logdir=$logroot/report
    [ -e $logdir ]||mkdir $logdir
    $bismarkdir/coverage2cytosine --NOMe-seq --gzip \
      --dir $reportdir --genome_folder $refdir \
      -o $base $reportdir/$base.bismark.cov.gz &> $logdir/$base.nome.log
  fi
  if [ "$1" == "cyto2bed" ];then
    [ -e $beddir ]||mkdir $beddir
    ~/Code/ilee/bsseq/scripts/cytoreport2bedgraph.sh \
      $reportdir/$base.NOMe.GpC_report.txt.gz | \
      gzip > $beddir/$base.cytoreport.gpc.bedGraph.gz
  fi

done 
