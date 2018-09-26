#!/bin/bash
rawdir=/dilithium/Data/NGS/Raw/180920_NA12878_methylation
outdir=/dilithium/Data/NGS/Aligned/180920_NA12878_methylation
logroot=$outdir/log
[ -e $logroot ]||mkdir $logroot
tmproot=$outdir/tmp
[ -e $tmproot ]||mkdir $tmproot
covdir=$outdir/covreport
[ -e $covdir ]||mkdir $covdir
freqdir=$outdir/methfreq
[ -e $freqdir ]||mkdir $freqdir
trimdir=$outdir/fqtrim
[ -e $trimdir ]||mkdir $trimdir
bamdir=$outdir/bam
[ -e $bamdir ]||mkdir $bamdir
extractdir=$outdir/extract
[ -e $extractdir ]||mkdir $extractdir
reportdir=$outdir/report
[ -e $reportdir ]||mkdir $reportdir
freqdir=$outdir/methfreq
[ -e $fredir ]||mkdir $freqdir
refdir=/mithril/Data/NGS/Reference/hg38_noalt

if [ "$1" == "fqtrim" ];then
  for fq in $(find $rawdir -name "*R1*.fastq.gz");do
    samp=$(basename "$fq")
    samp=${samp%%_R1*}
    echo $samp
    samples="$samples $samp"
  done
  fq1=$rawdir/{}_R1_001.fastq.gz
  fq2=$rawdir/{}_R2_001.fastq.gz
  log=$trimdir/{}.trim.log
  trim_galore -v > $logroot/trim.version
  com="trim_galore -q 28 --paired \
    $fq1 $fq2 -o $trimdir &> $log"
  parallel $com ::: $samples
fi

if [ "$1" == "align" ];then
  logdir=$logroot/align
  [ -e $logdir ]||mkdir $logdir
  for fq1 in $(find $trimdir -name "*val_1.fq.gz");do
    samp=$(basename "$fq1")
    samp=${samp%%_R1*}
    echo $samp
    fq2=$(find $trimdir -name "${samp}_*2.fq.gz")
    tmpdir=${tmproot}/align/$samp
    [ -e $tmpdir ]||mkdir -p $tmpdir
    log=$logdir/$samp.align.log
    bismark --version > $logdir/bismark.version
    bismark --bam --non_directional --bowtie2 --un -p 2 \
      --genome $refdir \
      -1 $fq1 \
      -2 $fq2 \
      --temp_dir $tmpdir \
      --output_dir $bamdir \
      -B $samp &> $log
  done
fi
if [ "$1" == "removedup" ];then
  tmpdir=$tmproot/removedup
  [ -e $tmpdir ]||mkdir $tmpdir
  logdir=$logroot/removedup
  [ -e $logdir ]||mkdir $logdir
  for bam in $(find $bamdir -name "*pe.bam");do
    samp=$(basename "$bam")
    samp=${samp%%_pe*}
    echo $samp
    samtools sort -@8 -m 2G -T $tmpdir/$samp.sorted -o $tmpdir/$samp.sorted.bam $bam
    echo "done sorting"
    log=$logdir/$samp.removedup.log
    com="picard -Djava.io.tmpdir=$tmpdir -Xmx5G MarkDuplicates \
      I=$tmpdir/$samp.sorted.bam \
      O=$tmpdir/$samp.nodup.bam \
      M=$bamdir/$samp.dupmetrics.txt \
      REMOVE_DUPLICATES=true &> $log"
    echo $com
    eval $com
    samtools sort -n -@ 8 -m2G -T $tmpdir/$samp.qsort \
      -o $bamdir/$samp.nodupqsort.bam $tmpdir/$samp.nodup.bam
  done
fi
if [ "$1" == "extract" ];then
  logdir=$logroot/extract
  [ -e $logdir ]||mkdir $logdir
  for bam in $(find $bamdir -name "*nodupqsort.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    log=$logdir/$samp.extract.log
    args="--ignore 2 --ignore_r2 2 --ignore_3prime 1 --ignore_3prime_r2 1"
    com="bismark_methylation_extractor -p --multicore 8 --gzip \
      --genome_folder $refdir --no_header $args \
      $bam -o $extractdir &> $log"
    echo $com
    eval $com
  done
fi
if [ "$1" == "report" ];then
  logdir=$logroot/report
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/{}.report.log
  for bam in $(find $bamdir -name "*nodupqsort.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    samps="$samps $samp"
  done
  com="bismark2bedGraph --counts --CX --dir $reportdir \
    -o {}.bedGraph $extractdir/*{}.*txt.gz &> $log"
  parallel $com ::: $samps
fi

if [ "$1" == "getfreq" ];then
  for bam in $(find $bamdir -name "*nodupqsort.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    samps="$samps $samp"
  done
  cov=$reportdir/{}.bismark.cov.gz
  logdir=$logroot/methfreq
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/{}.methfreq.log
  com="coverage2cytosine --NOMe-seq --gzip --dir $freqdir \
    --genome_folder $refdir \
    -o {} $cov &> $log"
  echo $com
  parallel "$com" ::: $samps
fi

if [ "$1" == "summarize" ];then
  summary=$freqdir/180922_methylation_summary.txt
  [ ! -e $summary ]||rm $summary
  for bam in $(find $bamdir -name "*nodupqsort.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    samps="$samps $samp"
  done
  for samp in $samps;do
    for mod in CpG GpC;do
      report=$freqdir/$samp.NOMe.${mod}_report.txt.gz
      sum=$( gunzip -c $report |\
        awk '{ meth+=$4;unmeth+=$5 }END{ print meth,unmeth,meth/(meth+unmeth) }')
      echo "$samp $mod $sum" | tr " " "\t" >> $summary
    done
  done
fi
  
if [ "$1" == "html" ];then
  find $bamdir -name "*report*" -exec cp \{\} $extractdir \;
  cd $extractdir
  bismark2report
fi

if [ "$1" == "mbias" ];then
  for bam in $(find $bamdir -name "*nodupqsort.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    samtools sort $bam
    python ../../../script/bsseq_mbias.py -v -t 8 -b $bam
  done
fi
# not yet updated
if [ "$1" == "tabix" ];then
  tmpdir=$tmproot/tabix
  for samp in sample control;do
    for mod in CpG GpC;do
      cyto=$(find $freqdir -name "*$samp.*$mod*report.txt.gz")
      pre=${cyto%%.*}
      lab=$(echo $mod | tr "[:upper:]" "[:lower:]")
      mkdir -p $tmpdir/$samp/$lab
      out=$pre.$lab.methfreq.txt.gz
      echo $out
      com="gunzip -c $cyto |\
        sort -T $tmpdir/$samp/$lab -k1,1 -k2,2n |\
        bgzip > $out && 
      tabix -b 2 -e 2 $out"
      eval "$com"
    done
  done
fi


