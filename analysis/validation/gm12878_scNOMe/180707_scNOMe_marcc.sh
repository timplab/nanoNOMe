#!/bin/bash
# from : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5487215/
root=/scratch/users/ilee29@jhu.edu/NGS/gm12878/nomeseq
logroot=$root/log
tmproot=$root/tmp
covdir=$root/covreport
freqdir=$root/methfreq
refdir=/scratch/groups/wtimp1/Reference/human/hg38_noalt
srcdir=/home-2/ilee29@jhu.edu/Code/ilee/methylation/bismark/slurm

fqdir=$root/fastq
trimdir=$root/fqtrim
[ -e $trimdir ]||mkdir $trimdir
if [ "$1" == "fqtrim" ];then
  tmpdir=$tmproot/fqtrim
  for fq1 in $(find $fqdir -name "*R1.fastq.gz");do
    samp=$(basename "$fq1")
    samp=${samp%%_R1*}
    tmp=$tmpdir/$samp
    [ -e $tmp ]||mkdir -p $tmp
    log=$trimdir/$samp.trim.log
    com="sbatch -e $log -o $log -J fqtrim-$samp \
      -p shared -c 1 -t 48:0:0 --mem=20GB \
      $srcdir/fastqTrim.scr -p $fqdir/$samp -o $trimdir \
      -t $tmp --marcc"
    echo $com
    eval $com
  done
fi

bamdir=$root/bam
[ -e $bamdir ]||mkdir $bamdir
if [ "$1" == "align" ];then
  logdir=$logroot/align
  [ -e $logdir ]||mkdir -p $logdir
  for fq1 in $(find $trimdir -name "*sample*rep1_*val_1.fq.gz");do
    samp=$(basename "$fq1")
    samp=${samp%%_R1*}
    echo $samp
    [ -e $tmpdir ]||mkdir -p $tmpdir
    log=$logdir/$samp.align.log
    com="sbatch -e $log -o $log -J align-$samp \
      -p shared -c 24 -t 24:0:0 --mem=120GB \
      $srcdir/bismarkAlign.scr -i $trimdir/$samp -o $bamdir/$samp \
      -r $refdir -t $tmproot --marcc"
    echo $com
    eval "$com"
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
if [ "$1" == "merge" ];then
  for samp in sample control;do
    bams=$(find $bamdir -name "*$samp*nodupqsort.bam")
    outbam=$bamdir/GM12878_$samp.pooled.bam
    echo $outbam
    samtools merge -fn $outbam $bams 
  done
fi
extractdir=$root/extract
[ -e $extractdir ]||mkdir $extractdir
if [ "$1" == "extract" ];then
  logdir=$logroot/extract
  [ -e $logdir ]||mkdir $logdir
  for bam in $(find $bamdir -name "*pooled.bam");do
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
reportdir=$root/report
[ -e $reportdir ]||mkdir $reportdir
if [ "$1" == "report" ];then
  logdir=$logroot/report
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/{}.report.log
  for bam in $(find $bamdir -name "*pooled.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    samps="$samps $samp"
  done
  com="bismark2bedGraph --counts --CX --dir $reportdir \
    -o {}.bedGraph $extractdir/*{}.*txt.gz &> $log"
  parallel $com ::: $samps
fi

freqdir=$root/methfreq
if [ "$1" == "getfreq" ];then
  for bam in $(find $bamdir -name "*pooled.bam");do
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
