#!/bin/bash
rawdir=/dilithium/Data/NGS/Raw/181003_chosigmaATACrep1_methylationStandards
outdir=/dilithium/Data/NGS/Aligned/181003_methylationStandards
trimdir=$outdir/fqtrim
[ -e $trimdir ]||mkdir -p $trimdir
tmproot=$outdir/tmp
[ -e $tmproot ]||mkdir -p $tmproot
logroot=$outdir/log
nref=/mithril/Data/NGS/Reference/hg38_noalt
eref=/mithril/Data/NGS/Reference/ecoli/bisulfite.mg1655

# get sample names
efqs=$(find $rawdir -name "Ecoli*R1*fastq.gz")
nfqs=$(find $rawdir -name "NA12878*R1*fastq.gz")
for fq in $efqs;do
  base=$(basename "$fq")
  samp=${base%%_*}
  ecoli="$ecoli $samp"
done
for fq in $nfqs;do
  base=$(basename "$fq")
  samp=${base%%_*}
  na="$na $samp"
done

if [ "$1" == "trim" ];then
  echo "trim"
  logdir=$logroot/fqtrim
  [ -e $logdir ]||mkdir -p $logdir
  log=$logdir/{}.trim.log
  fq1=$rawdir/{}_*R1*fastq.gz
  fq2=$rawdir/{}_*R2*fastq.gz
  com="trim_galore -v > $log &&\
    trim_galore --paired --fastqc \
    $fq1 $fq2 -o $trimdir &>> $log"
  parallel "$com" ::: $na $ecoli
fi

bamdir=$outdir/bam
[ -e $bamdir ]||mkdir $bamdir
if [ "$1" == "align" ];then
  logdir=$logroot/align
  [ -e $logdir ]||mkdir $logdir
  for samp in $na $ecoli;do
    echo $samp
    fq1=$(find $trimdir -name "${samp}_*1.fq.gz")
    fq2=$(find $trimdir -name "${samp}_*2.fq.gz")
    tmpdir=${tmproot}/align/$samp
    [ -e $tmpdir ]||mkdir -p $tmpdir
    log=$logdir/$samp.align.log
    cell=${samp%-*}
    if [ "$cell" == "Ecoli" ];then
      refdir=$eref
    elif [ "$cell" == "NA12878" ];then
      refdir=$nref
    fi
    bismark --version > $log
    bismark --bam --non_directional --bowtie2 --un -p 2 \
      --genome $refdir \
      -1 $fq1 \
      -2 $fq2 \
      --temp_dir $tmpdir \
      --output_dir $bamdir \
      -B $samp &>> $log
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
extractdir=$outdir/extract
[ -e $extractdir ]||mkdir $extractdir
if [ "$1" == "extract" ];then
  logdir=$logroot/extract
  [ -e $logdir ]||mkdir $logdir
  for bam in $(find $bamdir -name "*nodupqsort.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    cell=${samp%-*}
    if [ "$cell" == "Ecoli" ];then
      refdir=$eref
    elif [ "$cell" == "NA12878" ];then
      refdir=$nref
    fi
    log=$logdir/$samp.extract.log
    args="--ignore 2 --ignore_r2 2 --ignore_3prime 1 --ignore_3prime_r2 1"
    com="bismark_methylation_extractor -p --multicore 8 --gzip \
      --genome_folder $refdir --no_header $args \
      $bam -o $extractdir &> $log"
    echo $com
    eval $com
  done
fi
convertdir=$outdir/convert
[ -e $convertdir ]||mkdir $convertdir
if [ "$1" == "convert" ];then
  logdir=$logroot/convert
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/{}.convert.log
  for bam in $(find $bamdir -name "*nodupqsort.bam");do
    samp=$(basename "$bam")
    samp=${samp%%.*}
    samps="$samps $samp"
  done
  com="bismark2bedGraph --counts --CX --dir $convertdir \
    -o {}.bedGraph $extractdir/*{}.*txt.gz &> $log"
  parallel $com ::: $samps
fi
freqdir=$outdir/mfreq
[ -e $freqdir ]||mkdir $freqdir
if [ "$1" == "getfreq" ];then
  cov=$convertdir/{}.bismark.cov.gz
  logdir=$logroot/methfreq
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/{}.methfreq.log
  com="coverage2cytosine --NOMe-seq --gzip --dir $freqdir \
    --genome_folder $nref \
    -o {} $cov &> $log"
  parallel "$com" ::: $na
  com="coverage2cytosine --NOMe-seq --gzip --dir $freqdir \
    --genome_folder $eref \
    -o {} $cov &> $log"
  parallel "$com" ::: $ecoli
fi

if [ "$1" == "summarize" ];then
  summary=$freqdir/180922_methylation_summary.txt
  [ ! -e $summary ]||rm $summary
  for samp in $na $ecoli;do
    for mod in CpG GpC;do
      report=$freqdir/$samp.NOMe.${mod}_report.txt.gz
      sum=$( gunzip -c $report |\
        awk '{ meth+=$4;unmeth+=$5 }END{ print meth,unmeth,meth/(meth+unmeth) }')
      echo "$samp $mod $sum" | tr " " "\t" >> $summary
    done
  done
fi
  
reportdir=$outdir/report
[ -e $reportdir ]||mkdir $reportdir
if [ "$1" == "report" ];then
  find $bamdir -name "*txt" -exec cp \{\} $extractdir \;
  for samp in $na $ecoli;do
    pre=$extractdir/$samp
    log=$reportdir/$samp.report.log
    bismark2report --dir $reportdir \
      --alignment_report ${pre}_PE_report.txt \
#      --dedup_report ${pre}.dupmetrics.txt \
      --splitting_report ${pre}.*splitting_report.txt \
      --mbias_report ${pre}.*M-bias.txt 2> $log
  done
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


