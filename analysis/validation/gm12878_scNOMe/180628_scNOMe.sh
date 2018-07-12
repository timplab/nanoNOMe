#!/bin/bash
# from : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5487215/
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/validation/scNOMe
logroot=$root/log
tmproot=$root/tmp
covdir=$root/covreport
freqdir=$root/methfreq
meta="/home/isac/Dropbox (Timp Lab)/Data/nome-seq/validation/gm12878 scNOMe/180628_scNOMe_data.txt"
metainfo=$(awk 'OFS="\t"{ print $0,0,0 }' "$meta")
refdir=/mithril/Data/NGS/Reference/hg38_noalt

rawdir=$root/raw
#if [ "$1" == "datadl_cov" ];then
#  # deprecated?
#  echo "downloading data"
#  linkname="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2220nnn/accession/suppl/accession_rep_CpG_clean.cov.txt.gz"
#  log=$logroot/dl.log
#  [ -e $log ]&&rm $log
#  while IFS=$'\t' read -r -a line ;do
#    fileid=${line[0]}
#    samp=${line[1]}
#    treat=$(echo "${line[3]}"|tr "?" " " )
#    address=${linkname//accession/"$fileid"}
#    address=${address//rep/"$samp"}
#    for mod in CpG GpC;do
#      samppath="$covdir/${line[0]}.$mod.cov.txt.gz"
#      echo $samppath
#      address=${address//CpG/"$mod"}
#      wget "$address" -O $samppath 2>> $log
#    done
#  done < "$meta"
#fi

meta="/home/isac/Dropbox (Timp Lab)/Data/nome-seq/validation/gm12878 scNOMe/scNOMe_data.txt"
if [ "$1" == "datadl" ];then
  cd $rawdir 
  sras=$(awk 'FS="\t"{ if(NR>1) print $8 }' "$meta" | tr "\n" " ")
  echo $sras
  parallel fastq-dump --split-3 --gzip {} ::: $sras
fi

fqdir=$root/fastq
if [ "$1" == "rename" ];then
  awk 'NR>1' "$meta" | while IFS=$'\t' read -r -a line ;do
    cell=${line[2]}
    lab=${line[4]}
    rep=${line[5]}
    sra=${line[6]}
    for pair in 1 2;do
      fq=$(find $rawdir -name "${sra}_$pair.fastq.gz")
      newname=$fqdir/${cell}_${lab}_rep${rep}_R$pair.fastq.gz
      mv $fq $newname
    done
  done
fi

trimdir=$root/fqtrim
[ -e $trimdir ]||mkdir $trimdir
if [ "$1" == "fqtrim" ];then
  for fq in $(find $fqdir -name "*R1.fastq.gz");do
    samp=$(basename "$fq")
    samp=${samp%%_R1*}
    samples="$samples $samp"
  done
  fq1=$fqdir/{}_R1.fastq.gz
  fq2=$fqdir/{}_R2.fastq.gz
  log=$trimdir/{}.trim.log
  com="trim_galore -q 28 --paired \
    $fq1 $fq2 -o $trimdir &> $log"
  parallel $com ::: $samples
fi

bamdir=$root/bam
[ -e $bamdir ]||mkdir $bamdir
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
