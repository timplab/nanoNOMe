#!/bin/bash
codedir=`readlink -f $(dirname "$0")`
indexscr=$codedir/fast5index.scr
mcallscr=$codedir/callmeth.scr
batchscr=`readlink -f ${codedir}/../../slurm/batchcommand.scr`

module=$1
shift 1

while :
do
  case "$1" in
    -e | --environment ) # slurm environment
      env=$2
      shift 2
      ;;
    -i | --indir ) # root directory of raw tarball
      root=$2
      shift 2
      ;;
    -o | --outdir ) # root directory for output
      outdir=$2
      shift 2
      ;;
    -b | --base ) # basename for sample
      base=$2
      shift 2
      ;;
    -n | --npath ) #nanopolish path
      npath=$2
      shift 2
      ;;
    -m | --mod ) # meth calling mod
      mod=$2
      shift 2
      ;;
    -g | --genome ) # genome for meth calling
      ref=$2
      shift 2
      ;;
    --extra ) # extra arguments - for generating methylation bed files
      extarg=$2
      shift 2
      ;;
    * ) break
      ;;
  esac
done

if [ "$env" == "marcc" ];then
  sbarg="-p shared -t 24:0:0"
  cores=24
elif [ "$env" == "aws" ];then
  cores=36
fi

logroot=$outdir/log
[ -e $logroot ]||mkdir $logroot

if [ "$module" == "index" ];then
  fqs=`find $outdir/fastq -name "$base*fastq.gz"`
  fqnum=`echo $fqs | wc -w`
  if [ $fqnum -eq 1 ];then
    echo "only one fq for $base"
    rawdir=`readlink -f $root/fast5/$base`
    logdir=$logroot/index
    [ -e $logdir ]||mkdir $logdir
    log=$logdir/$base.index.log
    com="sbatch $sbarg -e $log -o $log \
      $indexscr -s $npath -d $rawdir -i $fqs"
    echo $com
    eval $com
    exit
  fi
  onepath=`find $root/fast5/$base -name "*.fast5" -type f -print -quit`
  onedir=$(dirname "$onepath")
  rawdir=$(dirname "$onedir")
  logdir=$logroot/index/$base
  [ -e $logdir ]||mkdir -p $logdir
  log="$logdir/$base.%a.index.log"
  fqprefix=$outdir/fastq/$base/$base
  arraynum=`echo $fqs | wc -w`
  array="0-$((${arraynum}-1))"
  com="sbatch $sbarg -e $log -o $log -a $array \
    $indexscr -s $npath -d $rawdir -i $fqprefix"
  echo $com
  eval $com
fi

if [ "$module" == "mcall" ];then
  fqs=`find $outdir/fastq -name "$base.*f*q.gz"`
  fqnum=`echo $fqs | wc -w`
  if [ $fqnum -eq 1 ];then
    logdir=$logroot/mcall-${mod}
    [ -e $logdir ]||mkdir -p $logdir
    log="$logdir/$base.mcall.log"
    fq=$fqs
    bam=$(find $outdir/bam -name "$base.sorted.bam")
    mcalldir=$outdir/mcall-$mod
    [ -e $mcalldir ]||mkdir -p $mcalldir
    mcall="$mcalldir/$base.$mod.meth.tsv"
    com="sbatch $sbarg -c $cores -e $log -o $log \
      $mcallscr -s $npath -i $fq -g $ref -b $bam \
      -m $mod -e $env -o $mcall"
    echo $com
    eval $com
    exit
  fi
  logdir=$logroot/mcall-${mod}/$base
  [ -e $logdir ]||mkdir -p $logdir
  log="$logdir/$base.%a.mcall.log"
  fq="$outdir/fastq/$base/$base.%a.fastq.gz"  
  bam="$outdir/bam/$base/$base.%a.sorted.bam"
  mcalldir=$outdir/mcall-${mod}/$base
  [ -e $mcalldir ]||mkdir -p $mcalldir
  mcall="$mcalldir/$base.%a.$mod.meth.tsv"
  arraynum=`echo $fqs | wc -w`
  array="0-$((${arraynum}-1))"
  com="sbatch $sbarg -c $cores -e $log -o $log -a $array \
    $mcallscr -s $npath -i $fq -g $ref -b $bam \
    -m $mod -e $env -o $mcall"
  echo $com
  eval $com
fi

if [ "$module" == "mbed" ];then
  scriptdir=`readlink -f $codedir/..`
  callnum=$(find $outdir/mcall-$mod -name "$base.*$mod.meth.tsv" | wc -w)
  if [ $callnum -eq 1 ];then
    logdir=$logroot/mbed-$mod
    [ -e $logdir ]||mkdir -p $logdir
    log="$logdir/$base.$mod.mbed.log"
    mcall=$(find $outdir/mcall-$mod -name "$base.*$mod.meth.tsv")
    echo $mcall
    beddir=$outdir/mbed
    [ -e $beddir ]||mkdir -p $beddir
    bed=$beddir/$base.$mod.meth.bed.gz
    subcom="python $scriptdir/mtsv2bedGraph.py -i $mcall -m $mod $extarg |\
      sort -T $beddir -k1,1 -k2,2n |\
      bgzip > $bed && tabix -p bed $bed"
    com="sbatch $sbarg -J mbed -e $log -o $log $batchscr $subcom"
    echo $com
    $com
    exit
  fi
  logdir=$logroot/mbed-$mod/$base
  [ -e $logdir ]||mkdir -p $logdir
  log="$logdir/$base.%a.$mod.mbed.log"
  mcall="$outdir/mcall-${mod}/$base/$base.%a.$mod.meth.tsv"
  arraynum=`find $outdir/mcall-${mod}/$base -name "$base*$mod.meth.tsv" | wc -w`
  array="0-$((${arraynum}-1))"
  tmpdir=$outdir/tmp/mbed-$mod/$base
  [ -e $tmpdir ]||mkdir -p $tmpdir
  tmpout="$tmpdir/$base.%a.$mod.meth.bed"
  pycom="python $scriptdir/mtsv2bedGraph.py -i $mcall -m $mod $extarg > $tmpout"
  com="sbatch $sbarg -J mbed -e $log -o $log -a $array $batchscr $pycom"
  echo $com
  job=$($com)
  echo $job
  jobnumber=${job##*" "}
  tmpbeds="$tmpdir/${base}*${mod}.meth.bed"
  beddir=$outdir/mbed
  [ -e $beddir ]||mkdir -p $beddir
  bedout=$beddir/$base.$mod.meth.bed.gz
  logdir=$logroot/bedmerge
  [ -e $logdir ]||mkdir -p $logdir
  log="$logdir/$base.$mod.merge.log"
  mergecom="cat $tmpbeds |\
    sort -T $tmpdir -k1,1 -k2,2n |\
    bgzip > $bedout && tabix -p bed $bedout"
#  com="sbatch $sbarg -J bedmerge -e $log -o $log $batchscr $mergecom"
  com="sbatch $sbarg -J bedmerge -e $log -o $log --depend=afterok:$jobnumber $batchscr $mergecom"
  echo $com
  $com
fi


