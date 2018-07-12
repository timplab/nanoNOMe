#!/bin/bash
codedir=`readlink -f $(dirname "$0")`
alignscr=$codedir/oxford_align.scr

while :
do
  case "$1" in
    -e | --environment ) # slurm environment
      env=$2
      shift 2
      ;;
    -d | --dir ) # root directory for data
      dir=$2
      shift 2
      ;;
    -b | --base ) # basename for sample
      base=$2
      shift 2
      ;;
    -a | --aligner ) # aligner
      aligner=$2
      shift 2
      ;;
    -r | --reference ) #reference fa
      ref=$2
      shift 2
      ;;
    * ) break
      ;;
  esac
done

if [ "$env" == "marcc" ];then
  sbarg="-p shared -t 12:0:0 -c 24"
elif [ "$env" == "aws" ];then
  sbarg="-c 36"
fi

fqs=`find $dir/fastq -name "$base.*f*q.gz"`
fqnum=`echo $fqs | wc -w`
echo $fqnum
if [ "$fqnum" -eq 1 ];then
  echo "only one fastq for this sample"
  # one fastq for this sample
  logdir=$dir/log/align
  [ -e $logdir ]||mkdir -p $logdir
  log=$logdir/$base.align.log
  fq=$fqs
  outdir=$dir/bam
  [ -e $outdir ]||mkdir -p $outdir
  com="sbatch $sbarg -e $log -o $log \
    $alignscr -r $ref -i $fq -a $aligner \
    -o $outdir -e $env --samarg \"-q 20\""
  echo $com
  eval $com
  exit
fi

logdir=$dir/log/align/$base
[ -e $logdir ]||mkdir -p $logdir
fqprefix=$dir/fastq/$base/$base
outdir=${dir}/bam/$base
[ -e $outdir ]||mkdir -p $outdir
log="$logdir/$base.%a.align.log"
arraynum=`echo $fqs | wc -w`
array="0-$((${arraynum}-1))"
echo "sbatch $sbarg -e $log -o $log -a $array \
  $alignscr -r $ref -i $fqprefix -a $aligner \
  -o $outdir -e $env --samarg "-q 20""
sbatch $sbarg -e $log -o $log -a $array \
  $alignscr -r $ref -i $fqprefix -a $aligner \
  -o $outdir -e $env --samarg "-q 20"
