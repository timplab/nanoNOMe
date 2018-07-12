#!/bin/bash
codedir=`readlink -f $(dirname "$0")`

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
    --tar ) # tar dir
      tardir=$2
      shift 2
      ;;
    --all ) # for fqcat : cat all of a sample into one file
      catarg="--all"
      shift 1
      ;;
    * ) break
      ;;
  esac
done

if [ "$env" == "marcc" ];then
  sbarg="-p shared"
fi

logroot=$outdir/log
[ -e $logroot ]||mkdir $logroot

if [ "$module" == "untar" ];then
  tar=`find $tardir -name "*$base*.t*gz"`
  if [ ! $tar ];then
    echo "no tar, exiting"
    exit
  fi
  rawdir=$root/fast5/$base
  [ -e $rawdir ]||mkdir -p $rawdir
  batchcom=$codedir/../../slurm/batchcommand.scr
  logdir=$logroot/untar
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/$base.untar.log
  if [ "$env" == "marcc" ];then
    sbarg="$sbarg -J untar -e $log -o $log"
    com="sbatch $sbarg $batchcom tar -vxzf $tar -C $rawdir"
  elif [ "$env" == "aws" ];then
    com="tar -vxzf $tar -C $rawdir > $log"
  fi
  echo $com
  eval $com
fi

if [ "$module" == "basecall" ];then
  rawdir=$root/fast5/$base
  calldir=$root/bcall/$base
  [ -e $calldir ]||mkdir -p $calldir
  configdir=`readlink -f $codedir/../config`
  if [ "$env" == "marcc" ];then
    config=$configdir/albacore_config_sv_marcc.cfg
  elif [ "$env" == "aws" ];then
    config=$configdir/albacore_config_sv_aws.cfg
  fi
  logdir=$logroot/bcall
  [ -e $logdir ]||mkdir -p $logdir
  log=$logdir/$base.bcall.log
  echo "$codedir/call_wrapper.sh -i $rawdir -o $calldir \
   -s $codedir -e $env -a "--config $config" &> $log"
  $codedir/call_wrapper.sh -i $rawdir -o $calldir \
   -s $codedir -e $env -a "--config $config" &> $log
fi

if [ "$module" == "fqcat" ];then
  calldir=$root/bcall/$base
  if [ "$catarg" == "--all" ];then
    fqdir=$outdir/fastq
  else
    fqdir=$outdir/fastq/$base
    [ -e $fqdir ]||mkdir -p $fqdir
  fi
  prefix=$fqdir/$base
  logdir=$logroot/fqcat
  [ -e $logdir ]||mkdir $logdir
  log=$logdir/$base.fqcat.log
  batchcom=$codedir/../../slurm/batchcommand.scr
  sbarg="$sbarg -J fqcat -e $log -o $log"
  com="sbatch $sbarg $batchcom $codedir/cat_bcall.sh -p $prefix -i $calldir $catarg"
  echo $com
  eval $com
fi
