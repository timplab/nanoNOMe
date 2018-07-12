#!/bin/bash

#argument parsing
while :
do
  case "$1" in 
    -i | --input) #path of the raw fast5 folders
      rawpath=`readlink -f $2`
      shift 2
      ;;
    -o | --outdir) # output dir
      outdir=$2
      shift 2
      ;;
    -a | --args) # arguments to basecaller
      args="$2"
      shift 2
      ;;
    -s | --src) #source code path
      srcpath=$2
      shift 2
      ;;
    --listpath) ##optional
      listpath=$2
      shift 2
      ;;
    -e | --environment ) #slurm env
      env=$2
      shift 2
      ;;
    *) break
      ;;
  esac
done

if [ "$env" == "marcc" ];then
  slurmopt="--partition=shared --time=4:0:0"
  slurmopt="${slurmopt} -c 24"
  slurmtype="--marcc"
elif [ "$env" == "aws" ];then
  slurmopt="-c 36"
  slurmtype="--aws"
fi

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

if [ $listpath ];then
  echo "path list provided in $listpath"
  arraylist=`cat $listpath`
else
  onepath=`find $rawpath -name "*.fast5" -type f -print -quit`
  onedir=$(dirname "$onepath")
  dir=$(dirname "$onedir")
  arraylist=`printf "%s\n" $dir/*`
fi
n=`echo $arraylist | wc -w`
array="1-$n"

#check
#echo "$arraylist"
#echo $rawpath

echo "following command is used:"
echo "sbatch -a $array -D $outdir $slurmopt \
  $srcpath/oxford_call.scr -i ${rawpath} -o $outdir \
  --arrayval "$arraylist" -a "$args" $slurmtype"
sbatch -a $array -D $outdir $slurmopt \
  $srcpath/oxford_call.scr -o $outdir \
  --arrayval "$arraylist" -a "$args" $slurmtype
  
