#!/bin/bash
env=aws
root=/shared/data
tarroot=$root/tar
rawroot=/shared/data/raw
outroot=$root/analysis
metadata=$outroot/nome-seq_data.csv
bcall=`readlink -f ../../../oxford/slurm/bcallWrapper.sh`  
align=../../../oxford/slurm/alignWrapper.sh
aligner=ngmlr
nwrapper=`readlink -f ../../../nanopolish/slurm/nanopolishWrapper.sh`
ref=$root/Reference/hg38_noalt/hg38_noalt.fa
if [ "$1" == "yeast" ];then
  echo yeast
  ref=$root/Reference/yeast/sacCer3.fa
fi
nanopolish=/shared/Code/nanopolish/nanopolish

awk 'NR>1' $metadata | while IFS=$',' read -r -a line
do
  samp=${line[1]}
  cell=${line[2]}
  tarflag=${line[8]}
  callflag=${line[9]}
  catflag=${line[10]}
  indexflag=${line[11]}
  alignflag=${line[12]}
  cpgflag=${line[13]}
  gpcflag=${line[14]}
  cpgmbedflag=${line[15]}
  gpcmbedflag=${line[16]}
  rawdir=$rawroot/$cell
  outdir=$outroot/$cell/ngmlr
  tardir=$tarroot/fast5/$cell
  [ -e $outdir ]||mkdir -p $outdir
  [ -e $rawdir ]||mkdir -p $rawdir
  if [ "$tarflag" == "n" ];then
    echo "untarring $samp"
    com="$bcall untar -e $env --tar $tardir -i $rawdir -o $outdir -b $samp"
    eval $com
    continue
  fi
  if [ "$callflag" == "n" ];then
    echo "bascalling $samp"
    $bcall basecall -e $env -i $rawdir -o $outdir -b $samp
    continue
  fi
  if [ "$catflag" == "n" ];then
    echo "concatenating $samp fastq files"
    $bcall fqcat -e $env -i $rawdir -o $outdir -b $samp
    continue
  fi
  if [ "$indexflag" == "n" ];then
    echo "indexing $samp"
    $nwrapper index -e $env -i $rawdir -o $outdir -b $samp -n $nanopolish
  fi
  if [ "$alignflag" == "n" ];then
    echo "aligning $samp"
    $align -e $env -d $outdir -b $samp -a $aligner -r $ref
    continue
  fi
  if [ "$cpgflag" == "n" ];then
    echo "methyation calling $samp"
    $nwrapper mcall -e $env -o $outdir -b $samp -n $nanopolish -g $ref -m cpg
  fi
  if [ "$gpcflag" == "n" ];then
    echo "methyation calling $samp"
    $nwrapper mcall -e $env -o $outdir -b $samp -n $nanopolish -g $ref -m gpc
    continue
  fi
done


