#!/bin/bash
root=/dilithium/Data/Nanopore/projects/gm12878/analysis
bedscr=~/Code/nanoNOMe/nanopolish/mtsv2bedGraph.py
mbeddir=$root/mbed
logdir=$root/log
tmproot=$root/tmp
[ -e $mbeddir ]||mkdir $mbeddir
[ -e $logdir ]||mkdir $logdir
[ -e $tmproot ]||mkdir $tmproot
mod=gpc
mcalldir=$root/mcall-$mod

if [ "$1" == "mbed" ];then
  bed=$mbeddir/{}.$mod.meth.bed.gz
  tsvs=$(find $mcalldir -name "*tsv.gz" -type f)
  for tsv in $tsvs;do
    base=$(basename "$tsv")
    samp=${base%.$mod*}
    ext=${base#$samp}
    bases="$bases $samp"
    tmpdir=$tmproot/$samp
    [ -e $tmpdir ]||mkdir $tmpdir
  done
  subcom="gunzip -c $mcalldir/{}$ext |\
    python $bedscr -m $mod |\
    sort -T $tmproot/{} -k1,1 -k2,2n |\
    bgzip > $bed && tabix -p bed $bed"
  log=$logdir/gm12878.mbed.log
  parallel "$subcom" ::: $bases 2> $log

  cd $mbeddir
  zcat *gz | sort -T ../pooled/ -k1,1 -k2,2n | bgzip > ../pooled/GM12878_wgs.gpc.meth.bed.gz
fi

mfreqdir=$root/mfreq
[ -e $mfreqdir ]||mkdir $mfreqdir
mfreq=$mfreqdir/GM12878_wgs.$mod.methfreq.txt.gz
parser=../../script/parseMethylbed.py
if [ "$1" == "mfreq" ];then
  bed=$(find $root/pooled/ -name "*bed.gz")
  com="python $parser frequency -i $bed -m $mod |\
    bgzip > $mfreq && tabix -b 2 -e 2 $mfreq"
  echo $com
  eval $com
fi
