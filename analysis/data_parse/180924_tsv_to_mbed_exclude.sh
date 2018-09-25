#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/gm12878/ngmlr
bedscr=~/Code/nanoNOMe/nanopolish/mtsv2bedGraph.py
mbeddir=$root/mbed_exclude
logdir=$mbeddir/log
tmproot=$mbeddir/tmp
[ -e $mbeddir ]||mkdir $mbeddir
[ -e $logdir ]||mkdir $logdir
[ -e $tmproot ]||mkdir $tmproot
mods="cpg gpc"

if [ "$1" == "mbed" ];then
  for mod in $mods;do
    mcalldir=$root/mcall-$mod
    outdir=$mbeddir/$mod
    [ -e $outdir ]||mkdir $outdir
    bed=$outdir/{}.$mod.meth.bed.gz
    tsvs=$(find $mcalldir -name "*tsv" -type f)
    bases=""
    for tsv in $tsvs;do
      base=$(basename "$tsv")
      samp=${base%.$mod*}
      ext=${base#$samp}
      bases="$bases $samp"
      tmpdir=$tmproot/$samp
      [ -e $tmpdir ]||mkdir $tmpdir
    done
    if [ "$mod" == "cpg" ];then
      exc="GC"
    elif [ "$mod" == "gpc" ];then
      exc="CG"
    fi
    subcom="python $bedscr -m $mod -e $exc -i $mcalldir/*/{}$ext |\
      sort -T $tmproot/{} -k1,1 -k2,2n |\
      bgzip > $bed && tabix -p bed $bed"
    log=$logdir/gm12878.$mod.mbed.log
    echo $subcom
    parallel "$subcom" ::: $bases 2> $log
  done
fi
if [ "$1" == "pool" ];then
  outdir="$mbeddir/{}"
  pooled="$mbeddir/GM12878_nobleedthrough.{}.meth.bed.gz"
  com="gunzip -c $outdir/*gz | sort -T $outdir -k1,1 -k2,2n |\
    bgzip > $pooled &&\
    tabix -p bed $pooled"
  log=$logdir/gm12878.pool.log
  parallel "$com" ::: $mods 2> $log

fi

mfreqdir=$mbeddir
[ -e $mfreqdir ]||mkdir $mfreqdir
mfreq=$mfreqdir/GM12878_nobleedthrough.{}.methfreq.txt.gz
parser=../../script/parseMethylbed.py
if [ "$1" == "mfreq" ];then
  com="python $parser frequency -i $mbeddir/GM12878_nobleedthrough.{}.meth.bed.gz -m {} |\
    bgzip > $mfreq && tabix -b 2 -e 2 $mfreq"
  echo $com
  parallel "$com" ::: $mods
fi
