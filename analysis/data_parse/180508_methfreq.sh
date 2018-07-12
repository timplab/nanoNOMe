#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation
mbeddir=$root/methbyread_all
mfreqdir=$root/mfreq_all
[ -e $mfreqdir ]||mkdir $mfreqdir
parser=/home/isac/Code/ilee/nanopolish/script/parseMethylbed.py

for cell in GM12878 MCF10A MCF7 MDAMB231;do
  echo $cell
  for mod in gpc;do
    echo $mod
    base=$cell.$mod
    mbed=`find $mbeddir -name "$base*gz"`
    echo $mbed
    mfreq=$mfreqdir/$base.methfreq.txt.gz
    if [ "$1" == "getfreq" ];then
      com="python $parser frequency -i $mbed -m $mod |\
        bgzip > $mfreq && tabix -b 2 -e 2 $mfreq"
      eval $com
    fi
    if [ "$1" == "getreadlevel" ];then
      python $parser readlevel -i $mbed > $m
    fi
  done
done

