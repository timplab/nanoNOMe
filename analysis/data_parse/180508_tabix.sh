#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation
mfreqdir=$root/mfreq
parser=/home/isac/Code/ilee/nanopolish/script/parseMethylbed.py

for cell in GM12878 MCF10A MCF7 MDAMB231;do
  echo $cell
  for mod in cpg gpc;do
    echo $mod
    base=$cell.$mod
    infile=$mfreqdir/$base.methfreq.txt
    echo $infile
    if [ "$1" == "bgzip" ];then
      bgzip $infile
    fi
    if [ "$1" == "tabix" ];then
      tabix -b 2 -e 2 $infile.gz
    fi
  done
done

