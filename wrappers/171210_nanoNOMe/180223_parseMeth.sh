#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171210_nomeseq
mcalldir=$root/methcall
beddir=$root/bed
[ -e $beddir ]||mkdir $beddir
scriptdir=~/Code/ilee/dnamods/scripts

if [ "$2" == "gpc" ];then
  mod=gpc
  rec=GC
elif [ "$2" == "cpg" ];then
  mod=cpg
  rec=CG
fi

samp=171210_gpc_nome_mda231_200U
base=$samp.$mod
methcall=$mcalldir/${base}-hg38.meth.tsv.gz
mbed=$beddir/${base}.meth.bedGraph
sortbed=$beddir/${base}.meth.sorted.bedGraph
freq=$beddir/${base}.meth.freq.bedGraph
sortfreq=$beddir/${base}.meth.freq.sorted.bedGraph.gz


if [ "$1" == "makebed" ];then
  gunzip -c $methcall |\
    $scriptdir/methylationToBed.py -c 2.5 -m $rec \ 
    > $mbed
fi

if [ "$1" == "sortbed" ];then
  sort -k1,1 -k2,2n $mbed > $sortbed
fi

if [ "$1" == "frequency" ];then
  $scriptdir/parseMethylbed.py -t frequency -i $sortbed -m $rec > $freq
fi

if [ "$1" == "sortfreq" ];then
  gunzip -c $freq | sort -k1,1 -k2,2n | gzip > $sortfreq
fi

