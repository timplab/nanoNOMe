#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
if [ ! -z $1 ];then
  root=$1
fi
dir=$root/database/gm12878/dnase
samp=GM12878
bigwig=$dir/*bigWig
converter=../../util/bigWigToBedGraph

bed=$dir/${samp}_dnase_signal.bedGraph
if [ ! -e $bed ];then
  $converter $bigwig $bed
  bgzip $bed;tabix -p bed $bed.gz
fi

