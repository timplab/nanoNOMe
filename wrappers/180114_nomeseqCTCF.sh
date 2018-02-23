#!/bin/bash
outroot=/home/isac/Dropbox/Data/nome-seq

# liftover
dbdir=$outroot/db
if [ 0 -eq 1 ];then
  echo "performing liftover"
  ctcfh19=$dbdir/nome-seq_CTCF_coordinates.bed
  chain=/mithril/Data/NGS/Reference/human/liftover/hg19ToHg38.over.chain
  ctcfhg38=$dbdir/nome-seq_CTCF_coordinates_hg38.bed
  unmapped=$dbdir/CTCF_liftover_unmapped.bed
  liftOver $ctcfh19 $chain $ctcfhg38 $unmapped
fi

# R code for plotting
if [ 1 -eq 1 ];then
  Rscript 180114_nomeseqCTCF.R
fi
