#!/bin/bash

# get TSS from gene bed file (sorted)

infile=$1

awk '{ OFS="\t" }{ if($6=="+"){ $2=$2;$3=$2+1 }else{ $3=$3;$2=$3-1 } print }' $infile
