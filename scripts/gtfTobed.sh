#!/bin/bash

# changing to 0-based start, 1-based end

infile=$1
ann=$2
ann=${ann:=gene}

pre="gunzip -c $1 | grep -v \"#\""
post="tr -d '\";' | sort -k1,1 -k2,2n"
if [ "$ann" == "gene" ];then
  eval $pre | awk '{ OFS="\t" }{ if($3=="gene") print $1,$4-1,$5,$10,".",$7,$14 }' |\
    eval $post
elif [ "$ann" == "transcript" ];then
  eval $pre | awk 'OFS="\t"{ if($3=="transcript") print $1,$4-1,$5,$12,".",$7,$10,$16 }' |\
    eval $post
fi

