#!/bin/bash
dir=/home/isac/Dropbox/Data/nome-seq/db/ctcf
inbed=$dir/nome-seq_CTCF_coordinates.bed
chain=/mithril/Data/NGS/Reference/human/liftover/hg19ToHg38.over.chain

liftOver $inbed $chain $dir/ctcf_hg38.bed $dir/ctcf_unmapped.bed

# sort
sort -k1,1 -k2,2n $dir/ctcf_hg38.bed > $dir/ctcf_hg38.sorted.bed
