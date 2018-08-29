#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
dir=$root/pooled/methylation/mfreq_all
cell=GM12878
cell=GM12878_BSseq_ENCLB794YYH
mfreq=$dir/$cell.gpc.methfreq.txt.gz
parser=../../nanopolish/parseMethylFreq.py
reg=/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38/hg38_repeats_LINE.bed
out=$root/tmp.bed

head $reg -n1000 | python -u $parser by-region -v -t 12 -i $mfreq #-o $out



