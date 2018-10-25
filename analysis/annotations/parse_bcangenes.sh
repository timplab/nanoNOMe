#!/bin/bash
dir=../../annotations/breast_cancer
genes=../../annotations/hg38/hg38_genes.bed

# convert to bed format
bcanbed=$dir/breast_cancer_genes.bed
names=$(awk 'NR>1{ print $1 }' $dir/*txt)
grep -F "$names" $genes > $bcanbed
