#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38
regs=$root/hg38_genes.TSS.2000bp.bed
shuffle=$root/hg38_genes.TSS.2000bp.shuffle.bed
genome=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa.fai

bedtools shuffle -i $regs -g $genome |\
  sort -k1,1 -k2,2n > $shuffle
