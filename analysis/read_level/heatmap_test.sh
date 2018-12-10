#!/bin/bash
reg=/dilithium/Data/Nanopore/projects/nomeseq/analysis/annotations/breastcancer/MCF10A_vs_MCF7_top_epigenetic_state_genes.TSS.2000bp.bed

head -n1 $reg | python -u /home/isac/Code/nanoNOMe/analysis/bcan/../../script/readlevelHeatmap_test.py -v -w 20 -c 5 -i /dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation/methbyread_all/MCF10A.cpg.pooled.meth.bed.gz -o ~/Dropbox/Data/tmp/heatmatp_test.pdf
