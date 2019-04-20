#!/bin/bash
hg38=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.bed
outdir=/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/hg38/bismap
bedtools subtract -a $hg38 -b $outdir/k100.bismap.bedgraph.gz > $outdir/100.bismap.unmappable.bed
