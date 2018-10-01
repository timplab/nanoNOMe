#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
methdir=$root/pooled/methylation/methbyread_all
dbdir=$root/database/gm12878/ctcf
cell=GM12878
db=$(find $dbdir -name "${cell}_ctcf.2000bp.bed")
cpg=$(find $methdir -name "$cell.cpg*bed.gz")
out=$root/intersect/$cell.ctcf.readnames.txt

bedtools intersect -a $cpg -b $db | cut -f4 > $out
