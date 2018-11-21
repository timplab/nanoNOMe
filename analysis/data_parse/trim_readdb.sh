#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled
dbdir=$root/fastq
idxdir=$root/index
cell=GM12878

dbs=$(find $dbdir -name "*$cell*readdb")
out=$idxdir/$cell.readname_to_filename.txt

cat $dbs |\
  awk '{ gsub(".*/","",$2); print }' > $out
