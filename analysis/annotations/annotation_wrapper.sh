#!/bin/bash
# this is a wrapper for downloading and parsing annotations
# you must change the root to the annotations 
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis/database"

./annotation_download.sh "$root"
./parse_repeats.sh "$root/hg38/hg38_repeats.txt.gz"
./parse_genes.sh "$root/hg38/hg38_genes.gtf.gz"
