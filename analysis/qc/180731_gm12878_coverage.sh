#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
datadir=$root/pooled/methylation/methbyread
dbdir=$root/database
covdir=$root/coverage
plotdir=$root/plots/coverage

cell=GM12878

dbname="$1"
# genes
if [ "$dbname" == "genes" ];then
  db=$dbdir/hg38/hg38_genes.bed
elif [ "$dbname" == "promoter" ];then
  db=$dbdir/hg38/hg38_genes.TSS.2000bp.bed
elif [ "$dbname" == "enhancer" ];then
  db=$dbdir/hg38/hg38_enhancers.bed
fi

data=$datadir/$cell.{}.pooled.methbyread.bed.gz
labs="cpg gpc"
if [ "$2" == "getcov" ];then
  out=$covdir/$cell.{}.$dbname.bedcov.bed
  parallel "bedtools coverage -sorted -hist -a $db -b $data > $out" ::: $labs
fi
if [ "$2" == "plothist" ];then
  plotout=$plotdir/$cell.$dbname.covhist.pdf
  for mod in $labs;do
    histdat="$histdat $covdir/$cell.$mod.$dbname.bedcov.bed"
  done
  Rscript plotCoverageHistogram.R $histdat $plotout
fi
