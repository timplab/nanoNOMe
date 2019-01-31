#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis
dir=$root/data/bsseq/mfreq
gs=$root/data/hg38/hg38_genomesize.txt
parser=~/Code/nanoNOMe/script/makeWig.py
samps="GM12878_BSseq_1"

for samp in $samps; do
  for mod in cpg gpc; do
    pre=$samp.$mod
    prefixes="$prefixes $pre"
  done
done

data=$dir/{}.methfreq.txt.gz
com="python $parser -i $data -v \
  -o $wigdir/{}.methylation.wig \
  -c $wigdir/{}.coverage.wig"
echo $com
parallel "$com" ::: $prefixes &> $log
    
for pre in $prefixes;do
  wigpre="$wigpre $pre.methylation $pre.coverage"
done
echo $wigpre
# make bigwig
parallel "wigToBigWig $wigdir/{}.wig $gs $wigdir/{}.bw" ::: $wigpre

#reg="chr1:1000000-1200000"
#
#out=$dir/$samp.$mod.wig
#echo "track type=wiggle_0" > $out
#echo "variableStep chrom=chr1" >> $out
#tabix $data $reg | awk 'OFS="\t"{ if($4+$5>2) print $2,$4/($4+$5)}' >> $out
#echo "variableStep chrom=chr2" >> $out
#tabix $data $reg | awk 'OFS="\t"{ if($4+$5>2) print $2,$4/($4+$5)}' >> $out
