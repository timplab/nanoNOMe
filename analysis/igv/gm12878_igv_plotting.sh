#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
bed=$root/annotations/breastcancer/bcan_10a_vs_231_repeats.bed
reg=bcan_sig_repeats
outdir=$root/igv
[ -e $outdir ]||mkdir $outdir

cell=GM12878
echo "$cell bsseq"
bam=$(find $root/data/gm12878 -name "GM12878*.bam")
out=$outdir/${cell}_BSseq_$reg.bam
echo $out
samtools view -@ 10 -L $bed -hb $bam > $out
samtools index $out

echo $cell
bam=$root/pooled/bam/$cell.pooled.bam
mbeddir=$root/pooled/methylation/methbyread_all
cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
script=$srcdir/../../script/convertBam.py
log=$outdir/${cell}_$reg.log
out=$outdir/${cell}_$reg.bam

com="python -u $script -v -t 10 \
  -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
  samtools sort -o $out"
echo $com
eval $com
samtools index $out


