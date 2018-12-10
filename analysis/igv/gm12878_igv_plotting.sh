#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
outdir=$root/igv
[ -e $outdir ]||mkdir $outdir
# change these
#reg=ctcf_paternal
#bed=$root/pooled/hap/hapcut_ase/GM12878.pooled.ase.intersected.CTCF.allcomp.center.noTSS.2000bp.paternal.bed
#bam=$root/pooled/hap/hapcut_ase/GM12878.pooled.ase.intersected.CTCF.allcomp.center.noTSS.2000bp.paternal.bam
reg=ctcf
bed=$root/annotations/gm12878/GM12878_CTCF_ctcfbsdb_allcomp.center.noTSS.2000bp.bed
cell=GM12878
bam=$root/pooled/bam/$cell.pooled.bam

if [ 1 -eq 0 ];then
  cell=GM12878
  echo "$cell bsseq"
  bam=$(find $root/data/gm12878 -name "GM12878*.bam")
  out=$outdir/${cell}_BSseq_$reg.bam
  echo $out
  samtools view -@ 10 -L $bed -hb $bam > $out
  samtools index $out
fi

echo $cell
mbeddir=$root/pooled/methylation/methbyread_all
cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
script=$srcdir/../../script/convertBam.py
log=$outdir/${cell}_$reg.log
out=$outdir/${cell}_$reg.bam

com="python -u $script -v -t 10 \
  -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
  samtools sort -o $out"
com="head -n100 $bed | python -u $script -v -t 10 \
  -b $bam -c $cpg -g $gpc 2> $log |\
  samtools sort -o $out"
echo $com
eval $com
samtools index $out


