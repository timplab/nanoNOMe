#!/bin/bash
cells="MCF10A MCF7 MDAMB231"
cells="$cells GM12878"
root=/kyber/Data/Nanopore/projects/nanonome/analysis
dbdir=$root/data/hg38
gs=$dbdir/hg38_genomesize.txt
#regbed=$dbdir/hg38_repeats_LINE_long.bed
#reg=LINE_long
#regbed=$dbdir/hg38_LINE_chr7_L1ME3B.bed
#reg=LINE_chr7_L1ME3B
regbed=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/sv/MDAMB231.sniffles.bed
reg=mda231sv
regbed=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/sv/MCF10A.sniffles.bed
reg=mcf10asv
if [ "$1" == "sine" ];then
  reg=SINE_onlyNanoporebam
  regbed=$dbdir/hg38_${reg}.bed
elif [ "$1" == "line" ];then
  reg=LINE
  regbed=$dbdir/hg38_${reg}_onlyNanopore.bed
fi
outdir=$root/igv

# first expand the region
echo "bed"
bedbase=${regbed%.bed}
bed=$bedbase.10kb.bed
bedtools slop -b 5000 -i $regbed -g $gs |\
  sort -k1,1 -k2,2n | uniq > $bed

if [ "bsseq" == "$2" ];then
  echo "bsseq"
  bam=$root/data/bsseq/GM12878_BSseq.nodup.bam
  base=$(basename "$bam")
  base=${base%%.*}
  out=$outdir/$base.$reg.bam
  echo $out
  samtools view -@ 10 -L $bed -hb $bam > $out
  samtools index $out
  exit
fi

pooldir=$root/data/nanonome/pooled
bamdir=$pooldir/bam
mbeddir=$pooldir/mbed

for cell in $cells;do
  base=${cell}_nanoNOMe.pooled
  bam=$bamdir/$base.bam
  cpg=$mbeddir/$base.cpg.meth.bed.gz
  gpc=$mbeddir/$base.gpc.meth.bed.gz
  script=../../script/convertBam.py
  log=$outdir/$base.$reg.log
  out=$outdir/$base.$reg.bam

  com="python -u $script -v -t 10 \
    -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
    samtools sort -o $out"
  echo $com
  eval $com
  samtools index $out
done


