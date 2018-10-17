#!/bin/bash
cells="MCF10A MCF7 MDAMB231"
cells="$cells GM12878"
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
bsseqdir=/dilithium/Data/NGS/projects/gm12878/bsseq/bamnodup
dbdir=$root/database/hg38
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
pooldir=$root/pooled
outdir=$pooldir/igv

# first expand the region
echo "bed"
bedbase=${regbed%.bed}
bed=$bedbase.10kb.bed
awk 'OFS="\t"{ if($2<5000){ $2=0 }else{ $2=$2-5000 };$3=$2+10000; if ($12>20 && $2 >=0 ) { print $1,$2,$3,$11 } }' $regbed |\
  sort -k1,1 -k2,2n | uniq > $bed
exit

if [ "bsseq" == "$2" ];then
  echo "bsseq"
  bam=$(find $bsseqdir -name "GM12878*nodup.bam")
  base=$(basename "$bam")
  base=${base%%.*}
  out=$outdir/$base.bsseq.$reg.bam
  echo $out
  samtools view -@ 10 -L $bed -hb $bam > $out
  samtools index $out
  exit
fi

for cell in $cells;do
  echo $cell
  bam=$pooldir/bam/$cell.pooled.bam
  mbeddir=$pooldir/methylation/methbyread_all
  cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
  script=../../script/convertBam.py
  log=$pooldir/igv/$cell.$reg.log
  outsort=$pooldir/igv/$cell.$reg.sorted.bam

  com="python -u $script -v -t 10 \
    -b $bam -c $cpg -g $gpc -r $bed 2> $log |\
    samtools sort -o $outsort"
  echo $com
  eval $com
  samtools index $outsort
done


