#!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
bamdir=$root/pooled/bam
mbeddir=$root/pooled/methylation/methbyread_all
svdir=$root/pooled/sv
parser=../../script/parse_sniffles.py
bammer=../../script/convertBam.py
cells="GM12878 MCF10A MCF7 MDAMB231"
#cells="MCF10A"

for cell in $cells;do
  echo $cell
  bam=$bamdir/$cell.pooled.bam
  sv=$svdir/$cell.sniffles.vcf
  cpg=$mbeddir/$cell.cpg.pooled.meth.bed.gz
  gpc=$mbeddir/$cell.gpc.pooled.meth.bed.gz
  out=$svdir/$cell.bcan_exclusive_TRA
  outbam=$svdir/$cell.TRAmultiregion.bam
  win=250
  svbed=$svdir/$cell.svregion.$win.bed
  svbed=$svdir/MDAMB231.svregion.250.bed

  if [ "$1" == "makebed" ];then
    log=$svdir/$cell.svbed.$win.log
    python -u $parser bed -v -t 8 -b $bam -s $sv -o $svbed -w $win 2> $log
  fi
  if [ "$1" == "bam" ];then
    bamout=$svdir/$cell.MDAMB231SVregion250bp.bam
    log=$bamout.log
    python -u $bammer -v -t 10 -b $bam -c $cpg -g $gpc -r $svbed 2> $log |\
      samtools sort -o $bamout
    samtools index $bamout
  fi
    

done

if [ "$1" == "svcompare" ];then
  python ./svoverlaps.py
fi

if [ "$1" == "meth" ];then
  Rscript ./180910_SVmethylation_comparison.R
fi
