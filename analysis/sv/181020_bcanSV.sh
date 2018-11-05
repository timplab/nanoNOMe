#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis
[ -z $1 ]||root="$1"
bamdir=$root/pooled/bam
mbeddir=$root/pooled/methylation/methbyread_all
svdir=$root/pooled/sv
parser=../../script/parse_sniffles.py
bammer=../../script/convertBam.py
cells="GM12878 MCF10A MCF7 MDAMB231"
db=$root/annotations/breastcancer/bcan_genes.bed
#cells="MCF10A"

# overlap SVs
comp=$svdir/SVcomparison.vcf
[ -e $comp ]||\
  python $srcdir/svoverlaps.py "$root"

# now make bed and find svs close to bcangenes
bcansv=$root/annotations/breastcancer/sv_bcangenes.bed
grep -E "oxx|oox|xoo|xxo" $comp |\
  python $srcdir/../../script/parse_sniffles.py bed -v -t 10 |\
  sort -k1,1 -k2,2n |\
  bedtools closest -D b -b $db -a stdin |\
  awk '{ if(sqrt($NF*$NF) <= 5000 && $11 != -1 && $3-$2 <= 100000 ) print }' \
  > $bcansv

exit




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
