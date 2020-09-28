#!/bin/bash
ref=/atium/Data/Reference/human/hg38_noalt/hg38_noalt.fa
dir=/uru/Data/Nanopore/projects/nanonome/haplotype
dars="/home/isac/Dropbox/Data/nome-seq/version3_guppy3/plots/haplotypes/200425_gm12878_haplotype_dars_sig.tsv"
vcf=/uru/Data/Nanopore/projects/nanonome/database/NA12878.vcf.gz
bed=$dir/gm12878_hap_dars.bed
slopbed=$dir/gm12878_hap_dars.1kbslop.bed
ovl=$dir/gm12878_hap_dars_with_het_snp.1kbslop.bed
atac=/dilithium/Data/NGS/projects/gm12878/atacseq/atacseq/bam/SRR891268.nodup.sorted.bam
sub=$dir/GM12878_ATAC_SRR891268_hap_dar_with_het_snp.bam
tagbam=$dir/GM12878_ATAC_SRR891268_hap_dar_with_het_snp.haplotag.bam
hap1bam=$dir/GM12878_ATAC_SRR891268_hap_dar_with_het_snp.hap1.bam
hap2bam=$dir/GM12878_ATAC_SRR891268_hap_dar_with_het_snp.hap2.bam
log=$dir/GM12878_ATAC_SRR891268_hap_dar_with_het_snp.haplotag.log
hap1cov=$dir/GM12878_ATAC_SRR891268_hap_dar_with_het_snp.hap1.cov.bed
hap2cov=$dir/GM12878_ATAC_SRR891268_hap_dar_with_het_snp.hap2.cov.bed

###################################################################################################################
#
# Step 1 : get DARs +- 1kb 
#
###################################################################################################################
if [ "$1" == "dars" ];then
  echo "getting dars"
  fai=$ref.fai
  awk -v OFS="\t" 'NR>1{ print $1,$2,$3 }' $dars |\
    sort -k1,1 -k2,2n -k3,3n > $bed
  bedtools slop -i $bed -b 1000 -g $fai |\
    sort -k1,1 -k2,2n -k3,3n > $slopbed
fi

###################################################################################################################
#
# Step 2 : Subset reads that are in the DARs +- 1kb 
#
###################################################################################################################
if [ "$1" == "sub" ];then
  echo "subsetting atac-seq data"
  samtools view -hb $atac -L $slopbed > $sub
  samtools index $sub
fi

###################################################################################################################
#
# Step 3 : Haplotype where possible
#
###################################################################################################################
if [ "$1" == "hap" ];then
  echo "haplotyping atac-seq data"
  whatshap haplotag --ignore-read-group -o $tagbam --reference $ref $vcf $sub 2> $log
  samtools index $tagbam
  samtools view -h $tagbam | awk '{if (/^[^@]/) { if (/HP:i:1/) {print}} else { print} }' |\
    samtools view -hb - > $hap1bam
  samtools index $hap1bam
  samtools view -h $tagbam | awk '{if (/^[^@]/) { if (/HP:i:2/) {print}} else { print} }' |\
    samtools view -hb - > $hap2bam
  samtools index $hap2bam
fi

###################################################################################################################
#
# Step 3 : get coverage inside the DARs
#
###################################################################################################################
if [ "$1" == "cov" ];then
  echo "getting haplotyped atac-seq coverae in dars"
  bedtools coverage -a $bed -b $hap1bam  | awk 'OFS="\t"{print $1,$2,$3,$(NF-3)}' > $hap1cov
  bedtools coverage -a $bed -b $hap2bam  | awk 'OFS="\t"{print $1,$2,$3,$(NF-3)}' > $hap2cov
fi

###################################################################################################################
#
# Step 4 : compare with nanoNOMe signal using R
#
###################################################################################################################
if [ "$1" == "comp" ];then
  Rscript 200720_gm12878_dars_atac_allelic_imbalance.R
fi

