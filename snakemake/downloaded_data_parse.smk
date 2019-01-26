#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

tb = pd.read_csv(config['codedir']+"/data_download_info.csv")
fnames = list(tb['filename'])
urls = list(tb['url'])

##################################################
# general 
##################################################
rule download_data:
	params:
		url = 
			lambda wildcards:
				urls[fnames.index(wildcards.sample)]
	output:
		"download/{sample}.download"
	shell:
		"wget -q {params.url} -O {output}"

rule get_bed_regions:
	input:
		bed = "{dir}/{sample}.{region}.bed",
		gs = "data/hg38/hg38_genomesize.txt"
	output:
		"{dir}/{sample}.{region}.{width}bp.bed"
	shell:
		"export LC_ALL=C && "
		"side=$(({wildcards.width}/2)) && "
		"bedtools slop -b $side -i {input.bed} "
		"-g {input.gs} | sort -k1,1 -k2,2n | "
		"awk '$3-$2>{wildcards.width}{{print}}' > {output}"
	
##################################################
# hg38 reference and annotations
##################################################

rule gunzip_chain:
	input: "download/{genome}Tohg38.chain.gz.download"
	output: "download/{genome}Tohg38.chain"
	shell: "gunzip -c {input} > {output}"
	
rule parse_hg38_genome:
	input:
		"download/hg38.genome.fa.gz.download"
	params:
		codedir = config['codedir']
	output:
		fa = "data/hg38/hg38_noalt.fa",
		gs = "data/hg38/hg38_genomesize.txt"
	shell:
		"gunzip -c {input} | "
		"python {params.codedir}/util/fasta_remove_alt.py "
		"> {output.fa} && "
		"samtools faidx {output.fa} && "
		"cut -f1,2 {output.fa}.fai > {output.gs}"

rule parse_hg38_gtf:
	input:
		"download/{sample}.gtf.gz.download"
	output:
		bed="data/hg38/{sample}.genebody.bed",
		tss="data/hg38/{sample}.TSS.bed"
	shell:
		"gunzip -c {input} | "
		"awk 'OFS=\"\t\"{{ if($3==\"gene\") "
		"print $1,$4-1,$5,$10,\".\",$7,$14 }}' | "
		"tr -d \'\";\' | sort -k1,1 -k2,2n "
		"> {output.bed} && "
		"awk 'OFS=\"\t\"{{ if($6==\"+\") "
		"{{ $2=$2;$3=$2+1 }}else{{ $3=$3;$2=$3-1 }} print }}' "
		"{output.bed} | sort -k1,1 -k2,2n > {output.tss} "

##################################################
# gm12878
##################################################

rule parse_ctcfbsdb:
	input:
		db="download/ctcfbsdb.txt.gz.download",
		chain="download/hg18Tohg38.chain"
	params:
		oldbed=temp("download/CTCF_bindingsites_hg18.bed"),
		unsorted=temp("data/hg38/CTCF_bindingsites_hg38.unsorted.bed")
	output:
		"data/hg38/CTCF_bindingsites_hg38.bed"
	shell:
		"gunzip -c {input.db} | "
		"awk '$2==\"Human\"{{ print $3 }}' | "
		"tr \":\" \"\t\" | tr \"-\" \"\t\" | "
		"awk 'OFS=\"\t\"{{ print $1,$2-1,$3 }}' "
		"> {params.oldbed} && "
		"liftOver {params.oldbed} {input.chain} "
		"{params.unsorted} {output}.unmapped && "
		"sort -k1,1 -k2,2n {params.unsorted} > {output}"

rule parse_chip:
	input:
		chip="download/{sample}_{context}_ChIP.bed.gz.download",
		ref="data/hg38/{context}_bindingsites_hg38.bed",
		tss="data/hg38/hg38_genes.TSS.bed"
	params: config['codedir']
	output:
		sites="{dir}/{sample}_{context}.bindingsites.bed",
		center="{dir}/{sample}_{context}.center.bed"
	shell:
		"gunzip -c {input.chip} | "
		"bedtools intersect -a {input.ref} -b stdin -u | " #intersect with reference
		"sort -k1,1 -k2,2n | "
		"bedtools closest -a stdin -b {input.tss} -d | "
		"awk 'OFS=\"\t\"{{ if($NF>2000)print $1,$2,$3 }}' " # remove features closer than 2kb
		"> {output.sites} && "
		"python {params}/util/bed_getcenter.py "
		"{output.sites} | sort -k1,1 -k2,2n > {output.center}"




