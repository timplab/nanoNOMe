#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

##################################################
#
# deletions
# 
##################################################
rule survivor_to_bed:
	input:
		"data/nanonome/pooled/sniffles/merged_cancer_1kbpdist.vcf.gz",
	params:
		config['codedir']
	output:
		"{dir}/{sample1}_vs_{sample2}_SVcomparison_{svtype}_survivor.bed"
	shell:
		"python {params}/scripts/parse_survivor.py {wildcards.svtype} "
		"-s {input} --one {wildcards.sample1} "
		"--two {wildcards.sample2} > {output}"

rule shuffle_bed:
	input:
		bed="{dir}/{sample}.bed",
		gs="data/hg38/hg38_genomesize.txt"
	output:
		"{dir}/{sample}.shuffle.bed"
	shell:
		"bedtools shuffle -i {input.bed} -g {input.gs} "
		"> {output}"

rule flank_bed:
	input:
		bed="{dir}/{sample}.bed",
		gs="data/hg38/hg38_genomesize.txt"
	params:
		config['codedir']
	output:
		"{dir}/{sample}.flank{width}bp.bed"
	shell:
		"bedtools flank -b {wildcards.width} -i {input.bed} " # flanking upstream and downstream
		"-g {input.gs} | awk 'OFS=\"\t\"{{ print $0,\".\",(NR+1)%2+1 }}' | " # add additional line for upstream vs downstream
		"sort -k1,1 -k2,2n > {output}"

rule get_TRA_region:
	input:
		bed="{dir}/{sample}.bed",
		gs="data/hg38/hg38_genomesize.txt"
	params:
		config['codedir']
	output:
		"{dir}/{sample}.TRAregion{width}bp.bed"
	shell:
		"bedtools slop -s -r 0 -l {wildcards.width} -i {input.bed} " # slopping upstream, strand-sensitive
		"-g {input.gs} | awk 'OFS=\"\t\"{{ print $0,(NR+1)%2+1 }}' | " # add additional line for upstream vs downstream
		"sort -k1,1 -k2,2n > {output}"

rule extract_mfreq_region:
	input:
		"data/bcan/{regname}.bed",
		"{dir}/pooled/mfreq/{sample}.mfreq.txt.gz"
	output:
		"{dir}/region/{sample}.region_{regname}.subset.bed"
	shell:
		"tabix -R {input} > {output}"
