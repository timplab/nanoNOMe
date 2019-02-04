#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

##################################################
#
# metaplots
# 
##################################################
def get_dir_from_samplename(sample):
	if "nanoNOMe" in sample :
		prefix = "data/nanonome/pooled/"
	elif "BSseq" in sample :
		prefix = "data/bsseq/"
	return prefix

rule make_gm12878_metaplot:
	input:
		cpg= 
			lambda wildcards :
				get_dir_from_samplename(wildcards.sample)+
					"mfreq/{sample}.cpg.mfreq.txt.gz",
		gpc= 
			lambda wildcards :
				get_dir_from_samplename(wildcards.sample)+
					"mfreq/{sample}.gpc.mfreq.txt.gz",
		region="data/gm12878/GM12878_{region}.bed"
	params:
		config['codedir']
	output:
		"plots/{sample}.metaplot.{region}.pdf"
	shell :
		"Rscript {params}/script/nanonome_plots.R "
		" metaplotByDistance -c {input.cpg} -g {input.gpc} "
		"-r {input.region} -o {output} 2> /dev/null"

rule overlap_bed_regions:
	input:
		data="data/{dir}/{sample}.bed.gz",
		region="data/gm12878/GM12878_{region}.bed",
	output:
		"data/{dir}/{sample}.distance.{region}.{width}bp.bed"
	shell:
		"gunzip -c {input.data} | "
		"bedtools closest -D b -b {input.region} -a stdin | "
		"awk 'sqrt($NF*$NF)*2<={wildcards.width}{{ print }}' "
		"> {output}"

rule make_mnase_metaplot:
	input:
		"data/mnase/{sample}_MNase.distance.{region}.bed"
	params:
		config['codedir']
	output:
		"plots/{sample}_MNase.metaplot.{region}.pdf"
	shell:
		"Rscript {params}/script/mnase_plot.R metaplotByDistance "
		"-i {input} -o {output} 2> /dev/null"

# need to work on
rule plot_methylation_chromatin_by_expression_quartile:
	input:
		"data/nanonome/pooled/mfreq/{sample}.cpg.mfreq.txt.gz",
		"data/nanonome/pooled/mfreq/{sample}.gpc.mfreq.txt.gz"
	params:
		config['codedir']
	shell:
		"Rscrpt {params}/analysis/frequency/methylation_chromatin_by_expression.R "
		"blah blah"
		
