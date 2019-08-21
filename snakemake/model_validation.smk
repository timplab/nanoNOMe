#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

##################################################
#
# roc curve
# 
##################################################
rule make_roc_curve:
	input:
		unmeth="data/methylation_standards/"+
			"{sample}_unmethylated.subset.{model}.meth.tsv.gz",
		meth="data/methylation_standards/"+
			"{sample}_{methylation}.subset.{model}.meth.tsv.gz"
	params:
		config['codedir']
	output:
		log="plots/{sample}_ROC_{methylation}_methylation_{model}_model.log",
		plot="plots/{sample}_ROC_{methylation}_methylation_{model}_model.pdf"
	shell:
		"python -u {params}/scripts/test_methylmodel.py roc -v "
		"--model {wildcards.model} --methylated {input.meth} "
		"--unmethylated {input.unmeth} -o {output.plot} 2> {output.log}"

		
		
