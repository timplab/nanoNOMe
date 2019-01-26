#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

rule mbed_to_mfreq:
	input:
		"{dir}/mbed/{sample}.{mod}.meth.bed.gz"
	params:
		config['codedir']
	output:
		mfreq="{dir}/mfreq/{sample}.{mod}.mfreq.txt.gz",
		log="{dir}/mfreq/{sample}.{mod}.mfreq.log"
	shell:
		"python -u {params}/script/parseMethylbed.py frequency -v "
		"-i {input} -m {wildcards.mod} 2> {output.log}| "
		"bgzip > {output.mfreq} && "
		"tabix -b 2 -e 2 {output.mfreq}"

