#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

rule mbed_to_mfreq:
	input:
		"{dir}/mbed/{sample}.{mod}.meth.bed.gz"
	params:
		config['codedir']
	output:
		"{dir}/mfreq/{sample}.{mod}.mfreq.txt.gz"
	shell:
		"python {params}/script/parseMethylbed.py frequency -v "
		"-i {input} -m {wildcards.mod} | "
		"bgzip > {output} && "
		"tabix -b 2 -e 2 {output}"
		
	


