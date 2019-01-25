#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

rule mbed_to_mfreq:
	input:
		"{dir}/mbed/{sample}.{mod}.meth.bed.gz"
	params:
		codedir=config['codedir'],
		log="{dir}/mfreq/{sample}.{mod}.mfreq.log"
	output:
		"{dir}/mfreq/{sample}.{mod}.mfreq.txt.gz"
	shell:
		"python -u {params.codedir}/script/parseMethylbed.py frequency -v "
		"-i {input} -m {wildcards.mod} 2> {params.log}| "
		"bgzip > {output} && "
		"tabix -b 2 -e 2 {output}"

