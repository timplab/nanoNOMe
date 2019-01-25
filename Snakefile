configfile: 
  "snakemake_config.yml"
workdir:
  config['workdir']
smkdir = config['codedir']+"/snakemake/"
include: 
  smkdir+"downloaded_data_parse.smk"
include: 
  smkdir+"frequency_analysis.smk"
include:
  smkdir+"nanonome_data_parse.smk"

rule parse_downloaded_data:
  input:
    "data/gm12878/GM12878_CTCF.center.2000bp.bed"

rule frequency_analysis:
  input:
    "plots/GM12878.CTCF.metaplot.pdf"

rule nanonome_data:
  input:
    expand("data/nanonome/pooled/mfreq/{sample}_nanoNOMe.pooled.{mod}.mfreq.txt.gz",
      sample=["GM12878"],mod=["cpg","gpc"])
