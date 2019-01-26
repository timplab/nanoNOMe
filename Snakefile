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
    "plots/GM12878_nanoNOMe.pooled.metaplot.CTCF.center.2000bp.pdf",
    "plots/GM12878_BSseq.metaplot.CTCF.center.2000bp.pdf",
    "plots/GM12878_MNase.metaplot.CTCF.center.2000bp.pdf"

rule nanonome_data:
  input:
    expand("data/nanonome/pooled/mfreq/{sample}_nanoNOMe.pooled.{mod}.mfreq.txt.gz",
      sample=["GM12878","MCF10A","MCF7","MDAMB231"],mod=["cpg","gpc"])
