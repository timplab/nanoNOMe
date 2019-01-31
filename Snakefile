configfile: 
  "snakemake_config.yml"
workdir:
  config['workdir']
smkdir = config['codedir']+"/snakemake/"
include: 
  smkdir+"downloaded_data_parse.smk"
include:
  smkdir+"model_validation.smk"
include: 
  smkdir+"frequency_analysis.smk"
include:
  smkdir+"nanonome_data_parse.smk"

rule parse_downloaded_data:
  input:
    "data/gm12878/GM12878_CTCF.noTSS.center.2000bp.bed",
    "data/hg38/hg38_genes.TSS.200bp.bed",
    "data/hg38/hg38_cgi.bed"

rule model_validation:
  input:
    "plots/NA12878_ROC_CpG_methylation_cpg_model.pdf",
    "plots/NA12878_ROC_GpC_methylation_gpc_model.pdf"

rule frequency_analysis:
  input:
    "plots/GM12878_nanoNOMe.pooled.metaplot.CTCF.noTSS.center.2000bp.pdf",
    "plots/GM12878_BSseq.metaplot.CTCF.noTSS.center.2000bp.pdf",
    "plots/GM12878_MNase.metaplot.CTCF.noTSS.center.2000bp.pdf"

rule nanonome_data:
  input:
    expand("data/nanonome/pooled/bigwig/{sample}_nanoNOMe.pooled.{mod}.{what}.bw",
      sample=["GM12878","MCF10A","MCF7","MDAMB231"],
        mod=["cpg","gpc"],
        what=["methylation","methcoverage"])
