#!/usr/bin/Rscript
library(tidyverse)
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(GenomicRanges)
library(GGally)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

root = commandArgs(trailingOnly=TRUE)[1]
# default path if not provided
if (is.na(root)){ root="/dilithium/Data/Nanopore/projects/nomeseq/analysis" }
cat(paste0("data root : ",root,"\n"))

# rnaseq data
exp.pre=file.path(root,"annotations/gm12878/GM12878_RNAseq_")
exp.fp=c(paste0(exp.pre,"1.tsv"),
         paste0(exp.pre,"2.tsv"))
exp.list=lapply(seq_along(exp.fp),function(i){
    read_tsv(exp.fp[i])%>%mutate(replicate=i)})
# gene data
gene.fp = file.path(root,"annotations/hg38/hg38_genes.bed")
genes = load_db(gene.fp,"genename")
# get quartiles
qtiles=seq(0,1,.25)
qtiles.list = lapply(exp.list,function(x){
    y = x %>% select(id=gene_id,
                     fpkm=FPKM_ci_upper_bound,
                     replicate=replicate)
    qs = quantile(y$fpkm,qtiles)
    y%>%mutate(qtile=ifelse(fpkm>=qs[2],2,1),
               qtile=ifelse(fpkm>=qs[3],3,qtile),
               qtile=ifelse(fpkm>=qs[4],4,qtile))
})
# using only genes that are consistent b/w the two replicates
qtiles = do.call(rbind,qtiles.list) %>%
    select(-fpkm) %>%
    spread(replicate,qtile) %>%
    filter(`1`==`2`)
# match with hg38 gene db
qtiles.id = sapply(strsplit(qtiles$id,"[.]"),"[[",1)
genes.id = sapply(strsplit(genes$id,"[.]"),"[[",1)
id.match = match(genes.id,qtiles.id)
genes$score = qtiles$`1`[id.match]
genes$score[which(is.na(genes$score))] = "."

# write output
bed = GRangesTobed(genes)
# in the final output, the score column has the quartile info
outpath = file.path(root,"annotations/gm12878/GM12878_gene_quartiles.bed")
cat(paste0("writing output to ",outpath,"\n"))
write_tsv(bed,outpath,col_names=F)
