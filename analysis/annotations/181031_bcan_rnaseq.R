#!/usr/bin/Rscript
library(tidyverse)
root = "/dilithium/Data/Nanopore/projects/nomeseq/analysis/annotations"
datfp = file.path(root,"breastcancer/bcan_rnaseq.txt")
genefp = file.path(root,"hg38/hg38_genes.bed")

db = read_tsv(genefp,col_names=c("chrom","start","end","id","score","strand","name"))
dat = read_tsv(datfp)

# match names
# remove versions in ids
db.id = sapply(strsplit(db$id,"[.]"),"[[",1)
dat.id = sapply(strsplit(dat$Ensembl_ID,"[.]"),"[[",1)

dat = bind_cols(db[match(dat.id,db.id),],dat) %>% na.omit()
dat$p7 = dat$p231 = 1

# get t-stat
mindiff = 500
for (i in seq(dim(dat)[1])){
    print(i)
    m10 = as.numeric(dat[i,grep("MCF10A",names(dat))])
    m7 = as.numeric(dat[i,grep("MCF7",names(dat))])
    mda = as.numeric(dat[i,grep("MDA-MB-231",names(dat))])
    if (abs(mean(mda)-mean(m10)) > mindiff) {
        dat$p231[i] = try(t.test(log(m10+1),log(mda+1))$p.value)
    }
    if (abs(mean(m7)-mean(m10)) > mindiff) {
        dat$p7[i] = try(t.test(log(m10+1),log(m7+1))$p.value)
    }
}
# thresholding
a = 0.01
dat.7 = dat %>%
    filter(p7<=a) %>%
    arrange(p7)

dat.231 = dat %>%
    filter(p231<=a) %>%
    arrange(p231) 

names(dat.231)[1] = names(dat.7)[1] = "#chrom"
outdir=file.path(root,"breastcancer")
out.7 = file.path(outdir,"MCF10A_vs_MCF7_genes.bed")
write_tsv(dat.7,out.7)
out.231 = file.path(outdir,"MCF10A_vs_MDAMB231_genes.bed")
write_tsv(dat.231,out.231)
