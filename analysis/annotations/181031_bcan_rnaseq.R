#!/usr/bin/Rscript
library(tidyverse)
root = "/dilithium/Data/Nanopore/projects/nomeseq/analysis/annotations"
datfp = file.path(root,"breastcancer/bcan_rnaseq.txt")
genefp = file.path(root,"hg38/hg38_genes.bed")

db = read_tsv(genefp,col_names=c("chrom","start","end","id","score","strand","name"))
dat = read_tsv(datfp)
dat$p7 = dat$p231 = 1

# get t-stat
mindiff = 500
for (i in seq(dim(dat)[1])){
    print(i)
    m10 = as.numeric(dat[i,c(2:4)])
    m7 = as.numeric(dat[i,c(5:7)])
    mda = as.numeric(dat[i,c(8:10)])
    if (abs(mean(mda)-mean(m10)) > mindiff) {
        dat$p231[i] = try(t.test(log(m10+1),log(mda+1))$p.value)
    }
    if (abs(mean(m7)-mean(m10)) > mindiff) {
        dat$p7[i] = try(t.test(log(m10+1),log(m7+1))$p.value)
    }
}
# thresholding
a = 0.01
dat.filt = dat %>%
    filter(p7<=a|p231<=a) %>%
    arrange(p7)

# match names
# remove versions in ids
db.id = sapply(strsplit(db$id,"[.]"),"[[",1)
dat.id = sapply(strsplit(dat$Ensembl_ID,"[.]"),"[[",1)

db.231 = db[match(dat.id[sig.231],db.id),] %>% na.omit()
db.7 = db[match(dat.id[sig.7],db.id),] %>% na.omit()
