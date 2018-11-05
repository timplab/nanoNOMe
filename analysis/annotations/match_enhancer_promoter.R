#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
root = commandArgs(TRUE)[1]
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis/annotations/breastcancer"
enh.fp = file.path(root,"MCF10A_enhancer_hg38.bed")
trans.fp = file.path(root,"../hg38/hg38_transcripts.bed")
trans.fp = file.path(root,"bcan_transcripts.bed")

# read files
bed.cnames = c("chrom","start","end","id","score","strand")
trans.cnames = c(bed.cnames,"geneid","name")
enh = read_tsv(enh.fp,col_names=bed.cnames)
trans = read_tsv(trans.fp,col_names=trans.cnames)

# match transcript id
trans.id = sapply(strsplit(trans$id,"[.]"),"[[",1)
dat.comb = bind_cols(enh,trans[match(enh$id,trans.id),]) %>% na.omit()

mat =  match(enh$id,trans.id)
trans.match = trans[mat,]%>%na.omit() 
enh.match = enh[which(!is.na(mat)),]
trans.match$matchid = enh.match$matchid = seq(dim(enh.match)[1])

# get promoter region and distance
prom.gr = promoters(GRanges(trans.match),1000)
TSS.gr = promoters(GRanges(trans.match),1,0)
comp.tb = bind_cols(enh.match,as.tibble(prom.gr))

reg.tb = trans.match
for (i in seq(dim(enh.match)[1])){
    if (enh.match$chrom[i] != trans.match$chrom[i]){
        reg.tb$chrom[i] = NA
        next
    }
    pos = c(enh.match$start[i],enh.match$end[i],
            trans.match$start[i],trans.match$end[i])
    reg.tb$start[i] = min(pos)
    reg.tb$end[i] = max(pos)
    reg.tb$strand = "."
}

maxw = 10000
reg.tb = reg.tb%>%
    mutate(width = end-start) %>%
    filter(width <= maxw)



dist.tb = tibble(start=enh.match$start-start(TSS.gr),
                 end=enh.match$end-start(TSS.gr))



distanceToNearest(enh.gr,prom.gr)
