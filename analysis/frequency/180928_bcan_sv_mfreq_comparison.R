#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(GenomicRanges)
library(parallel)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

# set this to TRUE to remove unnecessary objects throughout the process
limitedmem=TRUE
cores = detectCores()-2

# set directories
root=commandArgs(trailingOnly=TRUE)[1]
if (is.na(root)){
    root="/dilithium/Data/Nanopore/projects/nomeseq/analysis" # default
}
outdir=file.path(root,"annotations/breastcancer")
outpath=file.path(outdir,"bcan_10a_vs_231_promoters.bed")
datroot=file.path(root,"pooled/methylation/mfreq_all")
cells=c("MCF10A","MCF7","MDAMB231")
fpaths=tibble(cell=cells,
          cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")),
          gpc=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)

dbpath = file.path(root,"annotations/breastcancer/bcan_genes.bed")
db.gr = load_db(dbpath,extracols=c("genename","fxn"))
prom.gr = promoters(db.gr,200,200)

dat.list = mclapply(mc.cores=cores,seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],prom.gr)
})

regmeth.list = lapply(seq_along(dat.list),function(i){
    getRegionMeth(dat.list[[i]],db.gr) %>%
        mutate(cell=pd$cell[i],
               calltype=pd$calltype[i],
               calltype=ifelse(calltype=="cpg","CpG Methylation","GpC Accessibility"))
})
regmeth = do.call(rbind,regmeth.list)

meth.spread = regmeth %>% select(-totcov,-numsites) %>%
    spread(cell,freq) %>%
    replace(.,is.na(.),0)%>%
    mutate(del=MDAMB231-MCF10A) %>%
    arrange(desc(abs(del)))

topregs = prom.gr[unique(head(meth.spread$feature.index,n=20))]
bed.tb = as.tibble(topregs)
write_tsv(bed.tb,outpath,col_names=F)
