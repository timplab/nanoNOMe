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
plotpath=file.path(root,"plots/promoters")
datroot=file.path(root,"pooled/methylation/mfreq_all")
cells=c("MCF10A","MDAMB231")
fpaths=tibble(cell=cells,
          cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")),
          gpc=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)
cgi = file.path(root,"annotations/hg38/hg38_cgi.bed")

dbpath = file.path(root,"annotations/breastcancer/bcan_genes.TSS.400bp.bed")

dat.list = mclapply(mc.cores=cores,seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],dbpath)
})
db.gr = load_db(dbpath,"genename")
regmeth.list = lapply(seq_along(dat.list),function(i){
    getRegionMeth(dat.list[[i]],db.gr) %>%
        mutate(cell=pd$cell[i],
               calltype=pd$calltype[i])
})
regmeth = do.call(rbind,regmeth.list)

meth.spread = regmeth %>% select(-totcov,-numsites) %>%
    spread(cell,freq) %>%
    replace(.,is.na(.),0)%>%
    mutate(del=MDAMB231-MCF10A) %>%
    select(-MCF10A,-MDAMB231) %>%
    spread(calltype,del) %>%
    na.omit() %>%
    arrange(desc(abs(gpc)),desc(abs(cpg)))

topregs = db.gr[unique(head(meth.spread$feature.index,n=20))]
awidth = 5000
regs.adjust = resize(topregs,width=awidth,fix="center")
bed.tb = GRangesTobed(regs.adjust)
write_tsv(bed.tb,outpath,col_names=F)
