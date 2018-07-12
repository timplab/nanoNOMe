#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
source("../../plot/methylation_plot_utils.R")

# set this to TRUE to remove unnecessary objects throughout the process
limitedmem=TRUE

# set directories
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
plotdir=file.path(root,"plots/repeats")
if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
datroot=file.path(root,"pooled/methylation/mfreq_all")
regdir=file.path(root,"database/hg38")
# get file paths of data
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
fpaths=tibble(cell=cells,
          cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")),
          gpc=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)

# validation group (illumina)
if (F) {
    illpd=tibble(cell="GM12878_illumina",calltype=c("cpg","gpc"),
                 filepath=file.path("/dilithium/Data/Nanopore/projects/nomeseq/analysis/validation/scNOMe/methfreq",
                                    paste("GM12878_sample",calltype,"methfreq.txt.gz",sep=".")))
    pd=bind_rows(pd,illpd)
}

# regions
if (T){
    regnames=c("DNA","RTP","LINE","SINE")
    reg.info=tibble(regtype=regnames)%>%
        mutate(filepath=paste0(regdir,"/hg38_repeats_",regtype,".bed"))
    subset=FALSE
}
#read in regions
cat("reading in the region\n")
extracols="regtype"
#db=lapply(reg.info$filepath,function(x){
#    load_db(x,extracols)})
covthr=2
trinuc="GCG"

MethByRegion <- function(pd,reg){
    cat(paste0(reg,"\n"))
    dat.list=lapply(seq(dim(pd)[1]),function(i){
        pd.samp=pd[i,]
        cat(paste0(pd.samp$cell,":",pd.samp$calltype,"\n"))
        # read in the data
        dat=tabix_mfreq(pd.samp$filepath,dbpath=reg$filepath,
                        cov=covthr,trinuc_exclude=trinuc)
        # overlaps
        cat("getting methylation by region\n")
        db=load_db(reg,extracols)
        dat.ovl=getRegionMeth(dat,db)
        rm(dat);gc()
        #get labels
        cat("attaching labels\n")
        dat.ovl$feature.type=db$regtype[dat.ovl$feature.index]
        dat.ovl$samp=pd.samp$cell
        dat.ovl$calltype=pd.samp$calltype
        dat.ovl
    })
    dat.cat=do.call(rbind,dat.list)
})

dat.list=lapply(reg.info$filepath,function(x){
    MethByRegion(pd,x)})
dat.all=do.call(rbind,dat.list)

cpg.plt=dat.plt[which(dat.plt$calltype=="cpg"),]
gpc.plt=dat.plt[which(dat.plt$calltype=="gpc"),]

# boxplot
cat("plotting\n")
cpg.box=ggplot(cpg.plt,aes(x=feature.type,
                         y=freq,color=samp))+
    geom_boxplot(alpha=0.5)+
    theme_bw()
gpc.box=ggplot(gpc.plt,aes(x=feature.type,
                         y=freq,color=samp))+
    geom_boxplot(alpha=0.5)+
    theme_bw()

plotpath=file.path(plotdir,paste(regname,"boxplot_by_feature.pdf",sep="_"))
pdf(plotpath,width=15,height=5,useDingbats=F)
print(cpg.box)
print(gpc.box)
dev.off()


if (F) {
                                        # first start off iwth scatterplot
    dat.scatter = dat.filt%>% ungroup() %>%select(-cpgnum)%>%
        spread(sample,avgfreq)%>%select(-index.feature)
    g.duo = ggduo(dat.scatter,)
    g.pair = ggpairs(dat.scatter,)


    pdf(plotpath,width=10,height=10,useDingbats=F)
    print(g.box)
    print(g.duo)
    print(g.pair)
    dev.off()
}
