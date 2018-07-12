#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
source("/home/isac/Code/ilee/plot/methylation_utils.R")

limitedmem=TRUE
plotdir="~/Dropbox/Data/nome-seq/plots/frequency"
datroot="/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation/mfreq_all"
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
fpaths=tibble(cell=cells,
#          cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")))#,
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
    regname="repregions"
    reg.fp="/mithril/Data/NGS/Reference/human_annotations/nestedRepeats.bed"
    subset=FALSE
}
#read in regions
cat("reading in the region\n")
extracols="regtype"
db=load_db(reg.fp,extracols)

covthr=2
trinuc="GCG"
dat.list=lapply(seq(dim(pd)[1]),function(i){
    pd.samp=pd[i,]
    cat(paste0("sample : ",pd.samp$cell,"\n"))
    # read in the entire data first
    dat=tabix_mfreq(pd.samp$filepath,cov=covthr,trinuc_exclude=trinuc)
    # overlaps
    cat("finding overlaps\n")
    dat.ovl=getRegionMeth(dat,db)
    if (limitedmem) rm(dat);gc()
    #get labels
    cat("attaching labels\n")
    dat.ovl$feature.type=db$regtype[dat.ovl$feature.index]
    dat.ovl$samp=pd.samp$cell
    dat.ovl$calltype=pd.samp$calltype
    dat.ovl
})

dat.plt=do.call(rbind,dat.list)
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
