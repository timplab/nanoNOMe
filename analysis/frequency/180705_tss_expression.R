#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
source("/home/isac/Code/ilee/plot/methylation_utils.R")

# data
datroot="/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation/mfreq_all"
cell="GM12878"
cpgpath=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz"))
gpcpath=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz"))
pd=tibble(cell=cell,type=c("cpg","gpc"),filepath=c(cpgpath,gpcpath))
# plot dir
plotdir="~/Dropbox/Data/nome-seq/plots/aggregate"
plotpath=file.path(plotdir,paste0(cell,".tssByExp.pdf"))
# databases
if (T){
    tss.fp="/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.2kb.bed"
    subset=FALSE
    exp.dir="/home/isac/Dropbox (Timp Lab)/Data/nome-seq/db/gm12878/rnaseq"
    exp.fp=c(file.path(exp.dir,paste0(cell,"_rnaseq.1.tsv")),
             file.path(exp.dir,paste0(cell,"_rnaseq.2.tsv")))
}

# read in data around TSS
dat.tss=lapply(seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],tss.fp)})
#read in tss db
tss.gr=load_db(tss.fp)
qtiles=seq(0,1,.25)
# expression data parsing
exp.list=lapply(exp.fp,read_tsv)
exp.list=lapply(seq_along(exp.list),function(i){
    x=exp.list[[i]]
    id=sapply(strsplit(x$gene_id,"[.]"),"[[",1)
    y=tibble(id=id,fpkm=x$pme_FPKM)
    y=y[which(y$id %in% tss.gr$id),]
    qs=quantile(y$fpkm,probs=c(0.25,0.5,0.75))
    y%>%mutate(qtile=ifelse(fpkm>=qs[1],2,1),
               qtile=ifelse(fpkm>=qs[2],3,qtile),
               qtile=ifelse(fpkm>=qs[3],4,qtile),
               rep=i)
})
exp=do.call(rbind,exp.list)
exp.qtile=exp%>%select(-fpkm)%>%
    spread(rep,qtile)%>%filter(`1`==`2`)%>%mutate(qtile=`1`)
# merge with tss.gr
tss.gr=tss.gr[which(tss.gr$id %in% exp.qtile$id)]
tss.gr$qtile=exp.qtile$qtile[match(tss.gr$id,exp.qtile$id)]

# re-do the data extraction based on new subset of genes
dat.exp=lapply(dat.tss,function(x){
    x[overlapsAny(GRanges(x),tss.gr),]})
# get distances
tss.center=getCenter(tss.gr)
dist.list=lapply(dat.exp,function(x){
    bind_cols(x,getDistance(x,tss.center))})
# add quartile label
dist.list=lapply(dist.list,function(x){
    x$qtile=tss.gr$qtile[x$index.feature];x})

# aggregate by quantile
agg.list=lapply(seq(length(qtiles)-1),function(i){
    dist.q=lapply(dist.list,function(x){
        x[which(x$qtile==i),]})
    agg.q=lapply(dist.q,aggregate_methylation)
    # combine into one
    agg.q=lapply(seq_along(agg.q),function(i){
        x=agg.q[[i]]
        x%>%mutate(cell=pd$cell[i],
                   mcall=pd$type[i],
                   lab=paste(cell,mcall,sep="."))})
    agg.tb=do.call(rbind,agg.q)
    agg.tb$qtile=i;agg.tb
})
dat.plt=do.call(rbind,agg.list)
dat.plt$qtile=factor(dat.plt$qtile)

# plot
g=ggplot(dat.plt,aes(x=dist,y=freq,linetype=qtile,color=mcall))+
    theme_bw()+geom_line()+
    lims(y=c(0,1))+
    labs(x="Distance to TSS", y="Average methylation")

pdf(plotpath,width=6,height=4,useDingbats=F)
print(g)
dev.off()
