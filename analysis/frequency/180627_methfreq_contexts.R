#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
source("/home/isac/Code/ilee/plot/methylation_utils.R")

plotdir="~/Dropbox/Data/nome-seq/plots/frequency"
datroot="/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation/mfreq_all"
cell="GM12878"
cpgpath=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz"))
gpcpath=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz"))

# first the entire data
cpg.all=tabix_mfreq(cpgpath)

# filter out loci with <= cov coverage and GCG,CCG contexts
covthr=5
cpg = cpg.all%>%filter(cov>covthr &
                       trinuc != "GCG")
cpg$sample=cell
cpg.gr=GRanges(cpg)

# regions
if (T){
    regnames=c("tss2k","cgi","del400b")
    reg.fp=c("/mithril/Data/NGS/Reference/human_annotations/hg38.91.TSS.2kb.bed",
             "/mithril/Data/NGS/Reference/human_annotations/cpgIslandExtUnmasked.bed",
             "/dilithium/Data/Nanopore/projects/nomeseq/analysis/gm12878/ngmlr/sniffles/GM12878.sniffles.del.400b.region.bed")
    reg.info=tibble(region=regnames,filepath=reg.fp)
    subset=FALSE
}

#read in regions
db.list=lapply(seq(dim(reg.info)[1]),function(i){
    fp=reg.info$filepath[i]
    db=load_db(fp)
    db[which(width(db)<1000000)]
})
# overlap by the first region
db.list=lapply(db.list,function(x){
    x[overlapsAny(x,db.list[[1]])]})

# find overlaps
cpg.ovl=lapply(db.list,function(x){
    ovl=findOverlaps(cpg.gr,x)
    cpg[unique(queryHits(ovl)),]
})

# label each locus with feature index and distance to feature
dist.list = lapply(seq_along(db.list),function(i){
    dat=cpg.ovl[[i]]
    db=db.list[[i]]
    bind_cols(dat,getDistance(dat,db))})
# summarize by feature
sum.list=lapply(dist.list,function(dat){
    dat %>% group_by(index.feature,sample) %>%
        summarize(cpgnum=n(),avgfreq=mean(freq))
})
# get features that have > 10 cpg loci
sig.list=lapply(sum.list,function(dat){
    dat %>% ungroup() %>%
    group_by(index.feature) %>%
        summarize(cpgnum=min(cpgnum),
                  freqdiff=max(avgfreq)-min(avgfreq))%>%
        filter(cpgnum>10)
})
# subset if option is set
if (subset == TRUE) {
    sampid=sample(seq(dim(dat.num)[1]),1000)
    dat.num=dat.num[sampid,]
}

db.filt=lapply(seq_along(sig.list),function(i){
    sig.ind=sig.list[[i]]$index.feature
    db.list[[i]][sig.ind,]})

#get data in these features 
filt.list=lapply(seq_along(sig.list),function(i){
    dat=sum.list[[i]]
    dat[which(dat$index.feature %in% sig.list[[i]]$index.feature),]
})

# now label with feature name and cat the list
dat.cat = do.call(rbind,
                  lapply(seq_along(filt.list),function(i){
                      filt.list[[i]]$feature=reg.info$region[i]
                      filt.list[[i]]}))
# boxplot
g.box=ggplot(dat.cat,aes(x=feature,y=avgfreq,color=feature,fill=feature))+
    geom_boxplot(alpha=0.5)+
    theme_bw()
g.viol=g.box+geom_violin(alpha=0)
plotpath=file.path(plotdir,paste0(cell,"_boxplot_by_feature.pdf"))
pdf(plotpath,width=5,height=5,useDingbats=F)
print(g.box)
print(g.viol)
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
