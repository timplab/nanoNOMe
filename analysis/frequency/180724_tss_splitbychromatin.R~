#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
source("../../plot/methylation_plot_utils.R")

# data
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/mfreq_all")
cell="GM12878"
lowercell="gm12878"
cpgpath=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz"))
gpcpath=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz"))
pd=tibble(cell=cell,type=c("cpg","gpc"),filepath=c(cpgpath,gpcpath))
plotwhat="correlation"
if ( plotwhat == "aggregateByexpression" ){
    tss.fp=file.path(root,"database/hg38/hg38_genes.TSS.2000bp.bed")
    subset=FALSE
    # plot dir
    plotdir=file.path(root,"plots/aggregate")
    plotpath=file.path(plotdir,paste0(cell,"_TSSbyExp.pdf"))
} else if ( plotwhat == "correlation" ){
    tss.fp=file.path(root,"database/hg38/hg38_genes.TSS.400bp.bed")    
    subset=TRUE
    # plot dir
    plotdir=file.path(root,"plots/correlation")
    plotpath=file.path(plotdir,paste0(cell,"_CpGvsGpCvsExp_scatter.pdf"))
}
if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
#read in tss db
tss.gr=load_db(tss.fp)

# rnaseq database
exp.dir=file.path(root,"database",lowercell,"rnaseq")
exp.fp=c(file.path(exp.dir,paste0(cell,"_rnaseq.1.tsv")),
         file.path(exp.dir,paste0(cell,"_rnaseq.2.tsv")))
# expression data parsing
qtiles=seq(0,1,.25)
exp.list=lapply(exp.fp,read_tsv)
exp.list=lapply(seq_along(exp.list),function(i){
    x=exp.list[[i]]
    id=sapply(strsplit(x$gene_id,"[.]"),"[[",1)
    y=tibble(id=id,fpkm=x$FPKM_ci_upper_bound)
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
# using only genes that are consistent b/w the two replicates, taking average
exp.fpkm=exp[which(exp$id %in% exp.qtile$id),] %>%
    group_by(id)%>%summarize(fpkm=mean(fpkm))
# merge with tss.gr
tss.gr=tss.gr[which(tss.gr$id %in% exp.fpkm$id)]
tss.gr$qtile=exp.qtile$qtile[match(tss.gr$id,exp.qtile$id)]
tss.gr$fpkm=exp.fpkm$fpkm[match(tss.gr$id,exp.fpkm$id)]

# read in data around TSS
dat.tss=lapply(seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],tss.fp)})

if (plotwhat == "correlation"){
    regmeth.list=lapply(dat.tss,function(x){
        getRegionMeth(x,tss.gr)})
    regmeth=do.call(rbind,lapply(seq_along(regmeth.list),function(i){
        regmeth.list[[i]]%>%
            mutate(calltype=pd$type[i],
                   cell=pd$cell[i])}))
    dat.spread = regmeth%>%
        filter(numsites>10)%>%
        select(-totcov,-numsites)%>%
        spread(calltype,freq)%>%na.omit()
    dat.spread$fpkm=tss.gr$fpkm[dat.spread$feature.index]
    # quantile normalization
    qtiles=seq(0,1,0.02)
    breaks=quantile(dat.spread$fpkm,probs=qtiles)
    dat.spread$exp=as.numeric(as.character(
        cut(dat.spread$fpkm,
            breaks,
            include.lowest=T,labels=qtiles[-1])))
    # just log - offset by min
    minexp=min(dat.spread$fpkm[which(dat.spread$fpkm!=0)])
    dat.spread$exp=log(dat.spread$fpkm+minexp) 
    dat.plt=dat.spread %>%
        select(c("cpg","gpc","exp"))
    subn=1000
    # subset
    sub=sample(seq(dim(dat.plt)[1]),subn)
    dat.sub=dat.plt[sub,]
    # plot
    g.pairs=ggpairs(dat.plt,aes(alpha=0.5))
    g.duo=ggduo(dat.sub)
    pdf(plotpath,width=6,height=4,useDingbats=F)
    print(g.pairs)
    print(g.duo)
    dev.off()

}

if ( plotwhat == "aggregateByexpression" ){
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

}

