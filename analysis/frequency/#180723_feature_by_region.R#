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
plotdir=file.path(root,"plots/regions")
dbname="tss"
#dbname="ctcf"
dbname="tss_shuffle"
dbname="genebody"
if ( dbname=="tss"){
    db.fp=file.path(root,"database/hg38/hg38_genes.TSS.2000bp.bed")
    names.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names=read_tsv(names.fp,col_names=c("chrom","start","end","id","score","strand","name","fxn"))
    subset=FALSE
    # plot dir
    plotpath=file.path(plotdir,"GM12878_TSS_regions.pdf")
}else if (dbname=="ctcf"){
    db.fp=file.path(root,"database/gm12878/ctcf/GM12878_ctcf.2000bp.bed")
    plotpath=file.path(plotdir,"GM12878_CTCF_regions.pdf")
}else if (dbname == "tss_shuffle"){
    db.fp=file.path(root,"database/hg38/hg38_genes.TSS.2000bp.shuffle.bed")
    names.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names=read_tsv(names.fp,col_names=c("chrom","start","end","id","score","strand","name","fxn"))
}else if (dbname == "genebody") {
    db.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names=read_tsv(names.fp,col_names=c("chrom","start","end","id","score","strand","name","fxn"))
}

if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
#read in tss db
db.gr=load_db(db.fp)
db.gr$name=names[match(db.gr$id,names$id),]$name

# read in data around feature
dat.db=lapply(seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],db.fp)})

# get methylation by region
regmeth.list=lapply(dat.db,function(x){
    getRegionMeth(x,db.gr)})

# order based on CpG
cpgind=which(pd$type=="cpg")
cpg.regmeth=regmeth.list[[cpgind]] %>%
    mutate(cov=totcov/numsites) %>%
    filter(numsites>20,
           cov<=50,cov>10) %>%
    arrange(desc(numsites),desc(cov))
#gpc.regmeth=regmeth.list[[-cpgind]]%>%
#    mutate(cov=totcov/numsites)
#gpc.regmeth=gpc.regmeth[which(gpc.regmeth$feature.index %in% cpg.regmeth$feature.index),]
db.sig=db.gr[cpg.regmeth$feature.index]
sig.sub = db.sig[1:1000,]

# do TSS/binding region scatterplots
if ( TRUE ){
    if (dbname == "tss" | dbname == "tss_shuffle"){
        sub.tb = as.tibble(as.data.frame(sig.sub))
        window=1000
        tss.tb = transmute(sub.tb,
                           seqnames,strand,name,
                           start=ifelse(strand=="-",start+window,start+window+1),
                           end=start)
        snum=5
        sections=seq(0,window,length.out=snum)
        labels=c(1,2,3,4,5,6)#c("one","two","three","four","five","six")
        params = tibble(upstream=sections[-snum],
                        downstream=sections[-1])
        params$lab = labels[1:dim(params)[1]]
        tss.str = mutate(tss.tb,start = ifelse(strand=="-",-start,start))
        regions.list = lapply(seq(dim(params)[1]),function(i){
            up=tss.str$start+params$upstream[i]
            down=tss.str$start+params$downstream[i]
            y = tss.str %>% mutate(start=ifelse(strand=="-",-down,up),
                                   end=ifelse(strand=="-",-up,down),
                                   lab=params$lab[i])
            GRanges(y)
        })
        regions=do.call("c",regions.list)
        regmeth.list = lapply(seq_along(dat.db),function(i){
            x = dat.db[[i]]
            getRegionMeth(x,regions)%>%select(-totcov,-numsites)%>%
                mutate(calltype=pd$type[i])
        })
        regmeth = do.call(rbind,regmeth.list) %>%
            spread(calltype,freq)
        regmeth$lab = regions$lab[regmeth$feature.index]
        g = ggplot(regmeth,aes(x=gpc,y=cpg))+
            facet_grid(lab~.)+
            geom_point(size=0.5,alpha=0.5)+
            theme_bw()
        plotpath = file.path(plotdir,paste0("GM12878_",dbname,"_scatter_sections.pdf"))
        pdf(plotpath,width=4,height=10,useDingbats=F)    
        print(g)
        dev.off()
        # entire 1kb region
        tss.gr = GRanges(tss.tb)
        regions = promoters(tss.gr,upstream=0,downstream=1000)
        regmeth.list = lapply(seq_along(dat.db),function(i){
            x = dat.db[[i]]
            getRegionMeth(x,regions)%>%select(-totcov,-numsites)%>%
                mutate(calltype=pd$type[i])
        })
        regmeth = do.call(rbind,regmeth.list) %>%
            spread(calltype,freq)
        g = ggplot(regmeth,aes(x=gpc,y=cpg))+
            geom_point(size=0.5,alpha=0.5)+
            theme_bw()
        plotpath = file.path(plotdir,paste0("GM12878_",dbname,"_scatter.pdf"))
        pdf(plotpath,width=4,height=3,useDingbats=F)    
        print(g)
        dev.off()
    }
    if (dbname == "ctcf") {
        regions = resize(db.gr,width = 400, fix = "center")
        regmeth.list = lapply(seq_along(dat.db),function(i){
            x = dat.db[[i]]
            getRegionMeth(x,regions)%>%select(-totcov,-numsites)%>%
                mutate(calltype=pd$type[i])
        })
        regmeth = do.call(rbind,regmeth.list) %>%
            spread(calltype,freq)
        g = ggplot(regmeth,aes(x=gpc,y=cpg))+
            geom_point(size=0.5,alpha=0.5)+
            theme_bw()
        plotpath = file.path(plotdir,"GM12878_ctcf_scatter.pdf")
        pdf(plotpath,width=4,height=3,useDingbats=F)    
        print(g)
        dev.off()
    }
    if (dbname == "genebody") {
        # entire region
        regmeth.list = lapply(seq_along(dat.db),function(i){
            x = dat.db[[i]]
            getRegionMeth(x,db.gr)%>%select(-totcov,-numsites)%>%
                mutate(calltype=pd$type[i])
        })
        regmeth = do.call(rbind,regmeth.list) %>%
            spread(calltype,freq)
        g = ggplot(regmeth,aes(x=gpc,y=cpg))+
            geom_point(size=0.5,alpha=0.5)+
            theme_bw()
        plotpath = file.path(plotdir,paste0("GM12878_",dbname,"_scatter.pdf"))
        pdf(plotpath,width=4,height=3,useDingbats=F)    
        print(g)
        dev.off()
    }
}

# write out regions as tsv
if ( FALSE ){
    sig.bed=as.tibble(as.data.frame(db.sig))%>%
        select(-width,-regstart,-regend)%>%
        mutate(start=start-1)
    bedpath=file.path(plotdir,paste0("GM12878_",dbname,"_regions.bed"))
    write_tsv(sig.bed,bedpath,col_names=FALSE)
}

# plotting regions
if ( FALSE ){
    # re-do the data extraction based on new subset of genes
    ovl.list=lapply(dat.db,function(x){
        findOverlaps(GRanges(x),db.sig)})
    pdf(plotpath,width=5,height=5,useDingbats=F)    
    for (i in seq_along(sig.sub)){
        reg=db.sig[i]
        sub.list=lapply(seq_along(ovl.list),function(j){
            x=ovl.list[[j]]
            dat.db[[j]][queryHits(x)[which(subjectHits(x)==i)],]%>%
                mutate(calltype=pd$type[j])
        })
        dat.sub=do.call(rbind,sub.list)%>%
            select(start,cov,freq,calltype)%>%
            gather(what,value,-start,-calltype)
        g=ggplot(dat.sub,aes(x=start,y=value,color=calltype,group=calltype))+
            facet_grid(what~.,scales="free_y")+
            geom_smooth(se=F,method="loess",span=0.2)+
            geom_point(size=0.5,alpha=0.5)+
            labs(x=paste0("Distance along ",seqnames(reg)),
                 y="Methylation Frequency")+
            theme_bw()
        if(dbname=="tss"){
            if(as.character(strand(reg))=="-"){
                gene=tibble(start=start(reg),
                            end=reg$regstart)
            }else{
                gene=tibble(start=reg$regend,
                            end=end(reg))
            }
            genename=names[match(reg$id,names$id),]$name
            gene$Gene=genename
            gene$what="freq"
            g=g+
                geom_rect(inherit.aes=F,data=gene,
                          mapping=aes(xmin=start,xmax=end,
                                      ymin=0,ymax=1,fill=Gene),alpha=0.3)
        }else if (dbname=="ctcf"){
            center=tibble(center=reg$regstart)
            g=g+geom_vline(xintercept=center$center)
        }
        g=g+theme(legend.position="bottom")
        print(g)
    }
    dev.off()

}
