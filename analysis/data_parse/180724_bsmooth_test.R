#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
library(bsseq)
source("../../plot/methylation_plot_utils.R")

# data
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/mfreq_all")
cell1="GM12878"
cell2="GM12878_wgs"
lowercell="gm12878"
pd=tibble(cell=c(cell1,cell1,cell2),
          type=c("cpg","gpc","cpg"),
          filepath=file.path(datroot,paste(cell,type,"methfreq.txt.gz",sep=".")))

plotdir=file.path(root,"plots/smooth")
dbname="tss"
dbname="ctcf"
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
}

if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
#read in tss db
db.gr=load_db(db.fp)

# read in data
dat.db=lapply(seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],db.fp)})

# bismark
bs.list=lapply(seq_along(dat.db),function(i){
    x=dat.db[[i]]
    BSseq(chr=x$chrom,
          pos=x$start,
          M=as.matrix(x$meth),
          Cov=as.matrix(x$cov),
          sampleNames=pd$cell[i])
})
# combine just the cpg
bs.cpg=combineList(bs.list[which(pd$type=="cpg")])
# gpc only in nome-seq sample
bs.gpc=bs.list[[which(pd$type=="gpc")]]

# smooth
cpg.cov = getCoverage(bs.cpg,type="Cov",what="perBase")
keepi = which(rowSums(cpg.cov>2)==dim(bs.cpg)[2])
params=expand.grid(ns=c(10,20,50),
                   h=c(100,500,1000))
smooth.cpg = lapply(seq(dim(params)[1]),function(i){
    ns = params$ns[i]
    h = params$h[i]
    bs=BSmooth(bs.cpg[keepi,],mc.cores=8,parallelBy="chromosome",verbose=TRUE,ns=ns,h=h)
    pData(bs) = data.frame(sample=sampleNames(bs),h=h,ns=ns)
    bs
})

# gpc
gpc.cov = getCoverage(bs.gpc,type="Cov",what="perBase")
keepi = which(rowSums(gpc.cov>2)==dim(bs.gpc)[2])
params=expand.grid(ns=c(5,10,20),
                   h=c(10,50,100))
smooth.gpc = lapply(seq(dim(params)[1]),function(i){
    ns = params$ns[i]
    h = params$h[i]
    bs=BSmooth(bs.gpc[keepi,],mc.cores=8,parallelBy="chromosome",verbose=TRUE,ns=ns,h=h)
    pData(bs) = data.frame(sample=sampleNames(bs),h=h,ns=ns)
    bs
})

# plotting

db.center=getCenter(db.gr)
# testing for various combos
if (F) {
    plotpath=file.path(plotdir,paste(cell1,dbname,"smooth_test.pdf",sep="_"))
    pdf(plotpath,width=6,height=4,useDingbats=F)
    for ( mod in c("cpg","gpc") ){
        if ( mod == "cpg" ) {
            bs.list = smooth.cpg
        } else if ( mod == "gpc" ) {
            bs.list = smooth.gpc
        }
        for ( pr in c("raw","smoothed") ){
            if ( pr == "raw" ){
                bs.in = bs.list[1]
            }else { bs.in = bs.list }
            plt.list = lapply(bs.in,function(bs){
                coords=as.tibble(as.data.frame(granges(bs)))
                if ( pr == "raw" ){
                    meth=getMeth(bs,type="raw",what="perBase")
                } else if ( pr == "smoothed" ){
                    meth=getMeth(bs,type="smooth",what="perBase")
                }
                ag.list = lapply(seq(dim(bs)[2]),function(i){
                    dat.tb = bind_cols(coords,as.data.frame(meth[,i]))
                    names(dat.tb) = c(names(coords),"freq")
                    dist.tb = bind_cols(dat.tb,getDistance(dat.tb,db.center))
                    aggregate_methylation(dist.tb,50)%>%
                        mutate(sample = pData(bs)$sample[i],
                               h = pData(bs)$h[i],
                               ns = pData(bs)$ns[i])
                })
                do.call(rbind,ag.list)})
            plt.tb = do.call(rbind,plt.list)
            plt.tb = plt.tb%>%
                mutate(lab = paste0("ns = ",ns," h = ",h))
            g = ggplot(plt.tb,aes(x=dist,y=freq,color=lab))+
                facet_grid(sample~.)+
                theme_bw()+geom_line(size=0.4,alpha=0.7)+
                lims(y=c(0,1)) + ggtitle(paste(mod,pr))
            print(g)
        }
    }
    dev.off()
}
# pick one
if (T) {
    plotpath=file.path(plotdir,paste(cell1,dbname,"smooth_compare.pdf",sep="_"))
    pdf(plotpath,width=6,height=4,useDingbats=F)
    for ( mod in c("cpg","gpc") ){
        if ( mod == "cpg" ) {
            params = do.call(rbind,lapply(smooth.cpg,function(x){pData(x)[1,]}))
            bs = smooth.cpg[[which(params$h == 500 & params$ns == 10)]]
        } else if ( mod == "gpc" ) {
            params = do.call(rbind,lapply(smooth.gpc,function(x){pData(x)[1,]}))
            bs = smooth.gpc[[which(params$h == 50 & params$ns == 10)]]
        }
        coords=as.tibble(as.data.frame(granges(bs)))
        meth.raw=getMeth(bs,type="raw",what="perBase")[,which(pData(bs)$sample==cell1)]
        meth.smooth=getMeth(bs,type="smooth",what="perBase")[,which(pData(bs)$sample==cell1)]
        meth = tibble(raw=meth.raw,smooth=meth.smooth)
        ag.list = lapply(seq(dim(meth)[2]),function(i){
            dat.tb = bind_cols(coords,as.data.frame(meth[,i]))
            names(dat.tb) = c(names(coords),"freq")
            dist.tb = bind_cols(dat.tb,getDistance(dat.tb,db.center))
            aggregate_methylation(dist.tb,50)%>%
                mutate(sample = names(meth)[i],
                       h = pData(bs)$h[1],
                       ns = pData(bs)$ns[1])
        })
        plt.tb = do.call(rbind,ag.list)
        plt.tb = plt.tb%>%
            mutate(lab = paste0("ns = ",ns," h = ",h))
        g = ggplot(plt.tb,aes(x=dist,y=freq,color=sample))+
            theme_bw()+geom_line(size=0.4,alpha=0.7)+
            lims(y=c(0,1)) + ggtitle(paste(mod,pr))
        print(g)
    }
    dev.off()
}
