#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R") 
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/qc")
cell="GM12878"
pd = tibble(cell=c(rep(cell,2),"GM12878_wgs"),type=c("cpg","gpc","cpg"),
            filepath=file.path(datroot,paste(cell,type,"pooled.meth.bed.gz",sep=".")))
if (T) {
    pd = tibble(cell=c("GM12878_wgs","GM12878_wgs","GM12878_nome","GM12878_nome","GM12878_nome"),
                type=c("cpg","gpc","cpg","cpg-new","gpc"),
                filepath=c("/dilithium/Data/Nanopore/projects/gm12878/analysis/mcall-cpg/GM12878.8.cpg.meth.tsv",
                           "/dilithium/Data/Nanopore/projects/gm12878/analysis/mcall-gpc/GM12878.9.gpc.meth.tsv.subset",
                           file.path(root,"gm12878/ngmlr/mcall-cpg/180322_GM12878_NOMe_rep3/180322_GM12878_NOMe_rep3.23.cpg.meth.tsv"),
                           file.path(root,"gm12878/ngmlr/mcall-cpg-new/180322_GM12878_NOMe_rep3/180322_GM12878_NOMe_rep3.23.cpg.meth.tsv"),
                           file.path(root,"gm12878/ngmlr/mcall-gpc/180322_GM12878_NOMe_rep3/180322_GM12878_NOMe_rep3.23.gpc.meth.tsv")))
    dat.list=lapply(pd$filepath,read_tsv)
    filt.list = lapply(seq_along(dat.list),function(i){
        x=dat.list[[i]]
        x[!grepl("GCG",x$sequence),]%>%
            transmute(calltype=pd$type[i],
                      sample=pd$cell[i],
                      score=log_lik_ratio)%>%
            sample_n(10000)
    })
    plt.tb = do.call(rbind,filt.list)
    plotpath = file.path(plotdir,"methylation_call_distribution_newvsold.pdf")
}

if (F){
    # chr22
    n = 10000000
    starts=seq(0,n*9,by=n)
    ends=starts+n
    reg.gr = GRanges(seqnames="chr22",ranges=IRanges(start=starts,end=ends))
    # load data in chr22
    dat.list = lapply(pd$filepath,function(x){
        tabix_mbed(x,reg.gr,by="read")})
    n = 1000
    dat.sub = lapply(dat.list,function(x){
        sample_n(x,size=n)
    })
    dat.call = lapply(seq_along(dat.sub),function(i){
        mbedByCall(dat.sub[[i]])%>%
            transmute(mcall,score,context,sample=pd$cell[i],calltype=pd$type[i])
    })

    plt.tb = do.call(rbind,dat.call)
    plt.tb = plt.tb[!grepl("GCG",plt.tb$context),]
    plotpath = file.path(plotdir,"methylation_call_distribution.pdf")
}

# distribution of ratios

pdf(plotpath,useDingbats=F)
for (i in seq(dim(pd)[1])){
    p = plt.tb[which(plt.tb$calltype==pd$type[i] & plt.tb$sample==pd$cell[i]),]
    print(length(which(abs(p$score)>=2.5))/length(p$score))
    
    g = ggplot(p,aes(x=sign(score)*log(abs(score))))+
        geom_histogram()+
        theme_bw()+
        geom_vline(xintercept=log(2.5))+
        geom_vline(xintercept=-log(2.5))+
        ggtitle(paste0(pd$cell[i]," : ",pd$type[i]))
    print(g)
}
dev.off()
