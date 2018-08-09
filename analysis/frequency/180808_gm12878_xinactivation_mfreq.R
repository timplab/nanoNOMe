#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R") 
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
plotdir=file.path(root,"plots/readlevel")
cell="GM12878"
what = "average"

# x chrom vs chrom 21,22
start = 0
end = 1e9
chrom.gr = GRanges(c("chr21","chr22","chrX"),
                   IRanges(start=c(rep(start,2),67676484),
                               end=c(rep(end,2),76800000)))
chrom.gr = GRanges(c("chr21","chr22","chrX"),
                   IRanges(start=c(rep(start,3)),
                               end=c(rep(end,3))))

if (what == "average") {
    plotpre = file.path(root,"plots/readlevel/180808_GM12878_XvsAuto_mfreq")
    datroot=file.path(root,"pooled/methylation/mfreq_all")
    pd = tibble(cell=rep(c(cell,"GM12878_wgs","GM12878_BSseq_ENCLB898WPW","MCF10A"),each=2),type=rep(c("cpg","gpc"),4),
                filepath=file.path(datroot,paste(cell,type,"methfreq.txt.gz",sep=".")))
    dat.list = mclapply(mc.cores=4,pd$filepath,function(x){
        tabix_mfreq(x,chrom.gr)
    })
    call.reads = lapply(seq_along(dat.list),function(i){
        x = dat.list[[i]]
        if("tbl" %in% class(x)){
            x%>%mutate(mcall=freq,
                       cell=pd$cell[i],
                       calltype=pd$type[i])%>%
                filter(cov>10)}
    })
    calls.tb = do.call(rbind,call.reads)
}
if (T){
    # plot
    plotpath = paste0(plotpre,".pdf")
    pdf(plotpath,useDingbats=F)
    for ( mod in unique(pd$type) ){
        dat = calls.tb[which(calls.tb$calltype==mod),]
        g = ggplot(dat,aes(x=mcall,fill=chrom,group=chrom,color=chrom))+
            facet_grid(cell~.,scales="free")+
            geom_density(alpha=0.5)+ggtitle(paste0(mod," mcalls"))+
            theme_bw()
        print(g)
        if (what == "perread"){
            g = ggplot(dat,aes(x=sign(score)*log(abs(score)),fill=chrom,group=chrom,color=chrom))+
                geom_density(alpha=0.5)+ggtitle(paste0(pd$cell[i]," : ",pd$type[i]," mcalls"))+
                theme_bw()
            print(g)
        }        
    }
    dev.off()
    plotpath = paste0(plotpre,"_ecdf.pdf")
    pdf(plotpath,useDingbats=F)
    for ( mod in unique(pd$type)) {
        dat = calls.tb[which(calls.tb$calltype==mod),]
        g = ggplot(dat,aes(x=mcall,group=chrom,color=chrom))+
            facet_grid(cell~.,scales="free")+
            stat_ecdf()+ggtitle(paste0(mod," mcalls"))+
            theme_bw()
        print(g)
    }
    dev.off()
}

