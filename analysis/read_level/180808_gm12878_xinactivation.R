#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R") 
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
plotdir=file.path(root,"plots/readlevel")
cell="GM12878"
what = "perread"

# x chrom vs chrom 21,22
start = 0
end = 1e9
chrom.gr = GRanges(c("chr21","chr22","chrX"),
                   IRanges(start=c(rep(start,2),67676484),
                               end=c(rep(end,2),76800000)))
chrom.gr = GRanges(c("chr21","chr22","chrX"),
                   IRanges(start=c(rep(start,3)),
                               end=c(rep(end,3))))

if (what == "perread"){
    datroot=file.path(root,"pooled/methylation/methbyread_all")
    plotpre=file.path(root,"plots/readlevel/180808_GM12878_XvsAuto")
    pd = tibble(cell=rep(c(cell,"GM12878_wgs"),each=2),type=rep(c("cpg","gpc"),2),
                filepath=file.path(datroot,paste(cell,type,"pooled.meth.bed.gz",sep=".")))
#    chrom.gr = GRanges("chrM:0-10000000")
    # load data around these chroms
    cnames = c("chrom","start","end","qname","score","strand","mcall","sitenum")
    dat.list = mclapply(mc.cores=4,seq_along(pd$filepath),function(i){
        x = pd$cell[i]
        mod = pd$type[i]
        cat(paste0(x,"-",mod,"\n"))
        x.list = mclapply(mc.cores=3,seq_along(chrom.gr),function(j){
            y = as.character(seqnames(chrom.gr[j]))
            cat(paste0(x,"-",mod," : ",y,"\n"))
            command=paste("tabix",pd$filepath[i],y,"|",
                          "python ../../nanopolish/parseMethylbed.py per_read -m",mod,sep=" ")
            commandout = system(command,intern=TRUE)
            outin = as.tibble(do.call(rbind,strsplit(commandout,"\t")))
            numcol=c("V2","V3","V5","V7","V8")
            outin[numcol] = sapply(outin[numcol],as.numeric)
            cat(paste0(x,"-",mod," : ",y," done\n"))
            outin
        })
        cat(paste0(x,"-",mod," joining\n"))
        out = do.call(rbind,x.list)
        colnames(out) = cnames
        cat(paste0(x,"-",mod," done\n"))
        out
    })
    calls.list = mclapply(mc.cores=4,seq_along(dat.list),function(i){
        x = dat.list[[i]]
        if ("tbl" %in% class(x)){
            x%>%filter(sitenum>10)%>%
                mutate(cell=pd$cell[i],calltype=pd$type[i])
        }
    })
    calls.tb = do.call(rbind,calls.list)
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

