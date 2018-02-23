#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(tidyverse)

plotdir="/home/isac/Dropbox/Data/Nanopore/171205_nomeseq/plots"
methdir="/dilithium/Data/Nanopore/Analysis/171205_nomeseq/dam/methcall"
mod="dam"
if (F){
    # use meth tsvs
    suffix="meth.tsv"
    fpaths=system(intern=T,
                  command=paste0('find ',methdir," -maxdepth 1 -type f -name *",suffix))
    fnames=sapply(strsplit(fpaths,"/"),"[[",9)
    pd=as.tibble(do.call(rbind,strsplit(fnames,"[.]"))[,1:2])
    names(pd)=c("samp","calltype")
    methlist=lapply(fpaths,read_tsv)
    meth.tidy.list = lapply(seq_along(methlist),function(i){
        cbind(methlist[[i]],pd[i,])})
    meth.tidy=as.tibble(do.call(rbind,meth.tidy.list))
    # by read
    thr=2.5
    methfreq=group_by(meth.tidy,read_name,calltype)%>%
        summarize(meth=sum(log_lik_ratio>thr),
                  unmeth=sum(log_lik_ratio<(-thr)),
                  freq=meth/(meth+unmeth))%>%
        filter(meth>10|unmeth>10)
    methcomp=methfreq[,c("read_name","calltype","freq")]%>%
        spread(key=calltype,value=freq) %>%
        filter(!is.na(cpg)&!is.na(gpc))
    cor(methcomp$cpg,methcomp$gpc)
    # plot?
    g = ggplot(methcomp,aes(x=cpg,y=gpc))+theme_bw()+
        geom_point()

    pdf(file.path(plotdir,"cpg_vs_gpc_byread.pdf"),width=9,height=6)
    print(g)
    dev.off()
    # still bad
}

    
if (T){
    suffix="freq.tsv"
    fpaths=system(intern=T,
                  command=paste0('find ',methdir," -maxdepth 1 -type f -name *",suffix))
    fnames=sapply(strsplit(fpaths,"/"),"[[",9)

    # get pData from file name
    pd=as.tibble(do.call(rbind,strsplit(fnames,"[.]"))[,1:2])
    names(pd)=c("samp","calltype")

    methlist=lapply(fpaths,read_tsv)
    meth.tidy.list = lapply(seq_along(methlist),function(i){
        cbind(methlist[[i]],pd[i,])})
    meth.tidy=as.tibble(do.call(rbind,meth.tidy.list))

    # bin
    n=1000
    meth.bin=mutate(meth.tidy,
                    chr=chromosome,
                    start=round(start/n)*n,
                    end=round(end/n)*n)%>%
        group_by(calltype,samp,chr,start,end)%>%
        summarize(called=sum(called_sites),
                  meth=sum(called_sites_methylated),
                  freq=meth/called) %>%
        filter(called>10)

    # compare
    meth.comp=meth.bin[,c("calltype","chr","start","end","freq")]%>%
        spread(key=calltype,value=freq) %>%
        filter(!is.na(cpg)&!is.na(dam))
    # correlation
    cor(meth.comp$cpg,meth.comp$dam)
    #bad

    # plot?
    g = ggplot(meth.comp,aes(x=cpg,y=dam))+theme_bw()+
        geom_point()

    pdf(file.path(plotdir,"cpg_vs_dam.pdf"),width=9,height=6)
    print(g)
    dev.off()
}

