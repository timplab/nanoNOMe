#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(tidyverse)

plotdir="/home/isac/Dropbox/Data/Nanopore/171210_nomeseq/plots"
methdir="/dilithium/Data/Nanopore/Analysis/171210_nomeseq/171210_gpc_nome_mda231_200U/methcall"
if (T){
    # use meth tsvs
    suffix="*.tsv"
    region="wnt7bTSS"
    tss=45977128  
    fpaths=system(intern=T,
                  command=paste0('find ',methdir,
                                 " -maxdepth 1 -type f -name *",region,suffix))
    ctypes=list(
        chromosome=col_character(),
        start=col_integer(),
        end=col_integer(),
        readname=col_character(),
        loglikratio=col_double()
    )
    cnames=c("chromosome","start","end","readname","loglikratio")
    methlist=lapply(fpaths,read_tsv,col_names=F)
    reads=unique(methlist[[1]]$X4)

    pd=sapply(strsplit(sapply(strsplit(fpaths,"/"),"[[",9),"[.]"),"[[",3)
    
    meth.tidy.list=lapply(seq_along(pd),function(i){
        mutate(methlist[[i]],calltype=pd[i])})
    meth.tidy=do.call(rbind,meth.tidy.list)    

    meth.list=lapply(seq_along(reads),function(i){
        filter(meth.tidy,X4==reads[i])%>%
            transmute(chromosome=X1,
                      start=X2,
                      end=X3,
                      ind=i,
                      calltype=calltype,
                      loglikratio=X5,
                      seq=X10)
    })
    meth=do.call(rbind,meth.list)
    thr=2.5
    meth=filter(meth,loglikratio<(-thr)|loglikratio>thr)%>%
        filter(!grepl("GCG",seq))%>%
        mutate(methcall=ifelse(loglikratio>thr,"meth","unmeth"),
               callid=paste0(calltype,methcall))%>%
        mutate(y=ifelse(calltype=="cpg",ind+0.2,ind-0.2))
    xlims=c(min(meth$start),max(meth$start))
    ylims=c(min(meth$y),max(meth$y))
    chr=as.character(meth[1,1])

    g=ggplot(meth,aes(x=start,y=y,color=callid))+theme_bw()+
        geom_jitter(alpha=0.7,height=0.1)+geom_vline(xintercept=tss)+
        labs(title=region,x=paste0(chr," Coordinate"))
    g.cpg=ggplot(filter(meth,calltype=="cpg"),aes(x=start,y=y,color=callid))+
        theme_bw()+
        geom_point(alpha=0.7)+xlim(xlims)+ylim(ylims)
    g.gpc=ggplot(filter(meth,calltype=="gpc"),aes(x=start,y=y,color=callid))+
        theme_bw()+
        geom_point(alpha=0.7)+xlim(xlims)+ylim(ylims)
   
    pdf(file.path(plotdir,paste0(region,"gpcvsmeth_readlevel.pdf")),width=9,height=6)
    print(g)
    print(g.cpg)
    print(g.gpc)
    dev.off()
    
    
}

