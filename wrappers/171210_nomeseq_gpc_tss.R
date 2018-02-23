#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(tidyverse)

plotdir="/home/isac/Dropbox/Data/Nanopore/171210_nomeseq/plots"
methdir="/dilithium/Data/Nanopore/Analysis/171210_nomeseq/171210_gpc_nome_mda231_200U/methcall/gpc-hg38"
if (T){
    # use meth tsvs
    suffix="freq.tsv"
    fpaths=system(intern=T,
                  command=paste0('find ',methdir,
                                 " -maxdepth 1 -type f -name *",suffix))
    ctypes=list(
        chromosome=col_character(),
        start=col_integer(),
        end=col_integer(),
        num_recsites_in_group=col_integer(),
        called_sites=col_integer(),
        methylated_frequency=col_double(),
        group_sequence=col_character()
    )
    methlist=lapply(fpaths,read_tsv,col_types=ctypes)
    meth.all=do.call(rbind,methlist)
    meth=mutate(meth.all,
                key=paste(chromosome,
                          start,
                          end,
                          num_recsites_in_group,
                          group_sequence,
                          sep=":"))%>%
        group_by(key)%>%summarize(called=sum(called_sites),
                                  meth=sum(called_sites_methylated),
                                  freq=meth/called)
    meth=mutate(meth,chr=strsplit(key,":")[[1]][1],
                start=strsplit(key,":")[[1]][2],
                end=strsplit(key,":")[[1]][3])

    meth=read_tsv(fpaths[1],col_types=ctypes)
    # need to make chr names compatible
    #meth=mutate(meth,chromosome=paste0("chr",chromosome))
    dbpath="/mithril/Data/NGS/Reference/human_annotations/Homo_sapiens.GRCh38.89.bed"
    cnames=c("chr","start","end","strand")
    db=read_tsv(dbpath,col_names=cnames)
    db=db%>%
        mutate(tss=ifelse(strand=="+",start,end))

    require(GenomicRanges)
    meth.gr=GRanges(meth)
    tss=transmute(db,chr=chr,start=tss,end=tss,strand=strand)
    tss.gr=GRanges(tss)
    dist=distanceToNearest(meth.gr,tss.gr)
    dist.tb=as.tibble(dist) %>%
        mutate(methloc=meth$start[queryHits],
               tssloc=tss$start[subjectHits],
               strand=tss$strand[subjectHits],
               call=meth$called_sites[queryHits],
               methcall=meth$called_sites_methylated[queryHits])%>%
        mutate(dist=ifelse(strand=="+",tssloc-methloc,methloc-tssloc))%>%
        filter(abs(dist)<=3000)
    b=10
    n=10
    meth.sum=mutate(dist.tb,bindist=round(dist/b)*b)%>%
        group_by(bindist)%>%
        summarize(calls=sum(call),
                  methcalls=sum(methcall),
                  freq=methcalls/calls)%>%
        filter(calls>n)

    # plot
    g=ggplot(meth.sum,aes(x=bindist,y=freq))+theme_bw()+
        geom_point(size=0.5)

    pdf(file.path(plotdir,"gpc_methylation_by_distance_to_tss.pdf"),width=10,height=4)
    print(g)
    dev.off()
    
    

}

