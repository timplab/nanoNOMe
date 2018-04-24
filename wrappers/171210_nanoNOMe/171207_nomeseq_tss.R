#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(tidyverse)

plotdir="/home/isac/Dropbox/Data/Nanopore/171205_nomeseq/plots"
methdir="/dilithium/Data/Nanopore/Analysis/171205_nomeseq/dam/methcall"
if (T){
    # use meth tsvs
    suffix="freq.tsv"
    fpaths=system(intern=T,
                  command=paste0('find ',methdir,
                                 " -maxdepth 1 -type f -name *.dam.*",suffix))
    ctypes=list(
        chromosome=col_character(),
        start=col_integer(),
        end=col_integer(),
        num_recsites_in_group=col_integer(),
        called_sites=col_integer(),
        methylated_frequency=col_double(),
        group_sequence=col_character()
    )
    meth=read_tsv(fpaths,col_types=ctypes)
    # need to make chr names compatible
    #meth=mutate(meth,chromosome=paste0("chr",chromosome))
    db="genes"
    #db="atac"
    if (db=="atac"){
    # wig file
        dbpath="/dilithium/Data/Nanopore/Analysis/171205_nomeseq/database/GSM2439558_MDA-MB-231-DMSO-4w.bedGraph"
        cnames=c("chromosome","start","end","cov")
        db=read_tsv(dbpath,col_names=cnames)
    }else if (db=="genes"){
        dbpath="/mithril/Data/NGS/Reference/human_annotations/Homo_sapiens.GRCh38.89.bed"
        cnames=c("chr","start","end","strand")
        db=read_tsv(dbpath,col_names=cnames)
        db=db%>%
            mutate(tss=ifelse(strand=="+",start,end))
    }
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
        filter(abs(dist)<=2000)
    b=10
    n=1
    meth.sum=mutate(dist.tb,bindist=round(dist/b)*b)%>%
        group_by(bindist)%>%
        summarize(calls=sum(call),
                  methcalls=sum(methcall),
                  freq=methcalls/calls)%>%
        filter(calls>n)

    # plot
    g=ggplot(meth.sum,aes(x=bindist,y=freq))+theme_bw()+
        geom_point()

    pdf(file.path(plotdir,"methylation_by_distance_to_tss.pdf"),width=9,height=6)
    print(g)
    dev.off()
    
    

}

