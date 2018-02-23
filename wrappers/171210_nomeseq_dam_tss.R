#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(tidyverse)

plotdir="/home/isac/Dropbox/Data/Nanopore/171210_nomeseq/plots"
methdir="/dilithium/Data/Nanopore/Analysis/171210_nomeseq/171210_dam_nome_mda231_100U/methcall/dam-hg38"
if (T){
    # use meth tsvs
    suffix="freq.*tsv"
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
        group_sequence=col_character(),
        num_reads=col_integer()
    )
    meth=read_tsv(fpaths,col_types=ctypes)
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
        filter(abs(dist)<=2000)
    b=100
    n=10
    meth.sum=mutate(dist.tb,bindist=round(dist/b)*b)%>%
        group_by(bindist)%>%
        summarize(calls=sum(call),
                  methcalls=sum(methcall),
                  freq=methcalls/calls)%>%
        filter(calls>n)

    # plot
    g=ggplot(meth.sum,aes(x=bindist,y=freq))+theme_bw()+
        geom_point()

    pdf(file.path(plotdir,"dam100U_methylation_by_distance_to_tss.pdf"),width=9,height=6)
    print(g)
    dev.off()
    
    

}

