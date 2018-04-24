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
    # binning function
    binChr=function(x,bin=1000){
        mutate(x,
               start=floor(start/bin)*bin,
               end=ceiling(end/bin)*bin)
    }
    
    # let's find overlaps before binning
    require(GenomicRanges)
    db.gr=GRanges(db)
    meth.gr=GRanges(meth)
    #open.gr=db.gr[which(db.gr$cov==10)]
    ovl=findOverlaps(meth.gr,db.gr)
    open.meth=meth[queryHits(open.ovl),]
    summarize(open.meth,sum(called_sites_methylated)/sum(called_sites),n())
    # bin
    bin=100
    meth.bin.all=binChr(meth,bin)%>%
        group_by(chromosome,start,end)%>%
        summarize(calls=sum(called_sites),
                  meth=sum(called_sites_methylated),
                  freq=meth/calls)
    meth.bin=filter(meth.bin.all,calls>5)

    cov.bin.all=binChr(db,bin)%>%
        group_by(chromosome,start,end)%>%
        summarize(cov=sum(cov))
    cov.bin=filter(cov.bin.all,cov>=0)
    # tidying
    meth.td=ungroup(meth.bin)%>%
        transmute(chr=chromosome,
                  start=start,
                  end=end,
                  val=freq,
                  type="meth")
    cov.td=ungroup(cov.bin)%>%
        transmute(chr=chromosome,
                  start=start,
                  end=end,
                  val=cov,
                  type="cov")
    # combine
    dat.td=rbind(meth.td,cov.td)
    # spread
    dat.comp=spread(dat.td,key=type,value=val)%>%
        filter(!is.na(cov)&!is.na(meth))
    cor(dat.comp$cov,dat.comp$meth)
    cov.ecdf=ecdf(dat.comp$cov)
    meth.ecdf=ecdf(dat.comp$meth)
    dat.norm=mutate(dat.comp,
           covnorm=cov.ecdf(cov),
           methnorm=meth.ecdf(meth))
    cor(dat.norm$covnorm,dat.norm$methnorm)

    # density plot?
    g=ggplot(dat.comp)+theme_bw()
    gmeth=g+geom_density(aes(x=meth))
    gcov=g+geom_density(aes(x=cov))

    pdf(file.path(plotdir,"gpcvscov.pdf"),width=9,height=6)
    print(gmeth)
p    print(gcov)
    dev.off()
    
    
}

