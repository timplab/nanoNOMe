#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R") 
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/readlevel")
cell="GM12878"
pd = tibble(cell=cell,type=c("cpg","gpc"),
            filepath=file.path(datroot,paste(cell,type,"pooled.meth.bed.gz",sep=".")))
enhancer.fp = file.path(root,"database/hg38/hg38_enhancers.bed")
tss.fp = file.path(root,"database/hg38/hg38_genes.TSS.400bp.bed")
enhancer.gr = load_db(enhancer.fp)
tss.gr = load_db(tss.fp)

# load data around enhancers
cnames=c("chrom","start","end")
dat.list = lapply(pd$filepath,function(x){
    tabix_mbed(x,enhancer.fp,by="read")})

# overlap with tss
cpg = dat.list[[which(pd$type=="cpg")]]
gpc = dat.list[[which(pd$type=="gpc")]]
cpg.gr = GRanges(cpg)
tss.ovl = findOverlaps(cpg.gr,tss.gr)
# overlap with enhancers
en.ovl = findOverlaps(cpg.gr,enhancer.gr)
# combine ovl
ovl = tibble(readind = queryHits(tss.ovl),
             tssind = subjectHits(tss.ovl),
             enind = subjectHits(en.ovl)[match(readind,queryHits(en.ovl))])

# take n most frequent overlaps
ovl.sum = as.tibble(ovl)%>%
    group_by(tssind,enind)%>%
    summarize(num = n()) %>%
    arrange(desc(num))

n = 1
plotpath = file.path(plotdir,"GM12878_enhancer_promoter_readlevel.pdf")
pdf(plotpath,width=6,height=4,useDingbats=FALSE)
for (i in seq(n)){
    sum.sub = ovl.sum[i,]
    ovl.sub = ovl %>%
        filter(tssind == sum.sub$tssind &
               enind == sum.sub$enind)
    tss = tss.gr[sum.sub$tssind]
    enhancer = enhancer.gr[sum.sub$enind]
    tss.tb = as.tibble(as.data.frame(tss)) %>% select(-regstart,-regend)
    enhancer.tb= as.tibble(as.data.frame(enhancer))
    regs = bind_rows(tss.tb,enhancer.tb)
    plotreg = tibble(chrom=as.character(seqnames(tss)),
                 start=min(c(start(tss),start(enhancer))),
                 end=max(c(end(tss),end(enhancer))))
    cpg.reads = cpg[ovl.sub$readind,]
    gpc.reads = gpc[which(gpc$readname %in% cpg.reads$readname),]
    
    reads.list = list(cpg.reads,gpc.reads)
    dat.list = lapply(seq(2),function(i){
        if (i == 1){
            h = 50
        }else { h = 50 }
        mbedByCall(reads.list[[i]],smooth=TRUE,h=h)%>%
            mutate(calltype=pd$type[i])})
    dat.plt = do.call(rbind,dat.list)
    ylevels=sort(levels(factor(dat.plt$qname)))
    g = ggplot(dat.plt,aes(x=start,y=factor(qname),color=mcall))+
        facet_grid(calltype~.)+
        geom_point(alpha=0.5,size=0.5) +
        xlim(c(plotreg$start,plotreg$end))+
        geom_rect(inherit.aes=F,data=regs,alpha=0.2,
                  mapping=aes(xmin=start,xmax=end,
                              ymin=min(ylevels),ymax=max(ylevels),
                              fill=id))+
        theme_bw()+
        theme(axis.text.y=element_blank(),
              panel.grid.major = element_line(size=0.2),
              panel.grid.minor = element_line(size=0.1))
    print(g)
}
dev.off()

