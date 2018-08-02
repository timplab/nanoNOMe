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
enhancer.fp = file.path(root,"database/gm12878/enhancer/GM12878_enhancer_promoter_hg38.bed")
trans.fp = file.path(root,"database/hg38/hg38_transcripts.bed")
exp.fp = file.path(root,"database/gm12878/rnaseq/GM12878_rnaseq_transcript_quartiles.tsv")

# load data around enhancers
cnames=c("chrom","start","end")
dat.list = lapply(pd$filepath,function(x){
    tabix_mbed(x,enhancer.fp,by="read")})

# annotations
exp.tb = read_tsv(exp.fp)
enhancer.gr = load_db(enhancer.fp)
trans.gr = load_db(trans.fp,extracols=c("name","geneid"))
tss.gr = promoters(trans.gr,upstream=1,downstream=0)
# mach enhancer with transcript
enhancer.tb = as.tibble(as.data.frame(enhancer.gr))%>%
    mutate(seqnames=as.character(seqnames))
tss.tb = as.tibble(as.data.frame(tss.gr)) %>%
    mutate(seqnames=as.character(seqnames))
# filter for same chrom, close
interaction = bind_cols(enhancer.tb,
                        tss.tb[match(enhancer.tb$id,tss.tb$id),],
                        exp.tb[match(enhancer.tb$id,exp.tb$id),]) %>%
    na.omit() %>% filter(seqnames == seqnames1,
                         abs(end1-start) < 10000 )%>%
    distinct(geneid,.keep_all=TRUE)
# make granges for the entire window of interaction
interaction.gr = GRanges(interaction%>%
                         transmute(seqnames,strand,
                                   start=ifelse(start<start1,start,start1),
                                   end=ifelse(end<start1,start1,end)))

# split data
cpg = dat.list[[which(pd$type=="cpg")]]
gpc = dat.list[[which(pd$type=="gpc")]]
# overlap with interaction
cpg.gr = GRanges(cpg%>%select(chrom,start,end))
interaction.ovl = findOverlaps(cpg.gr,interaction.gr)
ovl.frac = pintersect(cpg.gr[queryHits(interaction.ovl)],interaction.gr[subjectHits(interaction.ovl)])
w=width(ovl.frac)/width(interaction.gr[subjectHits(interaction.ovl)])
ovl.all=interaction.ovl[which(w==1)]

# take n most frequent overlaps
ovl.sum = as.tibble(ovl.all) %>%
    group_by(subjectHits)%>%
    summarize(num = n()) %>%
    arrange(desc(num)) %>%
    filter(num>10)

n = 200
plotpath = file.path(plotdir,"GM12878_enhancer_promoter_readlevel.pdf")
pdf(plotpath,width=6,height=4,useDingbats=FALSE)
for (i in seq(n)){
    # subset regions
    sum.sub = ovl.sum[i,]
    ovl.sub = ovl.all[which(subjectHits(ovl.all)==sum.sub$subjectHits)]
    sub.int = interaction[sum.sub$subjectHits,]
    reglabs = c("enhancer",
                paste0(sub.int$name," : ",sub.int$quartile))
    regs = tibble(chrom=sub.int$seqnames[1],
                  start=c(sub.int$start,ifelse(sub.int$strand1=="-",sub.int$start1-1000,sub.int$start1)),
                  end=c(sub.int$end,ifelse(sub.int$strand1=="-",sub.int$start1,sub.int$start1+1000)),
                  strand=c(as.character(sub.int$strand),as.character(sub.int$strand1)),
                  name=factor(reglabs,levels=reglabs))
    sub.reg = interaction.gr[sum.sub$subjectHits]
    
    cpg.reads = cpg[queryHits(ovl.sub),]
    gpc.reads = gpc[which(gpc$readname %in% cpg.reads$readname),]
    cpg.reads = cpg.reads[which(gpc.reads$readname %in% gpc.reads$readname),]
    
    reads.list = list(cpg.reads,gpc.reads)
    mcall.list = lapply(seq(2),function(i){
        mbedByCall(reads.list[[i]],h=h)%>%
            mutate(calltype=pd$type[i])})
    dat.plt = do.call(rbind,mcall.list)
    ylevels=sort(levels(factor(dat.plt$qname)))
    offset=1000
    g = ggplot(dat.plt,aes(x=start,y=factor(qname),color=mcall))+
        facet_grid(calltype~.)+
        geom_point(alpha=0.5,size=0.5) +
        xlim(c(start(sub.reg)-offset,end(sub.reg)+offset))+
        geom_rect(inherit.aes=F,data=regs,alpha=0.2,
                  mapping=aes(xmin=start,xmax=end,
                              ymin=min(ylevels),ymax=max(ylevels),
                              fill=name))+
        geom_vline(xintercept=ifelse(regs$strand[2]=="-",regs$end[2],regs$start[2]))+
        theme_bw()+
        theme(axis.text.y=element_blank(),
              panel.grid.major = element_line(size=0.2),
              panel.grid.minor = element_line(size=0.1))
    print(g)
}
dev.off()
