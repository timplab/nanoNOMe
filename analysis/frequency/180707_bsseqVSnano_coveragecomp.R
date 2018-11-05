#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(GenomicRanges)
library(parallel)
library(ggridges)
library(GGally)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

# set this to TRUE to remove unnecessary objects throughout the process
limitedmem=TRUE
cores = detectCores()-2

# set directories
root=commandArgs(trailingOnly=TRUE)[1]
if (is.na(root)){
    root="/dilithium/Data/Nanopore/projects/nomeseq/analysis" # default
}
plotdir=file.path(root,"plots/coverage")
if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
datroot=file.path(root,"pooled/methylation/mfreq_all")
cells=c("GM12878","GM12878_BSseq_1_destrand")
fpaths=tibble(cell=cells,
          cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")),
          gpc=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)

cgi = file.path(root,"annotations/hg38/hg38_cgi.bed")
shuffle = file.path(root,"annotations/hg38/hg38_shuffle.bed")
repeats = file.path(root,"annotations/hg38/hg38_repeats_SINE.bed")

lab = "SINE"
dbpath = repeats

dat.list = mclapply(mc.cores=cores,seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],dbpath,cov=0)
})
sum.list = lapply(dat.list,function(x){
    x%>% summarize(meancov=mean(cov),
                   totcount=n())})
dat.list = lapply(seq_along(dat.list),function(i){
    dat.list[[i]]%>%mutate(
                        cell=pd$cell[i],
                        calltype=pd$calltype[i],
                        calltype=ifelse(calltype=="cpg","CpG Methylation","GpC Accessibility"),
                        sample=ifelse(cell=="GM12878","nanoNOMe","BSseq"),
                        cov=cov/sum.list[[i]]$meancov)
})
count.list = lapply(seq_along(dat.list),function(i){
    dat.list[[i]] %>% group_by(sample,calltype,cov)%>%
        summarize(count=n())%>%
        mutate(count=count/sum.list[[i]]$totcount)
})

dat = do.call(rbind,dat.list)
dat.count = do.call(rbind,count.list)
dat.pair = dat %>%
    select(chrom,start,cov,sample,calltype) %>%
    spread(sample,cov)%>%na.omit()
pair.sub = sample_n(dat.pair,100000)

thr = quantile(dat$cov,.999)

# plot
g = ggplot(dat.count,aes(x=cov,y=count,color=sample,group=sample))+
    facet_grid(calltype~.,scales="free")+
    theme_bw()+
    lims(x=c(0,thr))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"),
          strip.background = element_blank(),
          panel.border = element_rect(colour="black"))+
    labs(x="Normalized coverage",y="Density")

g.line = g + geom_line()
g.hist = g + geom_histogram(mapping=aes(fill=sample),
                            stat="identity",alpha=0.5,
                            position="identity")
g.both = g.hist + geom_line()

g.box = ggplot(dat,aes(x=calltype,y=cov,fill=sample,color=sample))+
    lims(y=c(0,thr))+coord_flip()+
    geom_boxplot(outlier.color=NA)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour="black"))+
    labs(y="Normalized coverage")

g.pair = ggplot(pair.sub,aes(x=nanoNOMe,y=BSseq))+
    facet_grid(calltype~.)+
    geom_point(alpha=0.2)+
    geom_abline(slope=1)+
    theme_bw()+
    lims(x=c(0,thr),y=c(0,thr))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"),
          strip.background = element_blank(),
          panel.border = element_rect(colour="black"))

plotpath=file.path(plotdir,paste0("BSseq_vs_nano_cov_",lab,".pdf"))
pdf(plotpath,useDingbats=F,height=6,width=4)
print(g.box)
print(g.both)
print(g.pair)
dev.off()
