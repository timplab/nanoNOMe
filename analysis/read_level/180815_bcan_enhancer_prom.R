#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R")
cores = round(detectCores()/4*3) # using 3/4 of available cores
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/readlevel")
cells=c("MCF10A","MCF7","MDAMB231")
pd = tibble(cell=rep(cells,2),
            calltype=rep(c("cpg","gpc"),times=length(cells)),
            filepath=file.path(datroot,paste(cell,calltype,"pooled.meth.bed.gz",sep=".")))
enhancer.fp = file.path(root,"database/hg38/hg38_enhancers.bed")
trans.fp = file.path(root,"database/hg38/hg38_transcripts.bed")
#exp.fp = file.path(root,"database/gm12878/rnaseq/GM12878_rnaseq_transcript_quartiles.tsv")

# load data around enhancers
dat.list = mclapply(mc.cores=cores,pd$filepath,function(x){
    tabix_mbed(x,enhancer.fp,by="read")})
dat.list = lapply(seq_along(dat.list),function(i){
    x = dat.list[[i]]
    x %>% mutate(cell=pd$cell[i],
                 calltype=pd$calltype[i])})
dat.all = do.call(rbind,dat.list)
# annotations
enhancer.gr = load_db(enhancer.fp)
trans.tb = as.tibble(as.data.frame((load_db(trans.fp,extracols=c("name","geneid"))))) %>%
    mutate(id=name)%>%select(-name,-geneid)%>%unique()
trans.gr = GRanges(trans.tb)
tss.gr = unique(promoters(trans.gr,upstream=1,downstream=0))
prom.gr = unique(promoters(trans.gr,upstream=500,downstream=0))
body.gr = unique(promoters(trans.gr,upstream=1,downstream=500))
# match enhancer with transcript
dists.list = mclapply(mc.cores=cores,seq_along(enhancer.gr),function(i){
    d = distanceToNearest(tss.gr,enhancer.gr[i])
    as.tibble(d) %>% mutate(subjectHits=i) %>%
        filter(distance<5000)
})
dists = do.call(rbind,dists.list) %>%
    rename(queryHits="tx",subjectHits="enhancer")
interactions.list = mclapply(mc.cores=cores,seq(dim(dists)[1]),function(i){
    interaction = c(enhancer.gr[dists$enhancer[i]],
      tss.gr[dists$tx[i]],
      prom.gr[dists$tx[i]],
      body.gr[dists$tx[i]])
    interaction$lab = c("Enhancer","TSS","Promoter","Body")
    interaction
})
windows.list = mclapply(mc.cores=cores,interactions.list,function(x){
    GRanges(seqnames(x)[1],IRanges(start=min(start(x))-500,end=max(end(x))+500),id = x$id[2])
})
windows = do.call(c,windows.list)
if (T) {
    windows.tb = as.tibble(as.data.frame(windows))
    bedpath = file.path(plotdir,"enhancerwindow.bed")
    windows.bed = windows.tb %>%
        transmute(chrom=seqnames,
               start=start,
               end=end,
               name=id,
               score=".",
               strand=".")
    write_tsv(windows.bed,bedpath,col_names=F)
    
}
# get reads in the interaction windows
ovl = findOverlaps(windows,GRanges(dat.all),type="within")
dat.ovl = dat.all[subjectHits(ovl),]%>%mutate(interaction = queryHits(ovl))

# filter out interactions with low read num
ovl.sum = dat.ovl %>%
    group_by(interaction,cell,calltype)%>%
    summarize(num=n())%>%
    filter(num>=10)%>%
    ungroup()%>%group_by(interaction)%>%
    summarize(samples=n(),num=mean(num))%>%
    filter(samples==dim(pd)[1])%>%
    arrange(desc(num))

ovl.list = lapply(seq(dim(pd)[1]),function(i){
    dat.ovl[which(dat.ovl$cell == pd$cell[i] & dat.ovl$calltype == pd$calltype[i]),]
})
calls.list = lapply(seq_along(ovl.list),function(i){
    mbedByCall(ovl.list[[i]],verbose=F)%>%
        mutate(cell=pd$cell[i],calltype=pd$calltype[i])
})

getfreq = function(x,i,w){
        x[which(x$interaction==i),]%>%
            mbedByCall(verbose=F)%>%
            na.omit()%>%
            group_by(chrom,start,end)%>%
            summarize(freq=mean(mcall),
                      cov=n())%>%
            filter(start>=start(windows[i]),
                   end<end(windows[i]))
}

freqregs.list = mclapply(mc.cores=cores,ovl.sum$interaction,function(i){
    print(i)
    int = interactions.list[[i]]
    # get freq within the window
    freq.list  = mclapply(mc.cores=cores,ovl.list,function(x)getfreq(x,i,windows))
    # get freq in specific parts
    regfreq.list = mclapply(mc.cores=cores,seq_along(freq.list),function(j){
        x=freq.list[[j]]
        freq.ovl = findOverlaps(GRanges(x),int)
        regfreq=x[queryHits(freq.ovl),]
        regfreq$reg = int$lab[subjectHits(freq.ovl)]
        y = regfreq%>%group_by(reg)%>%
            summarize(totcov=sum(cov),
                      freq=sum(freq*cov)/totcov) %>%
            add_row(reg="All",totcov=sum(x$cov),freq=sum(x$freq*x$cov)/sum(x$cov))%>%
            mutate(cell=pd$cell[j],calltype=pd$calltype[j])
    })
    do.call(rbind,regfreq.list)%>%mutate(idx=i)
})
freq.regs = do.call(rbind,freqregs.list)
freq.spread = freq.regs %>% select(-totcov)%>%
    spread(cell,freq)%>%mutate(del=abs(MDAMB231-MCF10A))%>%
    select(reg,calltype,idx,del)%>%
    spread(reg,del)%>%arrange(desc(Enhancer))

if (T) {
    sig.list = lapply(freq.spread$idx,function(i){interactions.list[[i]][1]})
    sig = do.call(c,sig.list)%>%unique()
    bedpath = file.path(plotdir,"significant_enhancers.bed")
    sig.bed = as.tibble(as.data.frame(sig)) %>%
        transmute(chrom=seqnames,
                  start=start,
                  end=end,
                  name=id,
                  score=score,
                  strand=".")
    write_tsv(sig.bed,bedpath,col_names=F)
}

plotpath = file.path(plotdir,"180815_bcan_enhancer_promoter_readlevel.pdf")
ymax=10
pdf(plotpath,width=6,height=4,useDingbats=FALSE)
for (i in seq(dim(freq.spread)[1])){
    print(i)
    # subset regions
    idx = freq.spread$idx[i]
    int = as.tibble(as.data.frame(interactions.list[[idx]]))
    int$lab[which(int$lab=="Body")] = int$id[2]
    reg.plt = windows.list[[idx]]
    calls.reg = lapply(calls.list,function(x){
        y = x[overlapsAny(GRanges(x),reg.plt),]%>%
            na.omit()
        y$y = factor(y$qname)
        newlevels = seq(0,ymax,length.out=length(levels(y$y)))
        levels(y$y) = newlevels
        y$y = as.numeric(as.character(y$y))
        y
    })
    plt.tb = do.call(rbind,calls.reg)
    g = ggplot(plt.tb,aes(x=start,y=y,color=mcall))+
        facet_grid(cell~calltype)+
        geom_point(alpha=0.5,size=0.5)+
        geom_rect(inherit.aes=F,data=int,alpha=0.2,
                  mapping=aes(xmin=start,xmax=end,
                              ymin=0,ymax=ymax,
                              fill=lab))+
        theme_bw()+
        theme(axis.text.y=element_blank(),
              panel.grid.major=element_line(size=0.2),
              panel.grid.minor=element_line(size=0.1))
    print(g)
}
dev.off()
