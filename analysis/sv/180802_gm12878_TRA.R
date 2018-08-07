#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R") 
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/sv")
cell="GM12878"
pd = tibble(cell=cell,type=c("cpg","gpc"),
            filepath=file.path(datroot,paste(cell,type,"pooled.meth.bed.gz",sep=".")))
sv.fp = file.path(root,paste0("pooled/sv/",cell,".multiregion.bed"))

# read SVs
sv.gr = load_db(sv.fp,"rnames")
sv.centers = getCenter(sv.gr)

# load data around svs
dat.list = lapply(pd$filepath,function(x){
    tabix_mbed(x,sv.fp,by="read")})

# split data
cpg = dat.list[[which(pd$type=="cpg")]]
gpc = dat.list[[which(pd$type=="gpc")]]
# only use reads that have both cpg and gpc
reads.int = intersect(gpc$readname,cpg$readname)
cpg = cpg[which(cpg$readname %in% reads.int),]
gpc = gpc[which(gpc$readname %in% reads.int),]
cpg.gr = GRanges(cpg%>%select(chrom,start,end))
gpc.gr = GRanges(gpc%>%select(chrom,start,end))

# let's try with just one read
summary.tb = tibble()
for (i in seq(68,max(sv.gr$id))){
#    i = 1
    print(i)
    sv = sv.gr[which(sv.gr$id==i)]
    sv.ind = which(overlapsAny(gpc.gr,sv))
    sv.rnames = strsplit(sv$rnames[1],",")[[1]]
    # reads supporting sv
    gpc.sv = gpc.gr[which(gpc$readname%in%sv.rnames)]
    if (length(gpc.sv)==0) next
    # determine which is the inserted vs inserting region
    id.num = nearest(gpc.sv,sv.centers)
    id.sum = as.tibble(id.num) %>%group_by(value)%>%
        summarize(num=n())%>%arrange(desc(num))
    if (dim(id.sum)[1] != 2 ) next
#    if (dim(id.sum)[1]>2 | max(id.sum$num)==min(id.sum$num)) next # skip svs that doesn't have a clear origin vs target
    newloc = nearest(gpc.sv[which.max(width(gpc.sv))],sv.centers)#id.sum$value[1]
    gpc.new = gpc.sv[which(id.num==newloc)]
    # most lenient SV size
    # bp is on the right if range for ends is < rnage for starts
#    if (max(start(gpc.new))-min(start(gpc.new)) > max(end(gpc.new))-min(end(gpc.new))){
#        bp.tb = tibble(chrom=as.character(seqnames(gpc.new)[1]),
#                       start=end(gpc.new)[1],
#                       end = start(gpc.new)[length(gpc.new)])
#    } else {
    bp.tb = as.tibble(as.data.frame(sv.centers[newloc]))%>%select(seqnames,start,end)
    bp.gr = resize(GRanges(bp.tb),width=500,fix="center")
    gpc.insert = gpc.sv[which(id.num!=newloc)]
    insert.tb = tibble(chrom=as.character(seqnames(gpc.insert)[1]),
                       start = min(start(gpc.insert)),
                       end = max(end(gpc.insert)))
    insert.gr = GRanges(insert.tb)
    flanking.gr = resize(c(bp.gr,insert.gr),width=width(c(bp.gr,insert.gr))+1000,fix="center")
    # iterate through reads to determine what type
    gpc.svreg = gpc[sv.ind,]
    sv.sum = tibble(readname=character(),mcall=numeric(),
                    loglik=double(),sitenum=numeric(),calltype=character(),
                    allele=character(),svtype=character())
    for (rname in unique(gpc.svreg$readname)){
        calls.list = lapply(seq_along(dat.list),function(i){
            x = dat.list[[i]]
            mbedByCall(x[which(x$readname==rname),],verbose=F)%>%
                mutate(calltype=pd$type[i])
        })
        calls = do.call(rbind,calls.list)
        calls.gr = GRanges(calls)
        flank.ovl = findOverlaps(calls.gr,flanking.gr)
        calls$ref = NA
        calls$ref[queryHits(flank.ovl)] = subjectHits(flank.ovl)    
        calls$insert = overlapsAny(calls.gr,insert.gr)
        summary = group_by(calls,calltype,insert,ref)%>%
            summarize(mcall=mean(na.omit(mcall)),loglik=mean(score),sitenum=n())%>%
            mutate(allele="sv",
                   svtype=ifelse(insert,"insert",NA),
                   svtype=ifelse(ref==2&!insert,"origin",svtype),
                   svtype=ifelse(ref==1,"target",svtype),
                   readname=rname)%>%na.omit()%>%
            ungroup()%>%
            select(readname,mcall,loglik,sitenum,calltype,allele,svtype)
        if (rname %in% sv.rnames){
        # SV read
            summary = summary
        } else if (overlapsAny(insert.gr,calls.gr)) {
        # this read is the "original" insert region
            summary = summary %>%
                mutate(allele="nonsv")%>%
                filter(svtype %in% c("origin","insert"))
        } else if (overlapsAny(bp.gr,calls.gr)){
        # this is the original "target" region without the SV
            summary = summary %>%
                mutate(allele="nonsv")%>%
                filter(svtype=="target")
        }
        sv.sum = bind_rows(sv.sum,summary)
    }
    summary.tb = bind_rows(summary.tb,sv.sum %>%
                                      group_by(allele,svtype,calltype)%>%
                                      summarize(mcall=mean(na.omit(mcall)),
                                                loglik=mean(loglik),
                                                numreads=length(unique(readname)))%>%
                                      mutate(svid=i))
}

plt.tb = summary.tb %>%
    filter(numreads>=3)%>%
    select(-loglik,-numreads)%>%
    spread(allele,mcall)%>% na.omit()

g = ggplot(plt.tb,aes(x=nonsv,y=sv))+
    facet_grid(svtype~calltype)+
    geom_point(size=0.5,alpha=0.5)+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")
plotpath = file.path(plotdir,"GM12878_TRA_scatter_3x.pdf")
pdf(plotpath,width=6,height=4,useDingbats=FALSE)
print(g)
dev.off()







if (F){
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
}

