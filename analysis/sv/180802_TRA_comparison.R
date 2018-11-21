#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../script/methylation_plot_utils.R")
cores=10
parser="../../script/SVmethylation.py"
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/sv")
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
pd = tibble(cell=cells,
            cpgpath=file.path(datroot,paste(cell,"cpg.pooled.meth.bed.gz",sep=".")),
            gpcpath=file.path(datroot,paste(cell,"gpc.pooled.meth.bed.gz",sep=".")),
            bampath=file.path(root,"pooled/bam",paste0(cell,".pooled.bam")),
            svpath=file.path(root,paste0("pooled/sv/",cell,".pooled*.vcf.gz")))

coms = lapply(seq(dim(pd)[1]),function(i){
    paste("gunzip -c ",pd$svpath[i]," | python",parser,"-v -t ",cores," -b",pd$bampath[i],
          "-c",pd$cpgpath[i],"-g",pd$gpcpath[i])
})
#coms = lapply(seq(dim(pd)[1]),function(i){
#    paste("tail -n100",pd$svpath[i],"| python",parser,"-v -t 10 -b",pd$bampath[i],
#          "-c",pd$cpgpath[i],"-g",pd$gpcpath[i])
#})
raw.list = lapply(coms,function(x){
    system(x,intern=T)})
cnames = c("chrom1","start1","end1","chrom2","start2","end2",
           "readname","svid","strand1","strand2",
           "tag","location","cpgcnt","gpccnt","cpgcov","gpccov")
dat.list = mclapply(mc.cores=cores,seq_along(raw.list),function(i){
    x = as.tibble(do.call(rbind,strsplit(raw.list[[i]],"\t")))
    colnames(x) = cnames
    x %>% mutate(cell=pd$cell[i])
})

plotpath = file.path(plotdir,"TRA_scatter_sig.pdf")
dat.all = do.call(rbind,dat.list) %>% type_convert()
# randomize check
random = TRUE
if (random == TRUE) {
    n = dim(dat.all)[1]
    idx = sample(n,n)
    dat.all$tag = dat.all$tag[idx]
    plotpath = file.path(plotdir,"TRA_scatter_shuffle.pdf")
}

dat.sum = dat.all%>%group_by(cell,svid,tag,location)%>%
    summarize(numreads=n(),
              totcpgcov=sum(cpgcov),
              cpg = sum(cpgcnt)/totcpgcov,
              totgpccov=sum(gpccov),
              gpc = sum(gpccnt)/totgpccov)%>%
    filter(totcpgcov>10)
dat.sum = dat.sum %>% ungroup() %>%
    mutate(lab=paste(tag,location,sep="."))

# compare destination-origin SV vs nonSV
sum.sv = dat.sum[dat.sum$tag %in% "SV",]%>%ungroup()
sum.nonsv = dat.sum[dat.sum$tag %in% "nonSV",]%>%ungroup()
sv.cpg = sum.sv%>%select(cell,svid,location,cpg)%>%
    spread(location,cpg)%>%na.omit()%>%
    transmute(cell,svid,del=destination-origin,calltype="cpg",tag="SV")
sv.gpc = sum.sv%>%select(cell,svid,location,gpc)%>%
    spread(location,gpc)%>%na.omit()%>%
    transmute(cell,svid,del=destination-origin,calltype="gpc",tag="SV")
nonsv.cpg = sum.nonsv%>%select(cell,svid,location,cpg)%>%
    spread(location,cpg)%>%na.omit()%>%
    transmute(cell,svid,del=destination-origin,calltype="cpg",tag="nonSV") 
nonsv.gpc = sum.nonsv%>%select(cell,svid,location,gpc)%>%
    spread(location,gpc)%>%na.omit()%>%
    transmute(cell,svid,del=destination-origin,calltype="gpc",tag="nonSV")
# thresholding loci based on nonSV del
cpg.thr = 0 #as.numeric(quantile(abs(nonsv.cpg$del),.5))
gpc.thr = 0 #as.numeric(quantile(abs(nonsv.gpc$del),.5))
nonsv.cpg.filt = nonsv.cpg %>%
    filter(abs(del)>cpg.thr)
nonsv.gpc.filt = nonsv.gpc %>%
    filter(abs(del)>gpc.thr)

del = do.call(rbind,list(sv.cpg,sv.gpc,nonsv.cpg.filt,nonsv.gpc.filt))
del.spread = del %>% spread(tag,del) %>% na.omit()

# destination nonSV - origin SV vs origin nonSV - destination SV
sum.x = dat.sum %>% 
    filter(lab == "nonSV.origin" | lab == "SV.destination") %>%
    select(cell,svid,lab,cpg,gpc)
sum.y = dat.sum %>% 
    filter(lab == "nonSV.destination" | lab == "SV.origin") %>%
    select(cell,svid,lab,cpg,gpc)
x.cpg = sum.x %>% select(-gpc) %>%
    spread(lab,cpg) %>% na.omit() %>%
    transmute(cell,svid,calltype="cpg",
              var="x",
              del=nonSV.origin-SV.destination)
x.gpc = sum.x %>% select(-cpg) %>%
    spread(lab,gpc) %>% na.omit() %>%
    transmute(cell,svid,calltype="gpc",
              var="x",
              del=nonSV.origin-SV.destination)
y.cpg = sum.y %>% select(-gpc) %>%
    spread(lab,cpg) %>% na.omit()%>%
    transmute(cell,svid,calltype="cpg",
              var="y",
              del=nonSV.destination-SV.origin)
y.gpc = sum.y %>% select(-cpg) %>%
    spread(lab,gpc) %>% na.omit() %>%
    transmute(cell,svid,calltype="gpc",
              var="y",
              del=nonSV.destination-SV.origin)
comp.tb = bind_rows(x.cpg,x.gpc,y.cpg,y.gpc)
comp.plt = spread(comp.tb,var,del) %>% na.omit()

# plotting
g.del.cpg = ggplot(del.spread[which(del.spread$calltype=="cpg"),],aes(x=nonSV,y=SV)) +
    geom_point()+
    lims(x=c(-1,1),y=c(-1,1))+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")+
    geom_hline(yintercept=0,size=0.3,alpha=0.7,linetype="dotted")+
    theme_bw()
g.del.gpc = ggplot(del.spread[which(del.spread$calltype=="gpc"),],aes(x=nonSV,y=SV)) +
    geom_point()+
    lims(x=c(-1,1),y=c(-1,1))+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")+
    geom_hline(yintercept=0,size=0.3,alpha=0.7,linetype="dotted")+
    theme_bw()
g.comp.cpg = ggplot(comp.plt[which(comp.plt$calltype=="cpg"),],aes(x=x,y=y)) +
    geom_point()+lims(x=c(-0.5,0.5),y=c(-0.5,0.5))+
    lims(x=c(-1,1),y=c(-1,1))+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")+
    geom_hline(yintercept=0,size=0.3,alpha=0.7,linetype="dotted")+
    theme_bw()
g.comp.gpc = ggplot(comp.plt[which(comp.plt$calltype=="gpc"),],aes(x=x,y=y)) +
    geom_point()+lims(x=c(-0.3,0.3),y=c(-0.3,0.3))+
    lims(x=c(-1,1),y=c(-1,1))+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")+
    geom_hline(yintercept=0,size=0.3,alpha=0.7,linetype="dotted")+
    theme_bw()



pdf(plotpath,width=4,height=4,useDingbats=FALSE)
print(g.del.cpg)
print(g.del.gpc)
print(g.comp.cpg)
print(g.comp.gpc)
dev.off()
