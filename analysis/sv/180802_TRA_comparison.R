#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R")
cores=10
parser="../../nanopolish/SVmethylation.py"
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/sv")
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
pd = tibble(cell=cells,
            cpgpath=file.path(datroot,paste(cell,"cpg.pooled.meth.bed.gz",sep=".")),
            gpcpath=file.path(datroot,paste(cell,"gpc.pooled.meth.bed.gz",sep=".")),
            bampath=file.path(root,"pooled/bam",paste0(cell,".pooled.bam")),
            svpath=file.path(root,paste0("pooled/sv/",cell,".sniffles.vcf")))

coms = lapply(seq(dim(pd)[1]),function(i){
    paste("python",parser,"-v -t 10 -b",pd$bampath[i],
          "-c",pd$cpgpath[i],"-g",pd$gpcpath[i],
          "-s",pd$svpath[i])
})
#coms = lapply(seq(dim(pd)[1]),function(i){
#    paste("tail -n100",pd$svpath[i],"| python",parser,"-v -t 10 -b",pd$bampath[i],
#          "-c",pd$cpgpath[i],"-g",pd$gpcpath[i])
#})
raw.list = lapply(coms,function(x){
    system(x,intern=T)})
cnames = c("chrom1","start1","end1","chrom2","start2","end2","readname","svid","strand1","strand2","tag","cpgcnt","gpccnt","cpgcov","gpccov")
dat.list = mclapply(mc.cores=cores,seq_along(raw.list),function(i){
    x = as.tibble(do.call(rbind,strsplit(raw.list[[i]],"\t")))
    colnames(x) = cnames
    x %>% mutate(cell=pd$cell[i])
})
dat.all = do.call(rbind,dat.list) %>% type_convert()

dat.sum = dat.all%>%group_by(cell,svid,tag)%>%
    summarize(numreads=n(),
              totcpgcov=sum(cpgcov),
              cpg = sum(cpgcnt)/totcpgcov,
              totgpccov=sum(gpccov),
              gpc = sum(gpccnt)/totgpccov)%>%
    filter(totcpgcov>5)

# compare sv vs non-sv
sum.insert = dat.sum[grepl("insert",dat.sum$tag),]
insert.cpg = sum.insert%>%select(cell,svid,tag,cpg) %>%
    spread(tag,cpg) %>% na.omit() %>% mutate(calltype="cpg")
insert.gpc = sum.insert%>%select(cell,svid,tag,gpc) %>%
    spread(tag,gpc) %>% na.omit() %>% mutate(calltype="gpc")
insert.plt = do.call(rbind,list(insert.cpg,insert.gpc))
sum.target = dat.sum[grepl("target",dat.sum$tag),]
target.cpg = sum.target%>%select(cell,svid,tag,cpg) %>%
    spread(tag,cpg) %>% na.omit() %>% mutate(calltype="cpg")
target.gpc = sum.target%>%select(cell,svid,tag,gpc) %>%
    spread(tag,gpc) %>% na.omit() %>% mutate(calltype="gpc")
target.plt = do.call(rbind,list(target.cpg,target.gpc))

# compare target-insert SV vs nonSV
sum.sv = dat.sum[!grepl("noSV",dat.sum$tag),]%>%ungroup()
sv.cpg = sum.sv%>%select(cell,svid,tag,cpg)%>%
    spread(tag,cpg)%>%na.omit()%>%
    transmute(cell,svid,del=target-insert,calltype="cpg",tag="SV")
sv.gpc = sum.sv%>%select(cell,svid,tag,gpc)%>%
    spread(tag,gpc)%>%na.omit()%>%
    transmute(cell,svid,del=target-insert,calltype="gpc",tag="SV")
sum.nosv = dat.sum[grepl("noSV",dat.sum$tag),]%>%ungroup()
nosv.cpg = sum.nosv%>%select(cell,svid,tag,cpg)%>%
    spread(tag,cpg)%>%na.omit()%>%
    transmute(cell,svid,del=target_noSV-insert_noSV,calltype="cpg",tag="noSV")
nosv.gpc = sum.nosv%>%select(cell,svid,tag,gpc)%>%
    spread(tag,gpc)%>%na.omit()%>%
    transmute(cell,svid,del=target_noSV-insert_noSV,calltype="gpc",tag="noSV")
del = do.call(rbind,list(sv.cpg,sv.gpc,nosv.cpg,nosv.gpc))
del.plt = del %>% spread(tag,del)%>%na.omit()
    

# plotting
g = ggplot(target.plt,aes(x=target_noSV,y=target))+
    facet_grid(cell~calltype)+
    geom_point(size=0.5,alpha=0.5)+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")+
    theme_bw()
g.insert = ggplot(insert.plt,aes(x=insert_noSV,y=insert))+
    facet_grid(cell~calltype)+
    geom_point(size=0.5,alpha=0.5)+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")+
    theme_bw()
g.del = ggplot(del.plt,aes(x=noSV,y=SV))+
    facet_grid(cell~calltype)+
    geom_point(size=0.5,alpha=0.5)+
    geom_abline(slope=1,intercept=0,size=0.3,alpha=0.7,linetype="dashed")+
    geom_hline(yintercept=0,size=0.3,alpha=0.7,linetype="dotted")+
    theme_bw()
box.plt = del %>% mutate(tag = paste(cell,tag,sep=":"))
g.del.box = ggplot(box.plt,aes(x=tag,y=del,color=calltype))+
    geom_boxplot()+
    theme_bw()


plotpath = file.path(plotdir,"TRA_scatter.pdf")
pdf(plotpath,width=6,height=4,useDingbats=FALSE)
print(g)
print(g.insert)
print(g.del)
print(g.del.box)
dev.off()
