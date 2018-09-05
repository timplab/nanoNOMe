#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R") 
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/sv")
cells=c("MCF10A","MCF7","MDAMB231")
pd = tibble(cell=rep(cells,each=2),type=rep(c("cpg","gpc"),times=length(cells)),
            filepath=file.path(datroot,paste(cell,type,"pooled.meth.bed.gz",sep=".")))
sv.fp = file.path(root,paste0("pooled/sv/",cells,".TRAmultiregion.bed"))

# read SVs
sv.list = lapply(sv.fp,function(x){load_db(x,"rnames")})
centers.list = lapply(sv.list,getCenter)
# concatenate the sv list
sv.list = lapply(seq_along(sv.fp),function(i){
    x = sv.list[[i]]
    x$cell=cells[i];x})
sv.all = do.call(c,sv.list)

# load data around svs
dat.list = lapply(seq_along(pd$filepath),function(i){
    x = pd$filepath[i]
    s = sv.fp[which(cells %in% pd$cell[i])]
    print(x)
    print(s)
    tabix_mbed(x,s,by="read")
})


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
