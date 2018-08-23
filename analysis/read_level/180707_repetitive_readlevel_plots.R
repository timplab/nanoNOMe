#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(ggridges)
library(GGally)
source("../../plot/methylation_plot_utils.R")

# set this to TRUE to remove unnecessary objects throughout the process
limitedmem=TRUE

# set directories
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
plotdir=file.path(root,"plots/repeats")
if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
datroot=file.path(root,"pooled/methylation/methbyread_all")
regdir=file.path(root,"database/hg38")
# get file paths of data
cells=c("MCF10A","MCF7","MDAMB231")
pd=tibble(cell=rep(cells,each=2),
              calltype=rep(c("cpg","gpc"),length(cells)),
              filepath=file.path(datroot,paste(cell,calltype,"pooled.meth.bed.gz",sep=".")))

dbpath=file.path(plotdir,"LINE_regions.bed")
db = load_db(dbpath,c("one","two","cpg","gpc"))

# plot these regions
if (T){
    win=1000
    pltwin.gr = resize(db,width=width(db)+win,fix="center")
    dat.list = lapply(pd$filepath,function(x){
        tabix_mbed(x,pltwin.gr)
    })
    calls.list = lapply(seq_along(dat.list),function(i){
        dat.list[[i]][[2]]%>%mutate(cell=pd$cell[i],
                                    calltype=pd$calltype[i])
    })
    calls = do.call(rbind,calls.list)
    calls.gr = GRanges(calls)
    ymax=10
    plotpath=file.path(plotdir,"LINE_regions_readlevel.pdf")
    pdf(plotpath,useDingbats=F)
    for (i in seq_along(pltwin.gr)){
        reg.plt = pltwin.gr[i]
        reg = as.tibble(db[i])
        calls.reg = lapply(calls.list,function(x){
            y = x[overlapsAny(GRanges(x),reg.plt),]%>%
                na.omit()
            y$y = factor(y$qname)
            newlevels=seq(0,ymax,length.out=length(levels(y$y)))
            levels(y$y)=newlevels
            y$y = as.numeric(as.character(y$y))
            y
        })
        plt.tb = do.call(rbind,calls.reg)
        g = ggplot(calls.reg,aes(x=start,y=factor(qname),color=mcall))+
            facet_grid(cell~calltype)+
            geom_point(alpha=0.5,size=0.5)+
            geom_rect(inherit.aes=F,data=reg,alpha=0.2,
                      mapping=aes(xmin=start,xmax=end,
                                  ymin=0,ymax=ymax,
                                  fill=id))+
            theme_bw()+
            theme(axis.text.y=element_blank(),
                  panel.grid.major=element_line(size=0.2),
                  panel.grid.minor=element_line(size=0.1))
        print(g)
    }
    dev.off()
    
}
