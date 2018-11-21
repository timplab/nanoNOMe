#!/usr/bin/Rscript
library(tidyverse)
library(ggridges)
# script is for plotting current distribution using subset eventalign files

root=commandArgs(trailingOnly=TRUE)[1]
if (is.na(root)){
    root="/dilithium/Data/Nanopore/projects/nomeseq/analysis" # default
}
plotdir = file.path(root,"plots/model")
eventdir = file.path(root,"data/ecoli_methylation")

mods=c("UM","CpG","GpC","CpGGpC")
pd=tibble(mod=mods,
          fpath=file.path(eventdir,paste0("ecoli_",mod,".eventalign")))
# read data    
eventlist=lapply(pd$fpath,read_tsv)
# tidy up 
event.tidy.list=lapply(seq_along(eventlist),function(i){
    transmute(eventlist[[i]],
              kmer=model_kmer,
              eventlevel=event_level_mean,
              mod=rep(pd$mod[i]))
})
event.tidy=do.call(rbind,event.tidy.list)


# don't use group aesthetic
g=ggplot(event.tidy,aes(x=eventlevel,y=kmer,fill=mod,color=mod))+
    geom_density_ridges(alpha=0.1,scale=1)+
    xlim(c(95,120))+
    labs(x="Event Level (pA)",
         y="Density Distribution") +
    theme_bw()+
    theme(legend.position="bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"))
          

pdf(file.path(plotdir,"ecoli_kmer_eventcomparison.pdf"),width=4,height=3)
print(g)
dev.off()
