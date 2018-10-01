#!/usr/bin/Rscript
library(tidyverse)
library(parallel)

script = "./summarize_mappability.py"
root = "/dilithium/Data/Nanopore/projects/nomeseq/analysis"
samps = c("GM12878","GM12878_BSseq_1")
pd = tibble(samp = rep(samps,each=2),
            calltype = rep(c("cpg","gpc"),2),
            fpath = file.path(root,"mappability_bismap",paste(samp,calltype,"bismap.bedGraph",sep=".")))

qnum = 4
com = lapply(pd$fpath,function(x){
    paste("grep chr22",x,"|","python",script,"-q",qnum)
})

#com = com[[1]]
raw.list = mclapply(mc.cores=4,com,function(x){
    system(x,intern=T)
})
dat.list = lapply(raw.list,function(x){
    as.tibble(do.call(rbind,strsplit(x,"\t")))
        
