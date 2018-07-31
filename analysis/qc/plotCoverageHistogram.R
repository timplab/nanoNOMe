#!/usr/bin/Rscript
testinput="/dilithium/Data/Nanopore/projects/nomeseq/analysis/coverage/GM12878.cpg.bedcov.bed"
input.fp=testinput
library(tidyverse)
args=commandArgs(TRUE)

input.fp = args[-length(args)]
outpath = args[length(args)]

cnames=c("loc","depth","count","regsize","frac")
dat = lapply(input.fp,function(x){
    read_tsv(x,col_names=cnames,comment="chr")})

# get numbers
tsvpath = paste0(outpath,".summary.tsv")
sum.list = lapply(dat,function(x){
    x = x%>%mutate(cumfrac=cumsum(x$frac))
    qtile.ind = sapply(seq(0.25,0.75,0.25),function(y){
        which(x$cumfrac>y)[1]})
    frac10 = x$cumfrac[which(x$depth==10)]
    qtile = x[qtile.ind,]$depth
    as.tibble(t(c(qtile,frac10)))
})
sum.tb = do.call(rbind,sum.list) %>%
    mutate(basename(input.fp))
names(sum.tb) = c("Q1","Median","Q3","Fraction<10x","Filename")
write_tsv(sum.tb,tsvpath)

# plot
plotHist <- function(data){
    # filter out outliers
    dat.filt = data %>%
        filter(depth < data$depth[which.max(data$count)] * 3)
    g = ggplot(dat.filt,aes(x=depth,y=frac))+
        geom_histogram(stat="identity")+
        theme_bw()
    g
}
pdf(outpath,width=5,height=4,useDingbats=F)
for (i in seq_along(input.fp)){
    print(plotHist(dat[[i]])+ggtitle(basename(input.fp[i])))
}
dev.off()
