#!/usr/bin/Rscript
library(tidyverse)
library(parallel)

root = "/dilithium/Data/Nanopore/projects/nomeseq/analysis"
plotdir = file.path(root,"plots/repeats")
samps = c("GM12878","GM12878_BSseq_1")
db="SINE"
pd = tibble(samp = samps,
            fpath=c("/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/bed/GM12878.pooled.SINE.bedGraph",
                    "/dilithium/Data/NGS/projects/gm12878/bsseq/bed/GM12878_BSseq_ENCLB794YYH.nodup.SINE.bedGraph"))
cores = dim(pd)[1]

cnames=c("chrom","start","end","name","score","strand","chrom_reg","start_reg","end_reg","regname","regid","strand","sine")
dat.list = mclapply(mc.cores=cores,seq_along(pd$fpath),function(i){
    read_tsv(pd$fpath[i],col_names=cnames)%>%
        mutate(samp=pd$samp[i])
})
dat.tb = do.call(rbind,dat.list)

dat.sum = dat.tb %>%
    group_by(samp,chrom_reg,start_reg,end_reg,regname)%>%
    summarize(cov=n())

dat.cov = dat.sum %>%
    spread(samp,cov) %>%
    replace(is.na(.),0.1)

# write the elemetns that only occur in nanopore seq into a bed file
dat.nano = dat.cov%>%
    filter(GM12878>10 & GM12878_BSseq_1 < 1)
dat.bed = dat.nano %>% ungroup() %>%
    select(chrom_reg,start_reg,end_reg,regname)%>%
    arrange(chrom_reg,start_reg,end_reg)
bedpath=file.path(root,paste("database/hg38/hg38",db,"onlyNanoporebam.bed",sep="_"))
write_tsv(dat.bed,bedpath,col_names=F)    

dat.pair = dat.cov[,5:6]
plot.labs = names(dat.pair)
names(dat.pair) = c("one","two")
normalizer <- function(data){
    mean(data)
}
dat.ma = dat.pair %>%
    mutate(one=one/normalizer(dat.pair$one),
           two=two/normalizer(dat.pair$two),
           m = log2(one/two) ,a = 1/2*log2(one*two))
ma.ordered = dat.ma%>%
    arrange(a)
ma.ordered = ma.ordered[round(seq(1,dim(ma.ordered)[1],len=5000)),]
fit = loess(ma.ordered$m~ma.ordered$a)
bias = predict(fit,newdata=data.frame(dat.ma$a))

dat.samp = dat.ma %>%
    sample_frac(1)

# plot scatter first
g = ggplot(dat.samp,aes(x=a,y=m))+
    geom_point()
g2 = ggplot(dat.samp,aes(x=one,y=two))+
    geom_point()+
    labs(x=plot.labs[1],y=plot.labs[2])


plotpath = file.path(plotdir,paste("gm12878",db,"bamcoverage_bsseq_vs_nanoNOMe.pdf",sep="_"))
pdf(plotpath,useDingbats=F)
print(g)
print(g2)
dev.off()
