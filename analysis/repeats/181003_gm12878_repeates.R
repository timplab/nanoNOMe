#!/usr/bin/Rscript
library(tidyverse)
library(parallel)

root = "/dilithium/Data/Nanopore/projects/nomeseq/analysis"
plotdir = file.path(root,"plots/repeats")
samps = c("GM12878","GM12878_BSseq_1")
db="SINE"
pd = tibble(samp = rep(samps,each=2),
            calltype = rep(c("cpg","gpc"),2),
            fpath = file.path(root,"intersect",paste(samp,calltype,db,"methfreq.bedGraph",sep=".")))
cores = dim(pd)[1]

cnames=c("chrom","start","end","meth","unmeth","chrom_reg","start_reg","end_reg","regname","regid","strand","sine")
dat.list = mclapply(mc.cores=cores,seq_along(pd$fpath),function(i){
    read_tsv(pd$fpath[i],col_names=cnames)%>%
        mutate(samp=pd$samp[i],calltype=pd$calltype[i])
})
dat.tb = do.call(rbind,dat.list)

dat.sum = dat.tb %>%
    group_by(samp,calltype,chrom_reg,start_reg,end_reg,regname)%>%
    summarize(methsum=sum(meth),
              cov=methsum+sum(unmeth),
              freq=methsum/cov)%>%
    na.omit()

dat.cov = dat.sum %>%
    select(-methsum,-freq)%>%
    spread(samp,cov) %>%
    replace(is.na(.),0.01)

# write the elemetns that only occur in nanopore seq into a bed file
dat.nano = dat.cov[which(dat.cov$calltype=="cpg"),]%>%
    filter(GM12878>10 & GM12878_BSseq_1 < 1)
dat.bed = dat.nano %>% ungroup() %>%
    select(chrom_reg,start_reg,end_reg,regname)%>%
    arrange(chrom_reg,start_reg,end_reg)
bedpath=file.path(root,paste("database/hg38/hg38",db,"onlyNanopore.bed",sep="_"))
write_tsv(dat.bed,bedpath,col_names=F)
    

dat.pair = dat.cov[which(dat.cov$calltype=="cpg"),6:7]
plot.labs = names(dat.pair)
names(dat.pair) = c("one","two")
normalizer <- function(data){
    1
    #mean(data)/10000
}
dat.ma = dat.pair %>%
    mutate(one=one/normalizer(dat.pair$one),
           two=two/normalizer(dat.pair$two),
           m = log2(one/two) ,a = 1/2*log2(one*two))


dat.samp = dat.ma %>%
    sample_frac(1)

# plot scatter first
g = ggplot(dat.samp,aes(x=a,y=m))+theme_bw()+
    geom_point()
g2 = ggplot(dat.samp,aes(x=one,y=two))+theme_bw()+
    geom_point()+
    labs(x=plot.labs[1],y=plot.labs[2])


plotpath = file.path(plotdir,"gm12878_sine_coverage_bsseq_vs_nanoNOme.pdf")
pdf(plotpath,useDingbats=F)
print(g)
print(g2)
dev.off()
