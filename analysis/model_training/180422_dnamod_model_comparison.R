#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(tidyverse)
library(ggridges)
library(ggplot2)

plotdir="/home/isac/Dropbox/Data/Nanopore/171120_ecoliStandards/plots"
fpath="/home/isac/Code/nanopolish-isac/etc/r9-models/r9.4_450bps.cpg.6mer.template.model"
#/dilithium/Data/Nanopore/Analysis/171120_ecoliStandards/train/cpg/171124_ecoliCpG.0/r9.4_450bps.cpg.6mer.template.model"

dat.raw=read_tsv(fpath,comment="#",col_names=c("kmer","level"))

# subset kmers that have exactly one occurence of CG motif
dat.raw$count=str_count(dat.raw$kmer,"M")
dat=dat.raw%>%
    filter(count==1)
# M at the last position or CG motif
dat.mg=dat[sapply(dat$kmer,function(x){strsplit(x,'')[[1]][6]=="M"}),]
dat.last=dat[grep("MG",dat$kmer),]
dat.filt=bind_rows(dat.mg,dat.last)%>%
    select(-count)%>%
    mutate(mind=regexpr("M",kmer))
# get the unmeth
dat.comp=dat.filt%>%
    mutate(unmeth=gsub("M","C",kmer),
           baselevel=dat.raw$level[match(unmeth,dat.raw$kmer)],
           diff=level-baselevel)
# summarize by methylation index
dat.sum=dat.comp%>%
    group_by(mind)%>%
    summarize(mean=mean(diff))

# plot
g=ggplot(dat.comp,aes(x=diff,y=factor(mind)))+
    geom_density_ridges(alpha=0.1)+
    theme_bw()+xlim(c(-7.5,7.5))+
    labs(x="Unmethylated-Methylated (pA)",
         y="Position of 5m-C")

pdf(file.path(plotdir,"180422_kmer_distro.pdf"),width=9,height=6)
print(g)
dev.off()
