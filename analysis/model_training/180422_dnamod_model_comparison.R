#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files
library(tidyverse)
library(ggridges)
library(ggplot2)

plotdir="/home/isac/Dropbox/Data/nome-seq/reproduce/plots"
pd = tibble(mod=c("cpg","gpc"),
            fpath=paste0("/home/isac/Code/nanopolish/etc/r9-models/r9.4_450bps.",
                         mod,".6mer.template.model"))

raw.list=lapply(seq_along(pd$fpath),function(i){
    read_tsv(pd$fpath[i],comment="#",col_names=c("kmer","level"))%>%
        mutate(mod=pd$mod[i],
               lab=paste0(kmer,mod))
})
dat.raw = do.call(rbind,raw.list)

# subset kmers that have exactly one occurence of CG motif
dat.raw$count=str_count(dat.raw$kmer,"M")
dat.meth=dat.raw%>%
    filter(count==1)

# cpg
dat.mg = dat.meth[grep("MG",dat.meth$kmer),]
dat.lastm = dat.meth[sapply(dat.meth$kmer,function(x){strsplit(x,'')[[1]][6]=="M"}),]
meth.cpg = bind_rows(dat.mg,dat.lastm)%>%
    filter(mod=="cpg")%>%
    select(-count)%>%
    mutate(mind=regexpr("M",kmer))
# gpc
dat.gm = dat.meth[grep("GM",dat.meth$kmer),]
dat.firstm = dat.meth[sapply(dat.meth$kmer,function(x){strsplit(x,'')[[1]][1]=="M"}),]
meth.gpc = bind_rows(dat.gm,dat.firstm)%>%
    filter(mod=="gpc")%>%
    select(-count)%>%
    mutate(mind=regexpr("M",kmer))
#combine
meth = bind_rows(meth.cpg,meth.gpc)

# get the unmeth
dat.comp=meth%>%
    mutate(unmeth=gsub("M","C",lab),
           baselevel=dat.raw$level[match(unmeth,dat.raw$lab)],
           diff=level-baselevel)

# summarize by methylation index
dat.sum=dat.comp%>%
    group_by(mind,mod)%>%
    summarize(mean=mean(diff))

# plot
g = ggplot(dat.comp)+theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"))
g.ridges = g + geom_density_ridges(aes(x=diff,y=factor(mind)),alpha=0.1)+
    facet_grid(mod~.)+
    xlim(c(-10,10))+
    labs(x="Unmethylated-Methylated (pA)",
         y="Position of 5m-C")
g.box = g + geom_boxplot(aes(x=factor(mind),y=diff,color=mod),outlier.size=0.5)+
    ylim(c(-10,10))+
    labs(y="Unmethylated-Methylated (pA)",
         x="Position of 5m-C")

pdf(file.path(plotdir,"methylation_model_kmer_distro.pdf"),useDingbats=T,width=4,height=4)
print(g.ridges)
print(g.box)
dev.off()
