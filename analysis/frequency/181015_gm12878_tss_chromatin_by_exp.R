#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(gridExtra)
library(GenomicRanges)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

root = commandArgs(trailingOnly=TRUE)[1]
# default path if not provided
if (is.na(root)){ root="/dilithium/Data/Nanopore/projects/nomeseq/analysis" }
cat(paste0("data root : ",root,"\n"))

# stratified by expression quartile
# data
datroot=file.path(root,"pooled/methylation/mfreq_all")
cell="GM12878"
lowercell="gm12878"
cpgpath=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz"))
gpcpath=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz"))
pd=tibble(cell=cell,type=c("cpg","gpc"),filepath=c(cpgpath,gpcpath))

# annotations
tss.fp=file.path(root,"annotations/hg38/hg38_genes.TSS.400bp.bed")    
cgi.fp=file.path(root,"annotations/hg38/hg38_cgi.txt.gz")

# output dir
plotdir=file.path(root,"plots/correlation")
if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)

#read in tss db
tss.gr=load_db(tss.fp)
tss.gr$id = sapply(strsplit(tss.gr$id,"[.]"),"[[",1)

# rnaseq database
exp.dir=file.path(root,"database",lowercell,"rnaseq")
exp.fp=c(file.path(exp.dir,paste0(cell,"_rnaseq.1.tsv")),
         file.path(exp.dir,paste0(cell,"_rnaseq.2.tsv")))

# expression data parsing
qtiles=seq(0,1,.25)
exp.list=lapply(exp.fp,read_tsv)
exp.list=lapply(seq_along(exp.list),function(i){
    x=exp.list[[i]]
    id=sapply(strsplit(x$gene_id,"[.]"),"[[",1)
    y=tibble(id=id,transcripts=x$"transcript_id(s)",fpkm=x$FPKM_ci_upper_bound)
    y=y[which(y$id %in% tss.gr$id),]
    qs=quantile(y$fpkm,probs=c(0.25,0.5,0.75))
    y%>%mutate(qtile=ifelse(fpkm>=qs[1],2,1),
               qtile=ifelse(fpkm>=qs[2],3,qtile),
               qtile=ifelse(fpkm>=qs[3],4,qtile),
               rep=i)
})
exp=do.call(rbind,exp.list)
exp.qtile=exp%>%select(-fpkm)%>%
    spread(rep,qtile)%>%filter(`1`==`2`)%>%mutate(qtile=`1`)
# using only genes that are consistent b/w the two replicates, taking average
exp.fpkm=exp[which(exp$id %in% exp.qtile$id),] %>%
    group_by(id,transcripts,qtile)%>%summarize(fpkm=mean(fpkm))

# merge with tss.gr
tss.gr=tss.gr[which(tss.gr$id %in% exp.fpkm$id)]
tss.gr$qtile=exp.fpkm$qtile[match(tss.gr$id,exp.fpkm$id)]
tss.gr$fpkm=exp.fpkm$fpkm[match(tss.gr$id,exp.fpkm$id)]

# read in data around TSS
dat.tss=lapply(seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],tss.fp)})

regmeth.list=lapply(seq_along(dat.tss),function(i){
    getRegionMeth(dat.tss[[i]],tss.gr)%>%
        mutate(calltype=pd$type[i],
               cell=pd$cell[i])
})
regmeth=do.call(rbind,regmeth.list)

dat.spread = regmeth%>%
    filter(numsites>10)%>%
    select(-totcov,-numsites)%>%
    spread(calltype,freq)%>%na.omit() 
dat.spread = bind_cols(dat.spread,
                       as.data.frame(tss.gr[dat.spread$feature.index]))

plotter <- function(dat.plt,plotpath) {
    # plot
    dat.plt = dat.plt %>% ungroup() %>%
        mutate(qtile=ifelse(qtile==2,1,qtile),
               qtile=ifelse(qtile==3,4,qtile))
    g = ggplot(dat.plt,aes(y=cpg,x=gpc))+theme_bw()+
        facet_wrap(~qtile)+lims(x=c(0,1),y=c(0,1))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"),
              strip.background = element_blank(),
              plot.title = element_text(size=10,hjust=0),
              panel.border = element_rect(colour = "black"))+
        labs(y="CpG methylation",x="GpC accessibility",
             title="Chromatin state by expression")
    # scatter plot
    g.scatter = g + geom_point(size=0.2,alpha=0.5) +
        geom_rug(size=0.1,alpha=0.5)
    g.2d = g + stat_density_2d(aes(fill=..level..),geom="polygon")
    # heatmap
    hist = bin_pairwise_methylation(dat.plt)
    g.heat = g +
        geom_tile(inherit.aes=F,data=hist,
                  mapping=aes(y=y,x=x,fill=count))+
        scale_fill_distiller(palette="Spectral")+
        theme(legend.position="right")
    # scatter plot in one
    g.one = ggplot(dat.plt,aes(y=cpg,x=gpc))+theme_bw()+
        geom_point(aes(color=factor(qtile)),size=0.5,alpha=0.3)+
        lims(x=c(0,1),y=c(0,1))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"),
              legend.position="bottom",
              plot.title = element_text(size=10,hjust=0))+
        labs(y="CpG methylation",x="GpC accessibility",
             title="Chromatin state by expression",
             color="quartile")
    # density
    dat.gather = gather(dat.plt,calltype,freq,-qtile)
    g.density = g +
        geom_density(inherit.aes=F,data=dat.gather,
                     mapping=aes(x=freq,y=..scaled..))+
        facet_grid(calltype~qtile) +
        labs(y="Density",x="Methylation frequency",
             title="Methylation density by expression")
    # print
    pdf(plotpath,useDingbats=F,height=3,width=6)
    print(g.one)
    print(g.2d)
    print(g.density)
    print(g.scatter)
    print(g.heat)
#    print(g.density)
    dev.off()
}
dat.plt=dat.spread %>%
    select(c("gpc","cpg","qtile"))

# subsample by min qtile
minnum = min(table(dat.plt$qtile))
dat.plt=dat.plt%>%group_by(qtile)%>%
    do(sample_n(.,minnum))
plotpath=file.path(plotdir,paste0(cell,"_CpGvsGpCvsExp_scatter.pdf"))
plotter(dat.plt,plotpath)

# cgi 
# cgi db
cgi.cnames=c("idx","chrom","start","end")
cgi = GRanges(read_tsv(cgi.fp,cgi.cnames))
tss.cgi = tss.gr[overlapsAny(tss.gr,cgi)]
cgi.id = tss.cgi$id

dat.plt=dat.spread[match(tss.cgi$id,tss.gr$id),] %>%
    na.omit()%>%
    select(c("gpc","cpg","qtile"))

# subsample by min qtile
minnum = min(table(dat.plt$qtile))
dat.plt=dat.plt%>%group_by(qtile)%>%
    do(sample_n(.,minnum))

plotpath=file.path(plotdir,"GM12878_CGI_CpGvsGpCvsExp_scatter.pdf")
plotter(dat.plt,plotpath)
# non-cgi
dat.plt = dat.spread[-match(tss.cgi$id,tss.gr$id),] %>%
    na.omit() %>%
    select(c("gpc","cpg","qtile"))
# subsample by min qtile
minnum = min(table(dat.plt$qtile))
dat.plt=dat.plt%>%group_by(qtile)%>%
    do(sample_n(.,minnum))
plotpath=file.path(plotdir,"GM12878_nonCGI_CpGvsGpCvsExp_scatter.pdf")
plotter(dat.plt,plotpath)
