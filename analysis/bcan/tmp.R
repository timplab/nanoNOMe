#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))
library(parallel)
cores=detectCores() - 2
# root
root=commandArgs(trailingOnly=TRUE)[1]
if (is.na(root)){
    root="/dilithium/Data/Nanopore/projects/nomeseq/analysis" # default
}

# data
datroot=file.path(root,"pooled/methylation/mfreq")
cells=c("MCF10A","MCF7","MDAMB231")
fpaths=tibble(cell=cells,
              cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")),
              gpc=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)
plotdir=file.path(root,"plots/bcan_expression")
annodir=file.path(root,"annotations/breastcancer")

# get the db
if (FALSE){
    dbpaths = system(paste0("readlink -f ",
                            annodir,
                            "/*by_expression*TSS.bed"),
                     intern=T)
    extracols=c("gene","one","two","direction")
    db.list = lapply(dbpaths,load_db,extracols)
    db.gr = do.call(c,db.list)
    plotpath = file.path(plotdir,"breastcancer_diff_expression_regions_2kb.pdf")
}
if (TRUE){
    dbpath = file.path(annodir,"bcan_diffexp_bothcancer_cgi.TSS.bed")
    extracols="name"
    db.gr = load_db(dbpath,extracols)
    plotpath = file.path(plotdir,"breastcancer_diff_expression_bothcancer_regions_average_10kb.pdf")
}

db.gr = db.gr[order(db.gr$score,decreasing=T)]
reg.gr = promoters(db.gr,upstream=5000,downstream=5000)



load.bed <- function(bedfile) {
    
    ##Load region of interest from bed file
    reg.bed=read.delim(gzfile(bedfile), sep="\t", header=F)
    reg.gr=GRanges(seqnames=reg.bed[,1], ranges=IRanges(start=reg.bed[,2]+1, end=reg.bed[,3]))
    
    return(reg.gr)
}

cg.gr=load.bed("/home/timp/cgisl_hg38.bed.gz")
    
plotter <- function(data,reg,tss,calltype){
    g=ggplot(data,aes(x=start,y=freq,color=cell,group=cell))+
        geom_point(size=0.5,alpha=0.5)+
        geom_rug(inherit.aes=F,mapping=aes(x=start))+
        geom_vline(xintercept=start(tss),linetype=2)+
        lims(y=c(0,1),x=c(start(reg),end(reg)))+
        labs(x=paste0("Coordinate on ",seqnames(reg)),
             y="Methylation Frequency",
             title=paste(calltype,reg$id,"(",strand(tss),")",reg$score))+
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"),
              legend.position="bottom")
    g
}    

good.cgi=subsetByOverlaps(cg.gr, reg.gr)    


islmeth=tibble(mcf10=rep(-1.0, length(good.cgi)), mcf7=rep(-1.0, length(good.cgi)),mdamb231=rep(-1.0, length(good.cgi)))

for (i in seq_along(good.cgi)){   
    islmeth$mcf10[i]=mean(tabix_mfreq(pd$filepath[(pd$cell=="MCF10A"& pd$calltype=="cpg")], good.cgi[i])$freq)
    islmeth$mcf7[i]=mean(tabix_mfreq(pd$filepath[(pd$cell=="MCF7"& pd$calltype=="cpg")], good.cgi[i])$freq)
    islmeth$mdamb231[i]=mean(tabix_mfreq(pd$filepath[(pd$cell=="MDAMB231"& pd$calltype=="cpg")], good.cgi[i])$freq)
}

islmeth=islmeth %>%
    mutate(diff7=mcf7-mcf10) %>%
    mutate(diff231=mdamb231-mcf10) %>%
    mutate(unmeth.normal=mcf10<.2) %>%
    mutate(big.diff7=diff7 > .2) %>%
    mutate(big.diff231=diff231 >.2) %>%
    mutate(sum.diff=abs(diff7)+abs(diff231))

best.isl=islmeth %>%
    mutate(idx=1:n()) %>%
    arrange(-sum.diff)

pdf("/home/timp/try.pdf", width=4, height=3, useDingbats=F)
for (i in best.isl$idx) {
    print(i)
    myisl=good.cgi[i]
    myisl=myisl+3e3

    dat.list=lapply(seq(dim(pd)[1]),function(i){
        tabix_mfreq(pd$filepath[i],myisl) %>%
            mutate(calltype=pd$calltype[i],
                   cell=pd$cell[i])
    })
    dat.sub = do.call(rbind,dat.list)

    dat.cpg = dat.sub[which(dat.sub$calltype=="cpg"),]
    dat.gpc = dat.sub[which(dat.sub$calltype=="gpc"),]

    g.cpg = plotter(dat.cpg,myisl,good.cgi[i],"cpg") +
        geom_smooth(se=F,method="loess",span=0.2)
    g.gpc = plotter(dat.gpc,myisl,good.cgi[i],"gpc") +
        geom_smooth(se=F,method="loess",span=0.1)
    print(g.cpg)
    print(g.gpc)
}

dev.off()

system("scp /home/timp/try.pdf duchess:/home/timp/Dropbox/Temp/try.pdf")



##ok, let's try to also look at chromatin in the +/- 200 from tss

sm.prom=db.gr+200
tss.chr=tibble(gpc.mcf10=rep(-1.0, length(sm.prom)), gpc.mcf7=rep(-1.0, length(sm.prom)),
               gpc.mdamb231=rep(-1.0, length(sm.prom)))

for (i in seq_along(sm.prom)){   
    tss.chr$gpc.mcf10[i]=mean(tabix_mfreq(pd$filepath[(pd$cell=="MCF10A"& pd$calltype=="gpc")], sm.prom[i])$freq)
    tss.chr$gpc.mcf7[i]=mean(tabix_mfreq(pd$filepath[(pd$cell=="MCF7"& pd$calltype=="gpc")], sm.prom[i])$freq)
    tss.chr$gpc.mdamb231[i]=mean(tabix_mfreq(pd$filepath[(pd$cell=="MDAMB231"& pd$calltype=="gpc")], sm.prom[i])$freq)
}
