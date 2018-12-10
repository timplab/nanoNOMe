#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(GGally)
source("../../script/methylation_plot_utils.R")

# data
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/mfreq_all")
cell="GM12878"
lowercell="gm12878"
cpgpath=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz"))
gpcpath=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz"))
pd=tibble(cell=cell,type=c("cpg","gpc"),filepath=c(cpgpath,gpcpath))
plotdir=file.path(root,"plots/regions")
dbname="tss"
dbname="ctcf"
#dbname="tss_shuffle"
#dbname="genebody"
if ( dbname=="tss"){
    db.fp=file.path(root,"database/hg38/hg38_genes.TSS.2000bp.bed")
    names.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names=read_tsv(names.fp,col_names=c("chrom","start","end","id","score","strand","name","fxn"))
    subset=FALSE
    # plot dir
    plotpath=file.path(plotdir,"GM12878_TSS_regions.pdf")
}else if (dbname=="ctcf"){
    db.fp=file.path(root,"annotations/gm12878/GM12878_CTCF_ctcfbsdb_allcomp.center.noTSS.2000bp.bed")
    plotpath=file.path(plotdir,"GM12878_CTCF_regions.pdf")
}else if (dbname == "tss_shuffle"){
    db.fp=file.path(root,"database/hg38/hg38_genes.TSS.2000bp.shuffle.bed")
    names.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names=read_tsv(names.fp,col_names=c("chrom","start","end","id","score","strand","name","fxn"))
}else if (dbname == "genebody") {
    db.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names.fp=file.path(root,"database/hg38/hg38_genes.bed")
    names=read_tsv(names.fp,col_names=c("chrom","start","end","id","score","strand","name","fxn"))
}

if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
#read in tss db
db.gr=load_db(db.fp)

# read in data around feature
dat.db=lapply(seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],db.gr[1:100]) %>%
        mutate(calltype=pd$type[i])
})
dat.tb = do.call(rbind,dat.db)
dat.gr = GRanges(dat.tb)
ovl = findOverlaps(dat.gr,db.gr)

win = 50

pdf(plotpath,width=4,height=3,useDingbats=F)
for (i in unique(subjectHits(ovl))){
    print(i)
    subidx = queryHits(ovl)[which(subjectHits(ovl)==i)]
    dat.sub = dat.tb[subidx,]
    dat.cpg = dat.sub[which(dat.sub$calltype=="cpg"),]
    dat.gpc = dat.sub[which(dat.sub$calltype=="gpc"),]
#    rm.cpg = calculate_rollmean(dat.cpg) %>% na.omit() %>%
#        mutate(calltype="cpg")
#    rm.gpc = calculate_rollmean(dat.gpc) %>% na.omit() %>%
#        mutate(calltype="gpc")
    dat.rm = do.call(rbind,list(rm.cpg,rm.gpc))
    g=ggplot(dat.sub,aes(x=start,y=freq))+
        facet_grid(calltype~.)+
        geom_line()+
        lims(y=c(0,1))+
        labs(x=paste0("Coordinate on ",dat.sub$chrom[1]),
             y="Methylation Frequency")+
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"))
#        geom_point(size=0.5,alpha=0.5)+        
    #        geom_smooth(se=F,method="loess",span=0.2)+
    print(g)
}
dev.off()
