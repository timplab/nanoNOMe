#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(GenomicRanges)
library(parallel)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

# set this to TRUE to remove unnecessary objects throughout the process
limitedmem=TRUE
cores = detectCores()-2

# set directories
root=commandArgs(trailingOnly=TRUE)[1]
if (is.na(root)){
    root="/dilithium/Data/Nanopore/projects/nomeseq/analysis" # default
}
outdir=file.path(root,"annotations/breastcancer")
datroot=file.path(root,"pooled/methylation/mfreq_all")
cells=c("MCF10A","MCF7","MDAMB231")
fpaths=tibble(cell=cells,
          cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")),
          gpc=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)

# hg38 gene db
dbpath = file.path(root,"annotations/hg38/hg38_genes.bed")
db.gr = load_db(dbpath,extracols=c("genename","fxn"))
prom.gr = promoters(db.gr,200,200)

# cgi db
cgipath = file.path(root,"annotations/hg38/hg38_cgi.bed")
cgi.gr = load_db(cgipath)

# rna-seq data
rnafp = file.path(root,"annotations/breastcancer/bcan_rnaseq.txt")
rna = read_tsv(rnafp)
rna$p7 = rna$p231 = 1
# get t-stat
mindiff = 100
for (i in seq(dim(rna)[1])){
    print(i)
    m10 = as.numeric(rna[i,c(2:4)])
    mean10 = mean(m10)
    m7 = as.numeric(rna[i,c(5:7)])
    mean7 = mean(m7)
    mda = as.numeric(rna[i,c(8:10)])
    mean231 = mean(mda)
    if (abs(mean231-mean10) > mindiff) {
        rna$p231[i] = try(t.test(log(m10+1),log(mda+1))$p.value)
    }
    if (abs(mean7-mean10) > mindiff) {
        rna$p7[i] = try(t.test(log(m10+1),log(m7+1))$p.value)
    }
}
# thresholding
a = 0.05
rna.filt = rna %>%
    filter(p7<=a|p231<=a) %>%
    arrange(p7)

# match names
# remove versions in ids
db.id = sapply(strsplit(db.gr$id,"[.]"),"[[",1)
rna.id = sapply(strsplit(rna.filt$Ensembl_ID,"[.]"),"[[",1)
matchidx = match(rna.id,db.id)
rna.filt$name = db.gr$genename[matchidx]

db.rna = db.gr[na.omit(matchidx)]
db.rna$p231 = rna.filt$p231[-which(is.na(matchidx))]
db.rna$p7 = rna.filt$p7[-which(is.na(matchidx))]
prom.rna = prom.gr[na.omit(matchidx)]

# looking into CGI promoters?
upstream.rna = promoters(db.rna)
ovl.idx = overlapsAny(upstream.rna,cgi.gr,minoverlap=500)
rnacgi.gr = upstream.rna[ovl.idx]
prom.rnacgi = prom.rna[ovl.idx]

# methyltion data
dat.list = mclapply(mc.cores=cores,seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],prom.rnacgi)
})

regmeth.list = lapply(seq_along(dat.list),function(i){
    getRegionMeth(dat.list[[i]],db.gr) %>%
        mutate(cell=pd$cell[i],
               calltype=pd$calltype[i])
})
regmeth = do.call(rbind,regmeth.list)

meth.spread = regmeth %>% select(-totcov,-numsites) %>%
    spread(cell,freq) %>%
    replace(.,is.na(.),0)%>%
    mutate(del231=MDAMB231-MCF10A,
           del7=MCF7-MCF10A)

# sort by methylation difference
del231 = meth.spread %>%
    select(feature.index,calltype,del231) %>%
    spread(calltype,del231) %>%
    mutate(delsum = abs(cpg)+abs(gpc)) %>%
    arrange(desc(abs(gpc)),desc(delsum)) %>% na.omit()
del7 = meth.spread %>%
    select(feature.index,calltype,del7) %>%
    spread(calltype,del7) %>%
    mutate(delsum = abs(cpg)+abs(gpc)) %>%
    arrange(desc(abs(gpc)),desc(delsum)) %>% na.omit()

top231 = prom.gr[unique(head(del231$feature.index,n=50))]
top231$comparison = "MDA231vsMCF10A"
top7 = prom.gr[unique(head(del7$feature.index,n=50))]
top7$comparison = "MCF7vsMCF10A"
top.gr = unique(c(top231,top7))
bed.tb = GRangesTobed(top.gr)
outpath=file.path(outdir,"bcan_rnaseq_methylation.CGI.TSS.400bp.bed")
write_tsv(bed.tb,outpath,col_names=F)
outpath=file.path(outdir,"bcan_rnaseq_methylation.CGI.TSS.5000bp.bed")
top.5kb = GRangesTobed(resize(top.gr,width=5000,fix="center"))
write_tsv(top.5kb,outpath,col_names=F)
