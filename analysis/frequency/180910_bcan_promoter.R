#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(ggridges)
library(GGally)
source("../../plot/methylation_plot_utils.R")
parser="../../nanopolish/parseMethylFreq.py"

# set this to TRUE to remove unnecessary objects throughout the process
limitedmem=TRUE
cores = detectCores()-2

# set directories
ref = "/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa"
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
plotdir=file.path(root,"plots/promoters")
if (!dir.exists(plotdir)) dir.create(plotdir,recursive=TRUE)
datroot=file.path(root,"pooled/methylation/mfreq_all")
regdir=file.path(root,"database/hg38/")
# get file paths of data
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
fpaths=tibble(cell=cells,
          cpg=file.path(datroot,paste0(cell,".cpg.methfreq.txt.gz")),
          gpc=file.path(datroot,paste0(cell,".gpc.methfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)

# regions
if (T){
    reg.fp = file.path(regdir,"hg38_genes.TSS.400bp.bed")
    regname="prom"
    subset=FALSE
}
#read in regions
cat("reading in the region\n")
extracols=c("TSS","TSSend")
db = load_db(reg.fp,extracols)

# getting regional methylation
# first test with 1k lines
coms = lapply(pd$filepath,function(x){
    paste0("head -n1000 ",reg.fp," | python ",parser," by-region -v -t 12 -i ",x)#," -r ",reg$filepath)
})
coms = lapply(pd$filepath,function(x){
    paste0("python ",parser," by-region -t ",cores," -i ",x," -r ",reg.fp)
})

cnames = c("chrom","start","end","id","score","strand",extracols,"freq","cov","numsites")
dat.raw = lapply(coms,function(x){
    print(x)
    y = system(x,intern=T)
    out = as.tibble(do.call(rbind,strsplit(y,"\t")))
    names(out) = cnames
    out
})

dat.named = lapply(seq_along(dat.raw),function(i){
    dat.raw[[i]]%>%type_convert()%>%
        mutate(cell=pd$cell[i],
               calltype=pd$calltype[i])
})

dat.all = do.call(rbind,dat.named)

# spread by cell 
dat.spread = dat.all %>%
    select(-TSSend,-cov,-numsites)%>%
    spread(cell,freq)
# get delta
dat.del = dat.spread %>%
    mutate(del=MCF10A-MCF7)

# spread by cpg and gpc
del.sum = dat.del %>%
    select(-GM12878,-MCF10A,-MCF7,-MDAMB231)%>%
    spread(calltype,del)%>%
    na.omit()

# order by cumulative difference (after -gpc)
sum.order = del.sum %>%
    mutate(gpc=-gpc,
           score=cpg+gpc)%>%
    arrange(desc(gpc))
    
# output bed
n = 100
bedout = file.path(plotdir,"mdaVS10a_promoters.bed")
write_tsv(sum.order[1:n,],bedout,col_names=F)
