#!/usr/bin/Rscript
library(tidyverse)
library(DESeq2)

# rna-seq data
rnafp = "/kyber/Data/Nanopore/projects/nanonome/analysis/data/bcan/bcan_rnaseq.txt"
rna = read_tsv(rnafp)
names(rna)[8:10] = c("MDAMB231_R1","MDAMB231_R2","MDAMB231_R3")
samples = c("MCF10A","MCF7","MDAMB231")
# let's use DESeq2
# gene names
rowdata <- rna[[1]]
# data values
data.mat <- as.data.frame(rna[,-1])
rownames(data.mat) <- rowdata
# remove weird stuff at the end
keepi <- grepl("ENS",rowdata)
rowdata <- rowdata[keepi]
data.mat <- data.mat[keepi,]
# coldata
ids <- names(data.mat)
cells <- factor(sapply(strsplit(ids,"_"),"[[",1))
reps <- factor(sapply(strsplit(ids,"_"),"[[",2))
coldata <- tibble(id = ids, cell = cells, rep = reps)
# make dds 
design <- ~ cell
dds <- DESeqDataSetFromMatrix(countData = data.mat, colData = coldata, design = design)
dds <- DESeq(dds)
resultsNames(dds)
res_mcf7 <- results(dds,name = "cell_MCF7_vs_MCF10A")
res_mda <- results(dds,name = "cell_MDAMB231_vs_MCF10A")

# make tibble and comibne
res_mcf7.tb <- as_tibble(res_mcf7,rownames = "id")
res_mda.tb <- as_tibble(res_mda,rownames = "id")
res.tb <- bind_rows(list(MCF7 = res_mcf7.tb, MDAMB231 = res_mda.tb),.id = "two")
out_path <- "/kyber/Data/Nanopore/projects/nanonome/analysis/data/bcan/bcan_rnaseq_diffexp_DESeq2.txt"
write_tsv(res.tb,out_path)

# below is old version : just using t-stats
#combos = as.tibble(t(combn(samples,2)))
#names(combos) = c("one","two")
## get t-stat
#rna.pval = tibble(geneid=character(),one=character(),
#                  two=character(),meanone = numeric(), 
#                  meantwo = numeric(),p=numeric(),
#                  direction=character())
#for (i in seq(dim(combos)[1])){
#    print(i)
#    one = combos$one[i]
#    two = combos$two[i]
#    cat(paste(one,two,"\n"))
#    dat.one = rna[,grep(one,names(rna))]
#    dat.two = rna[,grep(two,names(rna))]
#    for (j in seq(dim(rna)[1])){
#        print(j)
#        one.j = as.numeric(dat.one[j,])
#        two.j = as.numeric(dat.two[j,])
#        meanone <- mean(one.j)
#        meantwo <- mean(two.j)
#        logfold = log(meantwo/meanone)
#        pval <- try(t.test(log(one.j+1),log(two.j+1))$p.value,TRUE)
#        pval <- ifelse(is.numeric(pval),pval,NA)
#        direction = ifelse(logfold<0,"under","over") # two in comparison to one
#        rna.pval = bind_rows(rna.pval,
#                            tibble(geneid = rna$Ensembl_ID[j],
#                                   one = one,two = two, meanone = meanone, meantwo = meantwo,
#                                   p = pval,logfold = logfold,direction=direction))
#    }
#}
## thresholding
#rna.pval$adjusted.pval <- p.adjust(rna.pval$p,"BH")
#rna.pval
#
#out_path <- "/kyber/Data/Nanopore/projects/nanonome/analysis/data/bcan/bcan_rnaseq_diffexp.txt"
#write_tsv(rna.pval,out_path)
