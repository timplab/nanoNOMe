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
res_no10a <- results(dds,contrast = c("cell","MDAMB231","MCF7"))


# make tibble and comibne
res_mcf7.tb <- as_tibble(res_mcf7,rownames = "id")
res_mda.tb <- as_tibble(res_mda,rownames = "id")
res_no10a.tb <- as_tibble(res_no10a,rownames = "id")
res.tb <- bind_rows(list(MCF7 = res_mcf7.tb, MDAMB231 = res_mda.tb, no10A = res_no10a.tb),.id = "what")
res.tb <- res.tb %>%
  mutate(
    one = ifelse(what != "no10A", "MCF10A","MCF7"),
    two = ifelse(what != "no10A", what, "MDAMB231")) %>%
  dplyr::select(-what)
out_path <- "/kyber/Data/Nanopore/projects/nanonome/analysis/data/bcan/bcan_rnaseq_diffexp_DESeq2.txt"
write_tsv(res.tb,out_path)

# significant ones?
sig <- res.tb %>%
  filter(padj < 0.01, abs(log2FoldChange) > 1.5) %>%
  mutate(comparison = paste0(one,"_vs_",two) ) %>%
  dplyr::select(id,comparison,log2FoldChange) %>%
  spread(comparison,log2FoldChange) %>%
  replace(is.na(.),0) %>%
  mutate(id = sapply(strsplit(id,"[.]"),"[[",1))

# merge with annotation
genes_fp <- "/kyber/Data/Nanopore/projects/nanonome/analysis/data/hg38/hg38_genes.bed"
genes_all <- read_tsv(genes_fp, col_names = c("chrom","start","end","ensid","score","strand","hgnc")) %>%
  mutate(start = start + 1) # convert it from bed to normal format
genes <- genes_all %>%
  arrange(chrom,start,end) %>%
  distinct(ensid, .keep_all =T)
genes$ensid <- sapply(strsplit(genes$ensid,"[.]"),"[[",1)

# merge the data
genes.sig <- genes %>%
  bind_cols(sig[match(genes$ensid,sig$id),]) %>%
  na.omit() %>%
  dplyr::select(-id)

out_path <- "/uru/Data/Nanopore/projects/nanonome/bcan/bcan_diffexp_sig_genes.txt"
write_tsv(genes.sig,out_path)
