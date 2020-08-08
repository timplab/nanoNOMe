#!/usr/bin/Rscript
library(tidyverse)
library(DESeq2)
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R") 

outdir <- "~/Dropbox/Data/nome-seq/version3_guppy3/plots/bcan"
# rna-seq data
rnafp = "/kyber/Data/Nanopore/projects/nanonome/analysis/data/bcan/bcan_rnaseq.txt"
rna = read_tsv(rnafp)
names(rna)[8:10] = c("MDAMB231_R1","MDAMB231_R2","MDAMB231_R3")
rna <- rna %>%
  mutate(ensid = sapply(strsplit(Ensembl_ID,"[.]"),"[[",1))

# merge with annotation
genes_fp <- "/kyber/Data/Nanopore/projects/nanonome/analysis/data/hg38/hg38_transcripts.bed"
genes_all <- read_tsv(genes_fp, col_names = c("chrom","start","end","txid","score","strand","ensid","hgnc")) %>%
  mutate(start = start + 1) # convert it from bed to normal format
genes <- genes_all %>%
  arrange(chrom,start,end) %>%
  distinct(ensid, .keep_all =T)
genes$ensid <- sapply(strsplit(genes$ensid,"[.]"),"[[",1)

# merge the data
rna.tb <- genes %>%
  bind_cols(rna[match(genes$ensid,rna$ensid),] %>% 
    dplyr::select(-ensid,Ensembl_ID)) %>%
  na.omit()

# select ERBB2, ESR1, PGR
select.tb <- rna.tb %>%
  filter(hgnc %in% c("ERBB2","ESR1","PGR"))

# DMRs
dmrs.fp <- file.path(outdir,"200802_bcan_DMRs.tsv")
dars.fp <- file.path(outdir,"200802_bcan_DARs.tsv")
dmrs <- read_tsv(dmrs.fp)
dars <- read_tsv(dars.fp)

dmrs.gr <- GRanges(dmrs)
dars.gr <- GRanges(dars)
select.gr <- GRanges(select.tb)
select.prom <- promoters(select.gr,upstream = 5e3,downstream = 5e3)

another <- GRanges(tibble(chrom = "chr6",start =151805000, end = start + 15000))
dars[subjectHits(findOverlaps(another,dars.gr)),] %>%as.data.frame()

dmr.ovl <- findOverlaps(select.prom,dmrs.gr)
dar.ovl <- findOverlaps(select.prom,dars.gr)

dmrs.sel <- dmrs[subjectHits(dmr.ovl),] %>%
  dplyr::select(Chromosome,Start,End,Comparison, Mean_Difference) %>%
  mutate(hgnc = select.tb$hgnc[queryHits(dmr.ovl)],
    Mean_Difference = Mean_Difference
  ) 

dars.sel <- dars[subjectHits(dar.ovl),] %>%
  dplyr::select(Chromosome,Start,End,Comparison, Mean_Difference) %>%
  mutate(hgnc = select.tb$hgnc[queryHits(dar.ovl)],
    Mean_Difference = Mean_Difference
  ) 

diff.sel <- bind_rows(list(
    DMR = dmrs.sel,
    DAR = dars.sel),.id = "what")
    
# wrtiting tables out
out_path <- file.path(outdir,"bcan_ER_PGR_ESR_expression.tsv")
write_tsv(select.tb,out_path)
out_path <- file.path(outdir,"bcan_ER_PGR_ESR_diffregions.tsv")
write_tsv(diff.sel,out_path)

# plots
sel.gather <- select.tb %>%
  gather(what,count,-chrom,-start,-end,-ensid,-score,-strand,-hgnc,-Ensembl_ID) %>%
  mutate(sample = sapply(strsplit(what,"_"),"[[",1))
plotpath <- file.path(outdir,"200804_ER_PGR_ESR_expression.pdf")
pdf(plotpath)
ggplot(sel.gather,aes( x = sample, y = count)) +
  facet_wrap(~hgnc,scales = "free") +
  geom_point()
dev.off()
