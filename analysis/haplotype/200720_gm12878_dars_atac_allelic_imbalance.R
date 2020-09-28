source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
dir <- "/uru/Data/Nanopore/projects/nanonome/haplotype"
hap1cov.fp <- file.path(dir,"GM12878_ATAC_SRR891268_hap_dar_with_het_snp.hap1.cov.bed")
hap2cov.fp <- file.path(dir,"GM12878_ATAC_SRR891268_hap_dar_with_het_snp.hap2.cov.bed")
dars.fp <- "/home/isac/Dropbox/Data/nome-seq/version3_guppy3/plots/haplotypes/200425_gm12878_haplotype_dars_sig.tsv"

# read data
cnames <- c("chr","start","end","cov")
hap1cov <- read_tsv(hap1cov.fp,col_names = cnames)
hap2cov <- read_tsv(hap2cov.fp,col_names = cnames)
dars <- read_tsv(dars.fp)

# summary atac-seq imbalance %
atac <- bind_rows(list(hap1 = hap1cov,hap2 = hap2cov),.id = "hap") %>%
  distinct() %>%
  spread(hap,cov) %>%
  mutate(lab = paste(chr,start,end,sep="_"))

# subset matching dars
dars.sel <- dars %>%
  dplyr::select(chr = Chromosome,start = Start, end = End, Mean_Difference) %>%
  mutate(lab = paste(chr,start,end,sep="_")) %>%
  filter(lab %in% atac$lab)

# summary both
atac.sum <- atac %>%
  mutate(total = hap1 + hap2,
    comp = log((hap1+1)/(hap2+1))) %>%
  filter(total >= 10) %>% # at least 10 reads total
  dplyr::select(lab,comp)

dars.sum <- dars %>%
  dplyr::select(chr = Chromosome,start = Start, end = End, comp = Mean_Difference) %>%
  mutate(lab = paste(chr,start,end,sep="_")) %>%
  filter(lab %in% atac.sum$lab) %>%
  dplyr::select(lab,comp)

sum.tb <- bind_rows(list(`ATAC-seq`= atac.sum,nanoNOMe = dars.sum),.id = "method") %>%
  spread(method,comp)

cornum <- round(cor(sum.tb$`ATAC-seq`,sum.tb$nanoNOMe),2)

plotpath <- "/home/isac/Dropbox/Data/nome-seq/version3_guppy3/plots/haplotypes/200720_gm_dars_atac_comp.pdf"
pdf(plotpath, useDingbats = F, width = 4, height = 4)
ggplot(sum.tb,aes(x = `ATAC-seq`,y = nanoNOMe)) + 
  geom_point() +
  geom_text(label = paste("r =",cornum), x = -4, y = 0.4,hjust = 0) +
  geom_text(label = paste("Number of DARs =",nrow(sum.tb)), x = -4, y = 0.45,hjust = 0) +
  labs( x = "ATAC-seq Coverage log Fold Change\n[log(Paternal + 1 / Maternal + 1)]",
    y = "nanoNOMe Accessibility Difference\n[Paternal - Maternal]")
dev.off()
