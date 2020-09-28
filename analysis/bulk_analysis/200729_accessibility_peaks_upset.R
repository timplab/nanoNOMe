rm(list=ls(all=TRUE)); gc()
source("~/Code/ilee/plot/ilee_plot_utils.R")
source("~/Code/nanopore-methylation-utilities/methylation_R_utils.R")
library(bsseq)
library(UpSetR)

## prep ----
mfreqdir <- "/uru/Data/Nanopore/projects/nanonome/pooled/mfreq"
plotdir <- "/home/isac/Dropbox/Data/nome-seq/version3_guppy3/plots/peak_calling"
samp <- "GM12878_nanoNOMe"
outdir <- "/uru/Data/Nanopore/projects/nanonome/peaks"
regpath <- file.path(plotdir,"GM12878_nanoNOMe_accessibility_peaks.tsv")

regions <- read_tsv(regpath)

## overlap with atac-seq peaks ----
atacpath <- "/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/atac/GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.hg38.bed"
atac <- read_tsv(atacpath,col_names = c("chrom","start","end","id","score","strand"))
dnasepath <- "/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/dnase/ENCFF598KWZ.bed.gz"
dnase <- read_tsv(dnasepath,col_names = c("chrom","start","end","id","score","strand","signal","pval","qval","etc"))

# let's first remove chroms that don't exist in all three
chroms.keep <- as.character(unique(regions$Chromosome))

regions.gr <- GRanges(regions %>% filter( Chromosome %in% chroms.keep))
atac.gr <- GRanges(atac %>% filter(chrom %in% chroms.keep))
dnase.gr <- GRanges(dnase %>% filter(chrom %in% chroms.keep))

# overlaps
# first let's compare nanonome to the other two
nome_atac <- overlapsAny(regions.gr,atac.gr)
nome_dnase <- overlapsAny(regions.gr,dnase.gr)
nome_only <- !(nome_atac | nome_dnase)
# atac as reference
atac_dnase <- overlapsAny(atac.gr,dnase.gr)
atac_nome <- overlapsAny(atac.gr,regions.gr)
atac_only <-  !(atac_dnase | atac_nome)
# dnase as reference
dnase_atac <- overlapsAny(dnase.gr,atac.gr)
dnase_nome <- overlapsAny(dnase.gr,regions.gr)
dnase_only <-  !(dnase_atac | dnase_nome)
# all three
nome_both <- nome_atac & nome_dnase
atac_both <- atac_nome & atac_dnase
dnase_both <- dnase_nome & dnase_atac
# numbers in reference to nanonome
nome_only.num <- sum(nome_only)
atac_only.num <- sum(atac_only)
dnase_only.num <- sum(dnase_only)
all.num <- max(sum(nome_both),sum(atac_both),sum(dnase_both))
nome_atac.num <- max(sum(nome_atac),sum(atac_nome)) - all.num#min(sum(nome_atac),sum(atac_nome)) - all.num
nome_dnase.num <- max(sum(nome_dnase),sum(dnase_nome)) - all.num
atac_dnase.num <- min(sum(atac_dnase),sum(dnase_atac)) - all.num
nome_only.num + all.num + nome_atac.num + nome_dnase.num == length(regions.gr)

# let's plot upset
# order : atac, dnase, nanoNOMe
atac_only.mat <- t(replicate(atac_only.num,c(1,0,0)))
dnase_only.mat <- t(replicate(dnase_only.num,c(0,1,0)))
nome_only.mat <- t(replicate(nome_only.num,c(0,0,1)))
atac_dnase.mat <- t(replicate(atac_dnase.num,c(1,1,0)))
nome_atac.mat <- t(replicate(nome_atac.num,c(1,0,1)))
nome_dnase.mat <- t(replicate(nome_dnase.num,c(0,1,1)))
all.mat <- t(replicate(all.num,c(1,1,1)))
upset.mat <- as.data.frame(rbind(atac_only.mat, dnase_only.mat, nome_only.mat, 
  atac_dnase.mat, nome_atac.mat, nome_dnase.mat, all.mat))
colnames(upset.mat) <- c("ATAC-seq","DNAse-seq","nanoNOMe")

plotdir <- "/home/isac/Dropbox/Data/nome-seq/version3_guppy3/plots/peak_calling"
plotpath <- file.path(plotdir,"GM_peakcalling_upset.pdf")
pdf(plotpath,width = 3, height = 2.5,useDingbats = F)
upset(upset.mat, order.by = "freq",matrix.color = "black",set_size.angles = 30, main.bar.color = "black", mb.ratio = c(0.6,0.4))
dev.off()
