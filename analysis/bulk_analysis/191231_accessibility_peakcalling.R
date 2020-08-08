rm(list=ls(all=TRUE)); gc()
source("~/Code/ilee/plot/ilee_plot_utils.R")
source("~/Code/nanopore-methylation-utilities/methylation_R_utils.R")
library(bsseq)
library(UpSetR)

## prep ----
mfreqdir <- "/uru/Data/Nanopore/projects/nanonome/pooled/mfreq"
samp <- "GM12878_nanoNOMe"
outdir <- "/uru/Data/Nanopore/projects/nanonome/peaks"
plotdir <- "/home/isac/Dropbox/Data/nome-seq/version3_guppy3/plots/peak_calling"
#cpg_fp <- paste0(mfreqdir,"/",samp,".cpg.BSmooth.rds")
gpc_fp <- paste0(mfreqdir,"/",samp,".gpc.BSmooth.rds")

# read data
gpc <- readRDS(gpc_fp)

## filter by coverage ----
gpc.cov <- getCoverage(gpc,type="Cov",what="perBase")[,1]
keepi <- which(gpc.cov>5)
gpc.keep <- gpc[keepi,]

## from here will be fxn
##qcutoff <- 0.99
#gpc.meth <- getMeth(gpc.keep,type="smooth",what="perBase")[,1]
## remove NA
#keepi <- !is.na(gpc.meth)
#gpc.keep <- gpc.keep[keepi,]
#gpc.meth <- gpc.meth[keepi]
#gpc.cov <- getCoverage(gpc.keep,type="Cov",what="perBase")[,1]
#gpc.m <- getCoverage(gpc.keep,type="M",what="perBase")[,1]
### compare to baseline ----
#baseline <- median(gpc.meth,na.rm = T)
#gpc.diff <- gpc.meth - baseline
#cutoff <- 0 # quantile(gpc.diff,qcutoff,na.rm = T)
#gpc.direction <- ifelse(gpc.diff > cutoff, 1, -1) # cut off by qcutoff
#gpc.gr <- granges(gpc.keep)
#chrs <- as.character(seqnames(gpc.gr))
#pos <- start(gpc.gr)
### find regions of + ----
#regions <- bsseq:::regionFinder3(gpc.direction,chrs,pos)$up
#regions <- as_tibble(regions)
### then add average and peak accesibility, along with coverage, then binomial test ---
#regions <- regions %>%
#  rowwise() %>%
#  mutate(
#    coverage = sum(gpc.cov[idxStart:idxEnd]),
#    methylated = sum(gpc.m[idxStart:idxEnd]),
#    average = mean(gpc.meth[idxStart:idxEnd]),
#    peak = max(gpc.meth[idxStart:idxEnd]),
#    p.value = binom.test(methylated,coverage,baseline,alternative = "greater")$p.value
#  ) %>%
#  ungroup()
### multiple test adjusting using FDR
#regions$adjusted.pval  <- p.adjust(regions$p.value,"BH")
#
### significance test :
### 1. peak size,
### 2. width
### 3. alpha
#winsize <- 50
#a <- 0.01
#minpeak <- baseline + quantile(gpc.diff,0.99,na.rm = T)
#
#regions <- regions %>%
#  mutate(
#    width = end - start + 1,
#    sig = ifelse(
#      adjusted.pval <= a &
#        width >= winsize &
#        peak >= minpeak,"sig","insig"))
#table(regions$sig)

# use function
source("~/Code/nanopore-methylation-utilities/methylation_R_utils.R")
regions <- gpcPeakCaller(gpc.keep)
regions %>%
  filter(sig == "sig", peak >= 0.7) %>%
  summarize( min(average), min(peak), n())
table(regions$sig)



# start from here if redoing plots
## write out the regions as rds ---
rdspath <- file.path(outdir,"GM12878_nanoNOMe_accessibility_peaks.rds")
if (0) {
  saveRDS(regions,rdspath)
}

regions <- readRDS(rdspath)
regions.sig <- regions %>%
  filter(sig == "sig")
if (0) {
  # output to tsv
  out <- regions.sig %>%
    dplyr::select(-cluster,-sig) %>%
    dplyr::rename(Chromosome = "chr",
      Start = "start",
      End = "end",
      Num_Sites = "n",
      Observations = "coverage",
      Methylated = "methylated",
      Average_Accessibility = "average",
      Maximum_Accessibility = "peak",
      Region_Width = "width"
    )
  outpath <- file.path(plotdir,"GM12878_nanoNOMe_accessibility_peaks.tsv")
  write_tsv(out,outpath)
}
