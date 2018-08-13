library(bsseq)
library(tidyverse)
library(cowplot)

outdir="~/Dropbox/timplab_data/nanonome/180809_gm12878"

aligndir="/mithril/Data/NGS/Aligned/160314_na12878/bismark"

raw=read.bismark(files=file.path(aligndir, "rep1.bismark.cov.gz"), sampleNames="GM12878_rep1", strandCollapse=F)

meth=getMeth(raw, type="raw", what="perBase")

locs=granges(raw)

meth.plot=tibble(chr=as.character(seqnames(locs)), meth=as.numeric(meth))


meth.plot=meth.plot %>%
    filter(chr %in% c("chrX", "chr22", "chr21"))

pdf(file.path(outdir, "test.pdf"))

ggplot(meth.plot, aes(x=meth, color=chr, group=chr))+geom_density()
ggplot(meth.plot, aes(x=meth, color=chr, group=chr))+geom_freqpoly()
ggplot(meth.plot, aes(x=meth, color=chr, group=chr))+stat_ecdf()

dev.off()
