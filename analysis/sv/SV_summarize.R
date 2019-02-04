#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(GenomicRanges)
library(VariantAnnotation)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

cores=detectCores()
root="/kyber/Data/Nanopore/projects/nanonome/analysis"
datroot=file.path(root,"data/nanonome/pooled/mfreq"
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
svpath=file.path(root,"data/nanonome/pooled/sniffles/merged_cancer_1kbpdist.vcf.gz")
outpath=file.path(root,"data/nanonome/pooled/sniffles/svsummary.txt")

vcf = readVcf(svpath,"hg38")
vcf.tb = read_tsv(svpath,comment="##")

vcf.del = vcf[info(vcf)$SVTYPE=="DEL"]

# all possible combinations of supp vec
supall = tibble(vec="111")
sup10 = tibble(vec="100")
sup7 = tibble(vec="010")
sup231 = tibble(vec="001")
sup107 = tibble(vec="110")
sup10231 = tibble(vec="101")
sup7231 = tibble(vec="011")
supvec.tb = tibble(all=supall$vec,
                   MCF10A=sup10$vec,
                   MCF7=sup7$vec,
                   MDA231=sup231$vec,
                   MCF10A_MCF7=sup107$vec,
                   MCF10A_MDA231=sup10231$vec,
                   MCF7_MDA231=sup7231$vec) %>%
    type_convert() %>%
    gather()

lenthr = 1000 # length threshold to throw out short (possibly fp) svs

# what are sv types?
svtypes = info(vcf)$SVTYPE

# start with del
vcf.sub = vcf[svtypes=="DEL"]
coords = as.tibble(rowRanges(vcf.sub))
ends = info(vcf.sub)$END
widths = width(GRanges(
    seqnames=as.character(coords$seqnames),
    IRanges(start=coords$start,
            end=ends)))
vcf.filt = vcf.sub[which(widths>lenthr)]
del.counts = as.tibble(
    table(
        supvec.tb$key[match(info(vcf.filt)$SUPP_VEC,
                            supvec.tb$value)]))

# TRA
vcf.sub = vcf[svtypes=="TRA"]
tra.counts = as.tibble(
    table(
        supvec.tb$key[match(info(vcf.sub)$SUPP_VEC,
                            supvec.tb$value)]))
# INS
vcf.sub = vcf[svtypes=="INS"]
widths = info(vcf.sub)$SVLEN
vcf.filt = vcf.sub[which(widths>lenthr)]
ins.counts = as.tibble(
    table(
        supvec.tb$key[match(info(vcf.filt)$SUPP_VEC,
                            supvec.tb$value)]))
# DUP
vcf.sub = vcf[svtypes=="DUP"]
widths = info(vcf.sub)$SVLEN
vcf.filt = vcf.sub[which(widths>lenthr)]
dup.counts = as.tibble(
    table(
        supvec.tb$key[match(info(vcf.filt)$SUPP_VEC,
                            supvec.tb$value)]))

# INV
vcf.sub = vcf[svtypes=="INV"]
widths = info(vcf.sub)$SVLEN
vcf.filt = vcf.sub[which(widths>lenthr)]
inv.counts = as.tibble(
    table(
        supvec.tb$key[match(info(vcf.filt)$SUPP_VEC,
                            supvec.tb$value)]))

# summarize
sum.tb = tibble(cell=inv.counts$Var1,
                del=del.counts$n,
                tra=tra.counts$n,
                dup=dup.counts$n,
                inv=inv.counts$n,
                ins=ins.counts$n) %>%
    mutate(total=del+tra+dup+inv+ins)

write_tsv(sum.tb,outpath)
