#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(GenomicRanges)
library(VariantAnnotation)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

cores=detectCores()
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/mfreq_all")
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
svpath=file.path(root,"pooled/sv/merged_SURVIVOR_1kbp_typesafe_10_31_2018.sort.anno.vcf.gz")
outpath=file.path(root,"sv/svsummary.txt")

vcf = readVcf(svpath,"hg38")
vcf.tb = read_tsv(svpath,comment="##")

vcf.del = vcf[info(vcf)$SVTYPE=="DEL"]

# all possible combinations of supp vec
require(gtools)
perm = as.tibble(permutations(2,3,c(0,1),repeats.allowed=T)) %>%
    transmute(perm=paste0(V1,V2,V3))
perm=perm$perm
supall = tibble(front="111",back=perm) %>%
    transmute(vec=paste0(front,back))
sup10 = tibble(front="10",back=perm) %>%
    transmute(vec=paste0(front,back))
sup7 = tibble(front="100",back=perm) %>%
    transmute(vec=paste0(front,back))
sup231 = tibble(front="1",back=perm) %>%
    transmute(vec=paste0(front,back))
sup107 = tibble(front="110",back=perm) %>%
    transmute(vec=paste0(front,back))
sup10231 = tibble(front="11",back=perm) %>%
    transmute(vec=paste0(front,back))
sup7231 = tibble(front="101",back=perm) %>%
    transmute(vec=paste0(front,back))
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
widths = info(vcf.sub)$AVGLEN
vcf.filt = vcf.sub[which(widths>lenthr)]
ins.counts = as.tibble(
    table(
        supvec.tb$key[match(info(vcf.filt)$SUPP_VEC,
                            supvec.tb$value)]))
# DUP
vcf.sub = vcf[svtypes=="DUP"]
widths = info(vcf.sub)$AVGLEN
vcf.filt = vcf.sub[which(widths>lenthr)]
dup.counts = as.tibble(
    table(
        supvec.tb$key[match(info(vcf.filt)$SUPP_VEC,
                            supvec.tb$value)]))

# INV
vcf.sub = vcf[svtypes=="INV"]
widths = info(vcf.sub)$AVGLEN
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
