#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(GenomicRanges)
library(VariantAnnotation)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

cores=detectCores()
parser="../../script/SVmethylation.py"
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/sv")
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
pd = tibble(cell=cells,
            cpgpath=file.path(datroot,
                              paste(cell,"cpg.pooled.meth.bed.gz",sep=".")),
            gpcpath=file.path(datroot,
                              paste(cell,"gpc.pooled.meth.bed.gz",sep=".")),
            bampath=file.path(root,"pooled/bam",paste0(cell,".pooled.bam")))
pd$svpath = sapply(cells,function(x){
    system(paste0("readlink -f ",root,"/pooled/sv/",x,".pooled*vcf.gz"),T)})

svpath=file.path(root,"pooled/sv/merged_SURVIVOR_1kbp_typesafe_10_31_2018.sort.anno.vcf.gz")


# fix contigs if necessary
fixcontigs <- function(granges){
    contigs=seqlevels(granges)
    if (! grepl("chr",contigs[1])){
        contigs.fix = paste0('chr',contigs)
        contigs.fix[match("chrMT",contigs.fix)] = "chrM"
        seqlevels(granges)=contigs.fix
    }
    granges
}
getEndCoord <- function(vcf){
    GRanges(info(vcf)$CHR2,
            IRanges(start=info(vcf)$END,end=info(vcf)$END))
}
genocoTogranges <- function(coords){
    coords.list = strsplit(coords,"-")
    bp1 = sapply(coords.list,"[[",1)
    bp1.tb = as.tibble(do.call(rbind,strsplit(bp1,"_"))) %>% type_convert()
    bp1.gr = fixcontigs(GRanges(bp1.tb$V1,IRanges(bp1.tb$V2,bp1.tb$V2)))
    bp2 = sapply(coords.list,"[[",2)
    bp2.tb = as.tibble(do.call(rbind,strsplit(bp2,"_"))) %>% type_convert()
    bp2.gr = fixcontigs(GRanges(bp2.tb$V1,IRanges(bp2.tb$V2,bp2.tb$V2)))
    bp1.gr$bp=1;bp2.gr$bp=2
    list(bp1.gr,bp2.gr)
}
makeGenoCo <- function(vcf){
    bp1 = rowRanges(vcf)
    bp2 = getEndCoord(vcf)
    paste(sep="-",
          paste(seqnames(bp1),start(bp1),sep="_"),
          paste(seqnames(bp2),start(bp2),sep="_"))
}

    

vcf = readVcf(svpath,"hg38")

vcf.tra = vcf[info(vcf)$SVTYPE=="TRA"]
combos = tibble(one=c("MCF10A","MCF10A"),
                two=c("MCF7","MDAMB231"),
                vec=c(100000,001000),
                ind=c(1,3))
window=2000
    
i = 1
for (i in seq(dim(combos)[1])){
    oneidx = match(combos$one[i],pd$cell)
    twoidx = match(combos$two[i],pd$cell)
    
    genes.fp = file.path(root,"annotations/breastcancer",paste0("MCF10A_vs_",combos$two[i],"_genes.bed"))
    genes = read_tsv(genes.fp) %>%
        dplyr::rename(chrom=`#chrom`)
    genes.gr = GRanges(genes)
    
    vcf.sub = vcf[info(vcf)$SUPP_VEC==combos$vec[i]]
    bpall.list = genocoTogranges(geno(vcf.sub)$CO[,combos$ind[i]])
    ovl.list = lapply(bpall.list,function(x){
        dist = as.data.frame(distanceToNearest(x,genes.gr))
        dist[which(dist$distance<10000),]
    })
    ovl = do.call(rbind,ovl.list)
    ovlidx = unique(ovl$queryHits)

    vcf.gene = vcf.sub[ovlidx]
    
    coords = geno(vcf.gene)$CO[,combos$ind[i]]
    bp.list = genocoTogranges(coords)

    bprange.list = lapply(bp.list,function(x){
        resize(x,window,fix="center")})
    bpranges = do.call(c,bprange.list)
    
    # read in reads
    cpg.list = lapply(seq_along(bprange.list),function(i){
        x=bprange.list[[i]]
        cpg.one = tabix_mbed(pd$cpgpath[oneidx],x,by="read")%>%
            mutate(bp=i,comp=1,calltype="cpg")
        cpg.two = tabix_mbed(pd$cpgpath[twoidx],x,by="read")%>%
            mutate(bp=i,comp=2,calltype="cpg")
        bind_rows(cpg.one,cpg.two)
    })
    reads.cpg = do.call(rbind,cpg.list)
    gpc.list = lapply(seq_along(bprange.list),function(i){
        x=bprange.list[[i]]
        gpc.one = tabix_mbed(pd$gpcpath[oneidx],x,by="read")%>%
            mutate(bp=i,comp=1,calltype="gpc")
        gpc.two = tabix_mbed(pd$gpcpath[twoidx],x,by="read")%>%
            mutate(bp=i,comp=2,calltype="gpc")
        bind_rows(gpc.one,gpc.two)
    })
    reads.gpc = do.call(rbind,gpc.list)
    reads = bind_rows(reads.cpg,reads.gpc)
    mods=c("cpg","gpc")
    reads.list = list(reads.cpg,reads.gpc)
    
    # compare SV in two and nonSV in one
    sv.reads = lapply(reads.list,function(x){
        x[which(x$comp == 2),]})
    nonsv.reads = lapply(reads.list,function(x){
        x[which(x$comp == 1),]})
    sv.calls = mclapply(mc.cores=2,seq_along(sv.reads),function(i){
        mbedByCall(sv.reads[[i]]) %>%
        mutate(type="SV",calltype=mods[i])})
    nonsv.calls = mclapply(mc.cores=2,seq_along(nonsv.reads),function(i){
        mbedByCall(nonsv.reads[[i]]) %>%
        mutate(type="nonSV",calltype=mods[i])})
    calls = bind_rows(do.call(rbind,sv.calls),
                      do.call(rbind,nonsv.calls))
    calls.bp = calls[overlapsAny(GRanges(calls),bpranges),]

    freq = calls.bp %>% na.omit() %>%
        group_by(calltype,type,chrom,start) %>%
        summarize(cov=n(),meth=sum(mcall),freq=meth/cov)%>%
        mutate(end=start)

    plotpath = file.path(plotdir,paste0(combos$one[i],"_vs_",combos$two[i],"_allSV_near_genes_byregion.pdf"))
    pdf(plotpath,6,4,useDingbats=F)
    for (j in seq_along(vcf.gene)){
        print(j)
        bp1.range = bprange.list[[1]][j]
        bp2.range = bprange.list[[2]][j]
        bp1 = bp.list[[1]][j]
        bp2 = bp.list[[2]][j]
        freq.bp1 = freq[overlapsAny(GRanges(freq),bp1.range),] %>%
            mutate(bp = start(bp1),bpidx=1)
        freq.bp2 = freq[overlapsAny(GRanges(freq),bp2.range),] %>%
            mutate(bp = start(bp2),bpidx=2)
        freq.bp = bind_rows(freq.bp1,freq.bp2) %>%
            filter(cov>5)
        g = ggplot(freq.bp,aes(x=start,y=freq,color=type,group=type))+
            facet_grid(calltype~bpidx,scales="free")+
            lims(y=c(0,1))+
            geom_smooth(se=F)+
            geom_point(alpha=0.5)+
            geom_vline(aes(xintercept=bp))+
            geom_rug(sides="b",alpha=0.5)+
            theme_bw()+
            theme(legend.position="bottom",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            ggtitle(paste(info(vcf.gene[i])$SVTYPE,
                          "; bp1 :",toString(bp1),"; bp2 :",toString(bp2)))
        gcov = ggplot(freq.bp,aes(x=start,y=cov,color=type,group=type))+
            facet_grid(calltype~bpidx,scales="free")+
                        geom_smooth(se=F)+
            geom_point(alpha=0.5)+
            geom_vline(aes(xintercept=bp))+
            geom_rug(sides="b",alpha=0.5)+
            theme_bw()+
            theme(legend.position="bottom",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            ggtitle(paste(info(vcf.gene[i])$SVTYPE,
                          "; bp1 :",toString(bp1),
                          "; bp2 :",toString(bp2)))
        print(g)
    #    print(gcov)
    }
    dev.off()

}
