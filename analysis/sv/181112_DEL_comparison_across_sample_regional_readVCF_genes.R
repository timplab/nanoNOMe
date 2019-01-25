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
datroot=file.path(root,"pooled/methylation/mfreq_all")
plotdir=file.path(root,"plots/sv")
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
pd = tibble(cell=rep(cells,each=2),
            calltype=rep(c("cpg","gpc"),4),
            filepath=file.path(datroot,
                               paste(cell,calltype,"methfreq.txt.gz",
                                     sep=".")))

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
# DE genes
#genes.fp = file.path(root,"annotations/breastcancer",paste0("MCF10A_vs_",combos$two[i],"_genes.bed"))
#genes = read_tsv(genes.fp) %>%
#    dplyr::rename(chrom=`#chrom`)
#genes.gr = GRanges(genes)
    
vcf = readVcf(svpath,"hg38")
vcf.tb = read_tsv(svpath,comment="##")

vcf.del = vcf[info(vcf)$SVTYPE=="DEL"]
combos = tibble(one=c("MCF10A","MCF10A"),
                two=c("MCF7","MDAMB231"),
                ind=c(10,12))
writebed = TRUE
makeplots = FALSE
# all possible combinations of supp vec
require(gtools)
perm = as.tibble(permutations(2,3,c(0,1),repeats.allowed=T)) %>%
    transmute(perm=paste0(V1,V2,V3))
perm=perm$perm
sup7 = tibble(front=rep(c("100","101"),length(perm)),back=rep(perm,each=2)) %>%
    transmute(vec=paste0(front,back))
sup7=as.numeric(sup7$vec)
sup231 = tibble(front=rep(c("001","101"),length(perm)),back=rep(perm,each=2)) %>%
    transmute(vec=paste0(front,back))
sup231=as.numeric(sup231$vec)

sup.list = list(sup7,sup231)
suppvec = info(vcf)$SUPP_VEC

window=4000

plotpath = file.path(plotdir,"bcan_hom_DEL_bp1_vs_bp2_frequency_delta_comparison.pdf")
plotpath = file.path(plotdir,"bcan_hom_DEL_bp1_vs_bp2_frequency_region_plots.pdf")
bedout = file.path(root,"sv/bcan_hom_DEL.bed")
if (writebed){file.remove(bedout)}
if (makeplots){
    pdf(plotpath,useDingbats=F,width=3,height=5)
}
i = 1
for (i in seq(dim(combos)[1])){
    print(i)
    del.idx = which((info(vcf)$SUPP_VEC %in% sup.list[[i]])&
                    (info(vcf)$SVTYPE == "DEL")) # only dels with these SUPP vecs
    vcf.del = vcf[del.idx]
    vcftb.del = vcf.tb[del.idx,]
    geno = sapply(strsplit(unlist(vcftb.del[,combos$ind[i]]),":"),"[[",1)
    hom.idx = which(geno=="1/1") # filtering for hom dels only
    vcf.del=vcf.del[hom.idx]
    vcftb.del = vcftb.del[hom.idx,]

    # get coords of dels
    endcoord = as.tibble(info(vcf.del)[,c("CHR2","END")])
    coord.tb = bind_cols(as.tibble(rowRanges(vcf.del))[,c("seqnames","start")],
                         endcoord) %>%
        filter(seqnames==CHR2) %>%
        dplyr::select(seqnames,start,end=END)
    del.gr = GRanges(coord.tb)
    del.gr = del.gr[which(width(del.gr)>5000)] # filtering <5kb dels
    del.gr = fixcontigs(del.gr)
    bp1.gr = flank(del.gr,width=1) # bp1 5'
    bp2.gr = flank(del.gr,start=F,width=1) # bp2 3'
    reg1.gr = resize(bp1.gr,window,fix="center")
    reg2.gr = resize(bp2.gr,window,fix="center")
    bp.list = list(bp1.gr,bp2.gr)
    reg.list = list(reg1.gr,reg2.gr)
    regs = c(reg1.gr,reg2.gr)
    if (writebed) {
        bed.tb = as.tibble(del.gr)%>%
            mutate(cell=combos$two[i],score=".") %>%
            filter(width<100000) %>% # filter out > 100kb for visualization
            dplyr::select(seqnames,start,end,cell,score,strand)
        write_tsv(bed.tb,bedout,append=T,col_names=F)
    }
    if (! makeplots){ next }
    # read in methylation
    pd.sub = pd[which(pd$cell %in% c(combos$one[i],combos$two[i])),] %>%
        mutate(cellidx=ifelse(cell==combos$one[i],"one","two"))
    meth.list = mclapply(mc.cores=8,seq(dim(pd.sub)[1]),function(i){
        tabix_mfreq(pd.sub$filepath[i],regs) %>%
            mutate(cell=pd.sub$cell[i],
                   calltype=pd.sub$calltype[i],
                   cellidx=pd.sub$cellidx[i])
    })
    meth = do.call(rbind,meth.list)
    meth.gr = GRanges(meth)
    # plot regions
    plotter <- function(data,bps,what="freq"){
        names(data)[which(names(data)==what)]="val"
        bp.tb = as.tibble(bps) %>% mutate(bp=c(1,2))
        g=ggplot(data,aes(x=start,y=val,color=cell,group=cell))+
            geom_point(size=0.5,alpha=0.5)+
            geom_vline(data=bp.tb,mapping=aes(xintercept=start),linetype=2)+
            facet_wrap(bp~.,scales="free",nrow=2)+
            labs(x=paste0("Coordinate on ",seqnames(bps)[1]),
                 y="Methylation Frequency")+
            theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"),
                  legend.position="bottom")
        if (what=="freq") {
            g=g+lims(y=c(0,1))
        }
        g
    }
    for (j in seq_along(del.gr)){
        print(j)
        regs.sub = c(reg1.gr[j],reg2.gr[j])
        bp.sub = c(bp1.gr[j],bp2.gr[j])
        ovl = findOverlaps(meth.gr,regs.sub)
        meth.ovl = meth[queryHits(ovl),] %>%
            mutate(bp = subjectHits(ovl))
        cov.sub = meth.ovl[,c("start","cov","calltype","cellidx","bp")] %>%
            spread(cellidx,cov) %>%
            mutate(two = ifelse(is.na(two),0,two)) %>%
            na.omit() %>%
            gather(cellidx,cov,-start,-calltype,-bp) %>%
            mutate(cell = ifelse(cellidx=="one",combos$one[i],combos$two[i]))
        print(plotter(cov.sub[which(cov.sub$calltype=="cpg"),],bp.sub,what="cov")+
              geom_line() + labs(title="cpg coverage",y="Coverage"))
        print(plotter(meth.ovl[which(meth.ovl$calltype=="cpg"),],bp.sub)+
              geom_line() + labs(title="cpg meth"))
        print(plotter(cov.sub[which(cov.sub$calltype=="gpc"),],bp.sub,what="cov")+
              geom_line() + labs(title="gpc coverage",y="Coverage"))
        print(plotter(meth.ovl[which(meth.ovl$calltype=="gpc"),],bp.sub)+
              geom_line() + labs(title="gpc meth"))

    }
    
}
if (makeplots){ dev.off()}
if (FALSE) {
    meth.gpc = meth %>% filter(calltype=="gpc")
    deleted.tb = meth.gpc %>% dplyr::select(-meth,-unmeth,-freq,-cell)%>%
        spread(cellidx,cov) %>% mutate(two=ifelse(is.na(two),0,two)) %>%
        na.omit() %>% filter(two==0)
    deleted.gr = GRanges(deleted.tb)
    ovl=findOverlaps(deleted.gr,reg1.gr)
    deleted.tb[queryHits(ovl),] %>%
        mutate(featureidx=subjectHits(ovl))
        
    
    # compare deltas b/w MCF10A (no del) vs sample (del)
    meth.del = meth %>% dplyr::select(-numsites,-totcov) %>%
        spread(bp,freq) %>% na.omit() %>%
        mutate(del=`2`-`1`) %>%
        dplyr::select(-`1`,-`2`) %>%
        spread(cell,del) %>%na.omit()
    names(meth.del) = c(names(meth.del)[1:2],"one","two")

    makeplot = function(x){
        lim = max(abs(c(x$one,x$two)))
        g = ggplot(x,aes(x=one,y=two))+
            lims(x=c(-lim,lim),y=c(-lim,lim))+
            geom_point()+
            theme_bw()+
            theme(legend.position="bottom",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            ggtitle(paste(combos$one[i],"vs",combos$two[i],x$calltype[1]))
        g
    }
        
    # scatter
    g.cpg = makeplot(meth.del[which(meth.del$calltype=="cpg"),])
    g.gpc = makeplot(meth.del[which(meth.del$calltype=="gpc"),])
    print(g.cpg)
    print(g.gpc)
}

