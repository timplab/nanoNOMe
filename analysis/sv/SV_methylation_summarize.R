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

combos = tibble(one=c("MCF10A","MCF10A"),
                two=c("MCF7","MDAMB231"),
                ind=c(10,12))
writebed = TRUE
makeplots = TRUE
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

cpgwindow = 400
gpcwindow = 50
# deletions
plotpath = file.path(plotdir,"bcan_hom_DEL_bp1_vs_bp2_comparison.pdf")
if (writebed){file.remove(bedout)}
if (makeplots){
    pdf(plotpath,useDingbats=F,width=3,height=3)
}
i = 1
for (i in seq(dim(combos)[1])){
    print(i)
    del.idx = which((info(vcf)$SUPP_VEC %in% sup.list[[i]])&
                    (info(vcf)$SVTYPE == "DEL")) # only dels with these SUPP vecs
    vcf.del = vcf[del.idx]
    vcftb.del = vcf.tb[del.idx,]
    # filtering for hom dels only
    geno = sapply(strsplit(unlist(vcftb.del[,combos$ind[i]]),":"),"[[",1)
    hom.idx = which(geno=="1/1")
    vcf.del=vcf.del[hom.idx]
    vcftb.del = vcftb.del[hom.idx,]
    # filtering out <1kb dels
    widths = geno(vcf.del)[["LN"]][,combos$ind[i]-9]
    longidx = which(widths>1000)
    vcf.del = vcf.del[longidx]
    vcftb.del = vcftb.del[longidx,]

    # get coords of dels
    endcoord = as.tibble(info(vcf.del)[,c("CHR2","END")])
    coord.tb = bind_cols(as.tibble(rowRanges(vcf.del))[,c("seqnames","start")],
                         endcoord) %>%
        filter(seqnames==CHR2) %>%
        dplyr::select(seqnames,start,end=END)
    del.gr = GRanges(coord.tb)
    del.gr = fixcontigs(del.gr)
    bp1.gr = flank(del.gr,width=1) # bp1 5'
    bp2.gr = flank(del.gr,start=F,width=1) # bp2 3'
    bp.list = list(bp1.gr,bp2.gr)
    # different windows for cpg and gpc 
    reg1.cpg = resize(bp1.gr,cpgwindow,fix="end")
    reg2.cpg = resize(bp2.gr,cpgwindow,fix="start")
    reg1.gpc = resize(bp1.gr,gpcwindow,fix="end")
    reg2.gpc = resize(bp2.gr,gpcwindow,fix="start")
    reg.cpg.list = list(reg1.cpg,reg2.cpg)
    reg.gpc.list = list(reg1.gpc,reg2.gpc)
    regs.list.list = list(reg.cpg.list,reg.gpc.list)
    # read in methylation
    pd.sub = pd[which(pd$cell %in% c(combos$one[i],combos$two[i])),] %>%
        mutate(cellidx=ifelse(cell==combos$one[i],"one","two"))
    pd.sub = bind_rows(pd.sub%>%mutate(bp=1),
                       pd.sub%>%mutate(bp=2)) %>%
        mutate(calltypeidx = ifelse(calltype=="cpg",1,2))
    pd.cpg = pd.sub[which(pd.sub$calltype=="cpg"),]
    pd.gpc = pd.sub[which(pd.sub$calltype=="gpc"),]
    meth.list = mclapply(mc.cores=8,seq(dim(pd.sub)[1]),function(i){
        regs = regs.list.list[[pd.sub$calltypeidx[i]]][[pd.sub$bp[i]]]
        tabix_mfreq(pd.sub$filepath[i],regs) %>%
            getRegionMeth(regs)%>%
            mutate(cell=pd.sub$cell[i],
                   calltype=pd.sub$calltype[i],
                   cellidx=pd.sub$cellidx[i],
                   bp=pd.sub$bp[i])
    })
    # filter out anythign that has tot cov > 1000 -
    # these are probably low mappability regions with secondary alignments
    meth = do.call(rbind,meth.list) %>%
        filter(totcov<=1000 & totcov >= 5)
    meth.spread = meth %>%
        dplyr::select(-totcov,-numsites)%>%
        spread(bp,freq) %>% na.omit() %>%
        mutate(lab="sample")
    ids = unique(meth.spread$feature.index)
    randidx = sample(ids)
    meth.spread$randomidx = randidx[match(meth.spread$feature.index,ids)]
    # deldel 
    

    
    # boxplot of  vs random
    box = bind_rows(meth.spread,meth.random) %>%
        mutate(del = `2`-`1`,
               or = paste(cell,calltype))

    g.box = ggplot(box,aes(x=or,y=del,group=or,color=or)) +
        facet_wrap(~lab)+
        geom_boxplot()
    print(g.box)
    makeplot <- function(x){
        maxval=max(c(x$`1`,x$`2`))
        g = ggplot(x,aes(x=`1`,y=`2`))+
            lims(x=c(0,maxval),y=c(0,maxval))+
            facet_wrap(~cell,nrow=1)+
            geom_point()+
            theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            labs(x="BP1",y="BP2")+
            coord_fixed()
        g
    }
    g.cpg = makeplot(meth.spread[which(meth.spread$calltype=="cpg"),])+
        ggtitle(paste("CpG :",combos$one[i],"vs",combos$two[i]))
    g.gpc = makeplot(meth.spread[which(meth.spread$calltype=="gpc"),])+
        ggtitle(paste("GpC :",combos$one[i],"vs",combos$two[i]))
    print(g.cpg)
    print(g.gpc)
    # dels
    delplot <- function(x){
        maxval = max(c(x$one,x$two))
        g = ggplot(x,aes(x=one,y=two)) +
            lims(x=c(-maxval,maxval),y=c(-maxval,maxval))+
            geom_point() +
            theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            coord_fixed()
    }
    meth.del = meth.spread %>%
        mutate(del=`2`-`1`) %>%
        dplyr::select(-`1`,-`2`,-cell) %>%
        spread(cellidx,del) %>%na.omit()
    g.cpg = delplot(meth.del[which(meth.del$calltype=="cpg"),])+
        labs(x=combos$one[i],y=combos$two[i])+
        ggtitle(paste("CpG delta :","BP2 - BP1"))
    g.gpc = delplot(meth.del[which(meth.del$calltype=="gpc"),])+
        labs(x=combos$one[i],y=combos$two[i])+
        ggtitle(paste("GpC delta :","BP2 - BP1"))
    print(g.cpg)
    print(g.gpc)
}
if (makeplots){ dev.off()}

# ins
plotpath = file.path(plotdir,"bcan_hom_INS_bp1_vs_bp2_comparison.pdf")
if (makeplots){
    pdf(plotpath,useDingbats=F,width=3,height=3)
}
i = 1
for (i in seq(dim(combos)[1])){
    print(i)
    ins.idx = which((info(vcf)$SUPP_VEC %in% sup.list[[i]])&
                    (info(vcf)$SVTYPE == "INS")) # only dels with these SUPP vecs
    vcf.ins = vcf[ins.idx]
    vcftb.ins = vcf.tb[ins.idx,]
    # filtering for hom only
    geno = sapply(strsplit(unlist(vcftb.ins[,combos$ind[i]]),":"),"[[",1)
    hom.idx = which(geno=="1/1")
    vcf.ins=vcf.ins[hom.idx]
    vcftb.ins = vcftb.ins[hom.idx,]
    # filtering out <1kb 
    widths = geno(vcf.ins)[["LN"]][,combos$ind[i]-9]
    longidx = which(widths>200)
    vcf.ins = vcf.ins[longidx]
    vcftb.ins = vcftb.ins[longidx,]

    ins.gr = fixcontigs(rowRanges(vcf.ins))
    # different windows for cpg and gpc 
    reg1.cpg = resize(ins.gr,cpgwindow,fix="end")
    reg2.cpg = resize(ins.gr,cpgwindow,fix="start")
    reg1.gpc = resize(ins.gr,gpcwindow,fix="end")
    reg2.gpc = resize(ins.gr,gpcwindow,fix="start")
    reg.cpg.list = list(reg1.cpg,reg2.cpg)
    reg.gpc.list = list(reg1.gpc,reg2.gpc)
    regs.list.list = list(reg.cpg.list,reg.gpc.list)
    # read in methylation
    pd.sub = pd[which(pd$cell %in% c(combos$one[i],combos$two[i])),] %>%
        mutate(cellidx=ifelse(cell==combos$one[i],"one","two"))
    pd.sub = bind_rows(pd.sub%>%mutate(bp=1),
                       pd.sub%>%mutate(bp=2)) %>%
        mutate(calltypeidx = ifelse(calltype=="cpg",1,2))
    pd.cpg = pd.sub[which(pd.sub$calltype=="cpg"),]
    pd.gpc = pd.sub[which(pd.sub$calltype=="gpc"),]
    meth.list = mclapply(mc.cores=8,seq(dim(pd.sub)[1]),function(i){
        regs = regs.list.list[[pd.sub$calltypeidx[i]]][[pd.sub$bp[i]]]
        tabix_mfreq(pd.sub$filepath[i],regs) %>%
            getRegionMeth(regs)%>%
            mutate(cell=pd.sub$cell[i],
                   calltype=pd.sub$calltype[i],
                   cellidx=pd.sub$cellidx[i],
                   bp=pd.sub$bp[i])
    })
    # filter out anythign that has tot cov > 1000 -
    # these are probably low mappability regions with secondary alignments
    meth = do.call(rbind,meth.list) %>%
        filter(totcov<=1000 & totcov >= 5)
    meth.spread = meth %>%
        dplyr::select(-totcov,-numsites) %>%
        spread(bp,freq) %>% na.omit()
    makeplot = function(x){
        maxval=max(c(x$`1`,x$`2`))
        g = ggplot(x,aes(x=`1`,y=`2`))+
            lims(x=c(0,maxval),y=c(0,maxval))+
            facet_wrap(~cell,nrow=1)+
            geom_point()+
            theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            labs(x="BP1",y="BP2")+
            coord_fixed()
        g
    }
    g.cpg = makeplot(meth.spread[which(meth.spread$calltype=="cpg"),])+
        ggtitle(paste("CpG :",combos$one[i],"vs",combos$two[i]))
    g.gpc = makeplot(meth.spread[which(meth.spread$calltype=="gpc"),])+
        ggtitle(paste("GpC :",combos$one[i],"vs",combos$two[i]))
    print(g.cpg)
    print(g.gpc)
    # dels
    delplot <- function(x){
        maxval = max(c(x$one,x$two))
        g = ggplot(x,aes(x=one,y=two)) +
            lims(x=c(-maxval,maxval),y=c(-maxval,maxval))+
            geom_point() +
            theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            coord_fixed()
    }
    meth.del = meth.spread %>%
        mutate(del=`2`-`1`) %>%
        dplyr::select(-`1`,-`2`,-cell) %>%
        spread(cellidx,del) %>%na.omit()
    g.cpg = delplot(meth.del[which(meth.del$calltype=="cpg"),])+
        labs(x=combos$one[i],y=combos$two[i])+
        ggtitle(paste("CpG delta :","BP2 - BP1"))
    g.gpc = delplot(meth.del[which(meth.del$calltype=="gpc"),])+
        labs(x=combos$one[i],y=combos$two[i])+
        ggtitle(paste("GpC delta :","BP2 - BP1"))
    print(g.cpg)
    print(g.gpc)
}
if (makeplots){ dev.off()}
