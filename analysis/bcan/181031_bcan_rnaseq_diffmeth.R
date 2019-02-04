#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
library(tidyverse)
library(GenomicRanges)
library(parallel)
source(file.path(srcdir,"../../script/methylation_plot_utils.R"))

# set this to TRUE to remove unnecessary objects throughout the process
limitedmem=TRUE
cores = detectCores()-2

# set directories
root=commandArgs(trailingOnly=TRUE)[1]
if (is.na(root)){
    root="/kyber/Data/Nanopore/projects/nanonome/analysis"
}
datroot=file.path(root,"data/nanonome/pooled/mfreq")
cells=c("MCF10A","MCF7","MDAMB231")
fpaths=tibble(cell=cells,
          cpg=file.path(datroot,paste0(cell,"_nanoNOMe.pooled.cpg.mfreq.txt.gz")),
          gpc=file.path(datroot,paste0(cell,"_nanoNOMe.pooled.gpc.mfreq.txt.gz")))
pd=gather(fpaths,key=calltype,value=filepath,-cell)
plotdir=file.path(root,"plots")
annodir=file.path(root,"data/bcan")

# hg38 gene db
dbpath = file.path(root,"data/hg38/hg38_genes.bed")
db.gr = load_db(dbpath,extracols=c("genename","fxn"))
db.id = sapply(strsplit(db.gr$id,"[.]"),"[[",1)
db.gr$id = db.id
prom.gr = promoters(db.gr,200,100)
# cgi db
cgipath = file.path(root,"data/hg38/hg38_cgi.bed")
cgi.gr = load_db(cgipath)
cgi.large = cgi.gr[which(width(cgi.gr)>1000)]

# rna-seq data
rnafp = file.path(annodir,"bcan_rnaseq.txt")
rna = read_tsv(rnafp)
rna.id = sapply(strsplit(rna$Ensembl_ID,"[.]"),"[[",1)
rna$Ensembl_ID = rna.id
names(rna)[8:10] = c("MDAMB231_R1","MDAMB231_R2","MDAMB231_R3")
samples = c("MCF10A","MCF7","MDAMB231")
combos = as.tibble(t(combn(samples,2)))
names(combos) = c("one","two")
# get t-stat
mindiff = 5
rna.pval = tibble(geneid=character(),one=character(),
                  two=character(),p=numeric(),
                  direction=character())
require(gtools)
for (i in seq(dim(combos)[1])){
    one = combos$one[i]
    two = combos$two[i]
    cat(paste(one,two,"\n"))
    dat.one = rna[,grep(one,names(rna))]
    dat.two = rna[,grep(two,names(rna))]
    for (j in seq(dim(rna)[1])){
        one.j = as.numeric(dat.one[j,])
        two.j = as.numeric(dat.two[j,])
        fold = foldchange(mean(one.j),mean(two.j))
        if (! is.nan(fold) & is.finite(fold) & abs(fold) >= mindiff) {
            pval = try(t.test(log(one.j+1),log(two.j+1))$p.value,TRUE)
            if(is.numeric(pval)){
                print(j)
                direction = ifelse(fold>0,"under","over") # two in comparison to one
                rna.pval = bind_rows(rna.pval,
                                    tibble(geneid = rna$Ensembl_ID[j],
                                           one = one,two = two,
                                           p = pval,direction=direction))
            }
        }
    }
}
# thresholding
a = 0.01
rna.sig = rna.pval[which(rna.pval$p <= a),]
# match names
matchidx = match(rna.sig$geneid,db.gr$id)
rna.sig = bind_cols(rna.sig,as.tibble(db.gr)[matchidx,]) %>%na.omit()
sig.gr = GRanges(rna.sig)
sig.tss = promoters(sig.gr,0,2)
sig.window = resize(sig.tss,2002,fix="center")
sig.prom = promoters(sig.gr,200,200)

# which are under/over expressed in both 7 and 231
common = rna.sig %>%
    filter(one=="MCF10A") %>%
    select(seqnames,start,end,strand,geneid,genename,one,two,direction) %>%
    spread(two,direction) %>% na.omit() %>%
    filter(MCF7==MDAMB231)

# CGI
common.gr = GRanges(common)
common.tss = promoters(common.gr,0,2)
ovl = findOverlaps(common.tss,cgi.large)
common.cgi = common[queryHits(ovl),]
common.cgi.tss = as.tibble(promoters(GRanges(common.cgi),0,2))
if (TRUE) {
    outpath = file.path(annodir,"bcan_diffexp_bothcancer_cgi.TSS.bed")
    bed = common.cgi.tss %>%
        select(seqnames,start,end,genename,MCF7,strand,geneid)
    write_tsv(bed,outpath,col_names=F)
}


meth.list = list()
export=FALSE
cgi=FALSE
plot=TRUE
for (i in seq(dim(combos)[1])){
    onename = combos$one[i]
    twoname = combos$two[i]
    prefix = paste0(onename,"_vs_",twoname,
                    "_by_expression_foldchange_",mindiff,
                    "_alpha_",a)

    pd.sub = pd%>%filter(cell==onename | cell==twoname)
    pd.sub = pd.sub %>%
        mutate(cell=ifelse(cell==onename,"one","two"))
    subidx = which(rna.sig$one==onename &
                   rna.sig$two==twoname)
    prom.sub = sig.prom[subidx]
    window.sub = sig.window[subidx]
    tss.sub = sig.tss[subidx]
    # looking into CGI promoters?
    if (cgi){
        prefix=paste0(prefix,"_cgi")
        ovl.idx = overlapsAny(prom.sub,cgi.gr)
        prom.sub = prom.sub[ovl.idx]
        tss.sub = sig.tss[ovl.idx]
        window.sub = window.sub[ovl.idx]
    }
    # methyltion data
    dat.list = mclapply(mc.cores=cores,seq(dim(pd.sub)[1]),function(i){
        tabix_mfreq(pd.sub$filepath[i],prom.sub)
    })
    regmeth.list = lapply(seq_along(dat.list),function(i){
        getRegionMeth(dat.list[[i]],prom.sub) %>%
            mutate(cell=pd.sub$cell[i],
                   calltype=pd.sub$calltype[i])
    })
    regmeth = do.call(rbind,regmeth.list)
    meth.spread = regmeth %>% select(-totcov,-numsites) %>%
        spread(cell,freq)

    chrom.one = meth.spread%>% select(-two)%>%
        spread(calltype,one) %>% mutate(chrom = sqrt((1-cpg)**2+gpc**2))
    chrom.two = meth.spread%>% select(-one)%>%
        spread(calltype,two) %>% mutate(chrom = sqrt((1-cpg)**2+gpc**2))
    chrom = tibble(feature.index=chrom.two$feature.index,
                   calltype="chrom",
                   one=chrom.one$chrom,two=chrom.two$chrom) %>% na.omit()
    meth.spread = bind_rows(meth.spread,chrom)
    meth.spread = transmute(meth.spread,feature.index=feature.index,
                         del=two-one,calltype=calltype) %>%
        bind_cols(as.data.frame(window.sub[meth.spread$feature.index]))
    meth.plt = meth.spread %>%
        select(-one,-two)%>%
        spread(calltype,del)%>%na.omit() %>%
        mutate(onename=onename,twoname=twoname)
    ymax=as.numeric(quantile(abs(meth.plt$cpg),0.99))
    xmax=as.numeric(quantile(abs(meth.plt$gpc),0.99))
    xthr = xmax/10
    ythr = ymax/10
    meth.list[[i]] = meth.plt %>%
        filter((cpg > ythr & cpg < ymax )|
               (cpg < -ythr & cpg > -ymax )|
               (gpc > xthr & gpc < xmax )|
               (gpc < -xthr & gpc > -xmax))
    meth.list[[i]] = meth.plt

    if (export){
        # get top sites
        meth.sig = meth.plt %>% 
            filter((cpg > 0.5 & gpc < -0.3 & direction == "under")|
                   (cpg < -0.5 & gpc > 0.3 & direction == "over"))%>%
            arrange(seqnames,start)
        bed.tb = as.tibble(window.sub[meth.sig$feature.index,]) 
        bed.tb = bed.tb[,c(1:3,6,9,5,13,7,8,10)]
        bedout=file.path(annodir,
                         paste0(prefix,
                                "_top_epigenetic_state_genes.TSS.2000bp.bed"))
        write_tsv(bed.tb,bedout,col_names=F)
        bed.tb = as.tibble(tss.sub[meth.sig$feature.index,]) 
        bed.tb = bed.tb[,c(1:3,6,9,5,13,7,8,10)]
        bedout=file.path(annodir,
                         paste0(prefix,
                                "_top_epigenetic_state_genes.TSS.bed"))
        write_tsv(bed.tb,bedout,col_names=F)
    }
    if (plot){
    # plot scatter
        plotpath = file.path(plotdir,paste0(prefix,"_scatter.pdf"))
    # labels for number in each quadrant
        idxlist = list(
            which(meth.plt$cpg>ythr & meth.plt$gpc>xthr),
            which(meth.plt$cpg>ythr & meth.plt$gpc<(-xthr)),
            which(meth.plt$cpg<(-ythr) & meth.plt$gpc<(-xthr)),
            which(meth.plt$cpg<(-ythr) & meth.plt$gpc>xthr))
        counts.list = lapply(seq_along(idxlist),function(i){
            meth.plt[idxlist[[i]],]%>%group_by(direction)%>%
                summarize(cnt=n())%>%mutate(quadrant=i)
        })
        counts = do.call(rbind,counts.list) %>%
            spread(direction,cnt) %>%
            mutate(text=paste0("over = ",over,"\nunder = ",under),
                   x=0,y=0)
        laby=ymax
        labx=xmax
        counts[1,] = counts[1,] %>% mutate(x=labx,y=laby)
        counts[2,] = counts[2,] %>% mutate(x=-labx,y=laby)
        counts[3,] = counts[3,] %>% mutate(x=-labx,y=-laby)
        counts[4,] = counts[4,] %>% mutate(x=labx,y=-laby)

        g = ggplot(meth.plt,aes(x=gpc,y=cpg,color=direction))+
            geom_point(alpha=0.5,size=0.5)+
            lims(y=c(-laby,laby),x=c(-labx,labx))+
            geom_hline(yintercept=0,linetype=2,size=0.5)+
            geom_vline(xintercept=0,linetype=2,size=0.5)+
            geom_text(inherit.aes=F,data=counts,size=3,
                      mapping=aes(x=x,y=y,label=text),
                      vjust="inward",hjust="inward",alpha=1)+
            theme_bw()+
            theme(legend.position="bottom",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color="black"))+
            labs(x="Accessibility",y="CpG methylation")+
            coord_fixed(labx/laby)
        

        pdf(plotpath,useDingbats=F,width=3,height=3)
        print(g)
        dev.off()
    }
}
# boxplot?
meth.allcomp = do.call(rbind,meth.list) %>%
    mutate(comp = paste0(onename,".vs.",twoname))
meth.gather = meth.allcomp %>%
    select(-seqnames,-start,-end,-width,-strand,
           -geneid,-p,-id,-score,-genename) %>%
    gather(key=feature,value=del,-feature.index,
           -onename,-twoname,-comp,-direction)

meth.over = meth.allcomp[which(meth.allcomp$direction == "over" &
                               meth.allcomp$comp == "MCF10A.vs.MCF7"),]
meth.under = meth.allcomp[which(meth.allcomp$direction == "under" &
                                meth.allcomp$comp == "MCF10A.vs.MCF7"),]

mean(meth.over$chrom)-mean(meth.under$chrom)
mean(meth.over$cpg)-mean(meth.under$cpg)
mean(meth.over$gpc)-mean(meth.under$gpc)

ttest.chrom = t.test(meth.over$chrom,meth.under$chrom)
ttest.cpg = t.test(meth.over$cpg,meth.under$cpg)
ttest.gpc = t.test(meth.over$gpc,meth.under$gpc)

ttest.chrom$p.value
ttest.gpc$p.value
ttest.cpg$p.value

plotpath = file.path(plotdir,
                     paste0("bcancomparison_by_expression_foldchange_",mindiff,
                            "_alpha_",a,"_boxplot.pdf"))
g = ggplot(meth.gather,aes(x=direction,
                           y=del,
                           color=comp))+
    facet_grid(.~feature)+
    geom_boxplot(outlier.size=-1)

pdf(plotpath,useDingbats=F)
print(g)
dev.off()
