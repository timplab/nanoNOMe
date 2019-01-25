#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(VariantAnnotation)
source("../../script/methylation_plot_utils.R")
cores=10
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

window=2000
    
i = 1
for (i in seq(dim(combos)[1])){
    oneidx = match(combos$one[i],pd$cell)
    twoidx = match(combos$two[i],pd$cell)
    
    genes.fp = file.path(root,"annotations/breastcancer",paste0("MCF10A_vs_",combos$two[i],"_genes.bed"))
    genes = read_tsv(genes.fp) %>%
        dplyr::rename(chrom=`#chrom`)
    genes.gr = GRanges(genes)
    
    vcf.sub = vcf[suppvec %in% sup.list[[i]]]
    bpall.list = genocoTogranges(geno(vcf.sub)$CO[,combos$ind[i]])
    lapply(bpall.list,function(x){
        dist = as.data.frame(distanceToNearest(x,genes.gr))
        dist[which(dist$distance<10000),]
    })
    
    vcf.tra = vcf.sub[info(vcf.sub)$SVTYPE=="TRA"]

    coords = geno(vcf.tra)$CO[,combos$ind[i]]
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
    
    # determine sv vs nonsv
    bp1names = reads$readname[which(reads$bp==1)]
    bp2names = reads$readname[which(reads$bp==2)]
    svnames = bp2names[bp2names %in% bp1names]
    reads.list = lapply(reads.list,function(x){x %>%
        mutate(type=ifelse(readname %in% svnames,"SV","nonSV"))})

    # compare SV in two and nonSV in one
    sv.reads = lapply(reads.list,function(x){
        x[which(x$type=="SV" & x$comp == 2),]})
    nonsv.reads = lapply(reads.list,function(x){
        x[which(x$type=="nonSV" & x$comp == 1),]})
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
    # scatterplot
    for (methwidth in c(100,400)){
        methranges = lapply(bprange.list,function(x){
            resize(x,width=methwidth,fix="center")})
        regmeth.list = lapply(seq_along(methranges),function(i){
            freq%>%group_by(calltype,type) %>%
            getRegionMeth(methranges[[i]]) %>%
            dplyr::select(-totcov,-numsites) %>%
                spread(type,freq)%>%mutate(bp=i)})
        regmeth = do.call(rbind,regmeth.list)

        plotscatter <- function(dat,title=title){
            ggplot(dat,aes(x=nonSV,y=SV))+
                geom_point()+
                lims(x=c(0,1),y=c(0,1))+
                theme_bw()+
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text = element_text(color="black"),
                    plot.title = element_text(size=10))+
                ggtitle(title)+
                coord_fixed()
        }
        mytheme <- function(g){
            g+
                theme_bw()+
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text = element_text(color="black"),
                    plot.title = element_text(size=10))+
                coord_fixed()
        }
        cpgfrac = round(dim(
            regmeth[which(regmeth$calltype=="cpg"),] %>%
            filter(SV<nonSV))[1]/dim(regmeth[which(regmeth$calltype=="cpg"),])[1],2)
        gpcfrac = round(dim(
            regmeth[which(regmeth$calltype=="gpc"),] %>%
            filter(SV<nonSV))[1]/dim(regmeth[which(regmeth$calltype=="gpc"),])[1],2)
        g.cpg = plotscatter(regmeth[regmeth$calltype=="cpg",],
                            paste0("cpg ; SV hypometh : ",cpgfrac))
        g.gpc = plotscatter(regmeth[regmeth$calltype=="gpc",],
                            paste0("gpc ; SV hypometh : ",gpcfrac))
        del.noSV = regmeth %>%
            dplyr::select(-SV) %>%
            spread(bp,nonSV) %>% na.omit() %>%
            mutate(del=`2`-`1`,type="noSV")
        del.SV = regmeth %>%
            dplyr::select(-nonSV) %>%
            spread(bp,SV) %>% na.omit() %>%
            mutate(del=`2`-`1`,type="SV")
        dels.tb = bind_rows(del.noSV,del.SV)
        dels.cpg = dels.tb[which(dels.tb$calltype=="cpg"),]
        dels.gpc = dels.tb[which(dels.tb$calltype=="gpc"),]
        dels.spread = dels.tb %>%
            dplyr::select(-`1`,-`2`) %>%
            spread(type,del) %>% na.omit()
        # plotting the bp comparison
        maxval = max(c(dels.cpg$`1`,dels.cpg$`2`))
        g.bpcomp.cpg = mytheme(ggplot(dels.cpg,aes(x=`1`,y=`2`))+
                               lims(x=c(0,maxval),y=c(0,maxval))+
                               facet_wrap(~type,ncol=2)+
                               geom_point()+
                               labs(x="BP1",y="BP2",title="CpG"))
        maxval = max(c(dels.gpc$`1`,dels.gpc$`2`))
        g.bpcomp.gpc = mytheme(ggplot(dels.gpc,aes(x=`1`,y=`2`))+
                               lims(x=c(0,maxval),y=c(0,maxval))+
                               facet_wrap(~type,ncol=2)+
                               geom_point()+
                               labs(x="BP1",y="BP2",title="GpC"))
        # plotting del
        delspread.cpg = dels.spread[which(dels.spread$calltype == "cpg"),]
        delspread.gpc = dels.spread[which(dels.spread$calltype == "gpc"),]
        maxval = max(abs(c(delspread.cpg$SV,delspread.cpg$noSV)))
        g.del.cpg = mytheme(ggplot(delspread.cpg,aes(x=noSV,y=SV))+
                            lims(x=c(-maxval,maxval),y=c(-maxval,maxval))+
                            geom_point()+
                            labs(x="MCF10A",y=combos$two[i],
                                 title="CpG BP2-BP1"))
        maxval = max(abs(c(delspread.gpc$SV,delspread.gpc$noSV)))
        g.del.gpc = mytheme(ggplot(delspread.gpc,aes(x=noSV,y=SV))+
                            lims(x=c(-maxval,maxval),y=c(-maxval,maxval))+
                            geom_point()+
                            labs(x="MCF10A",y=combos$two[i],
                                 title="GpC BP2-BP1"))
        plotpath = file.path(plotdir,paste0(combos$one[i],"_vs_",combos$two[i],"width",methwidth,"_scatter.pdf"))
        pdf(plotpath,3,3,useDingbats=F)
        print(g.bpcomp.cpg)
        print(g.bpcomp.gpc)
        print(g.del.cpg)
        print(g.del.gpc)
        dev.off()
    }
    if(FALSE){
        plotpath = file.path(plotdir,paste0(combos$one[i],"_vs_",combos$two[i],"TRA_byregion.pdf"))
        pdf(plotpath,3,3,useDingbats=F)
        for (j in seq_along(vcf.tra)){
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
                ggtitle(paste("bp1 :",toString(bp1),"; bp2 :",toString(bp2)))
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
                ggtitle(paste("bp1 :",toString(bp1),"; bp2 :",toString(bp2)))
            print(g)
                                        #    print(gcov)
        }
        dev.off()
    }
}
