#!/usr/bin/Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
source(file.path(srcdir,"../../scripts/ggplot_theme.R"))
library(GenomicRanges)
library(VariantAnnotation)
source(file.path(srcdir,"../../scripts/methylation_plot_utils.R"))

cores=detectCores()
root="/kyber/Data/Nanopore/projects/nanonome/analysis"
setwd(root)

## data paths
pd.all = tibble(cell=c("MCF10A","MCF10A","MCF7","MCF7","MDAMB231","MDAMB231"),
                calltype=c("cpg","gpc","cpg","gpc","cpg","gpc"),
                fpath=c(paste0("data/nanonome/pooled/mfreq/",cell,
                               "_nanoNOMe.pooled.",calltype,".mfreq.txt.gz")))

## sv beds
for (i in c(1,2)){
    for (j in c(1,2,3)){
        
        combos = tibble(onename="MCF10A",twoname=c("MCF7","MDAMB231"))
        svtypes = c("del","ins","tra")
        pd = pd.all[pd.all$cell %in% c(combos$onename[i],combos$twoname[i]),] %>%
            mutate(cell = ifelse(cell==combos$onename[i],"sampleone","sampletwo"))

        bedfp = paste0("data/bcan/MCF10A_vs_",combos$twoname[i],
                       "_SVcomparison_",svtypes[j],"_survivor.flank200bp.bed")
        shufflefp = paste0("data/bcan/MCF10A_vs_",combos$twoname[i],
                           "_SVcomparison_",svtypes[j],"_survivor.shuffle.flank200bp.bed")
        if (svtypes[j] == "tra"){
            bedfp = paste0("data/bcan/MCF10A_vs_",combos$twoname[i],
                           "_SVcomparison_",svtypes[j],"_survivor.TRAregion200bp.bed")
            shufflefp = paste0("data/bcan/MCF10A_vs_",combos$twoname[i],
                               "_SVcomparison_",svtypes[j],"_survivor.shuffle.TRAregion200bp.bed")
        }

        ## output plot path
        plotpath = paste0("plots/MCF10A_vs_",combos$twoname[i],
                          "_SVcomparison_",svtypes[j],"_survivor.flank200bp.plot.pdf")

        pd.bed = pd %>% mutate(regname="bed",regfp=bedfp)
        pd.shuffle = pd %>% mutate(regname="shuffle",regfp=shufflefp)
        pd = bind_rows(pd.bed,pd.shuffle)

        ## read in data
        dat.list = mclapply(mc.cores=cores,seq(dim(pd)[1]),function(i){
            tabix_mfreq(pd$fpath[i],pd$regfp[i]) %>%
                mutate(cell=pd$cell[i],
                       calltype=pd$calltype[i],
                       reg=pd$regname[i])
        })

        dat = do.call(rbind,dat.list)
        dat.gr = GRanges(dat)
        bed.gr = load_db(bedfp,"updown")
        bed.gr$reg = "bed"
        shuffle.gr = load_db(shufflefp,"updown")
        shuffle.gr$reg = "shuffle"
        reg.gr = c(bed.gr,shuffle.gr)

        ## distance
        reg.up = reg.gr[which(reg.gr$updown == 1)]
        reg.down = reg.gr[which(reg.gr$updown == 2)]
        idx = match(paste(reg.down$id,reg.down$score,reg.down$reg),paste(reg.up$id,reg.up$score,reg.up$reg))
        reg.upmatch = reg.up[idx]
        dists = tibble(id = reg.upmatch$id,
                       cell = reg.upmatch$score,
                       dist = distance(reg.upmatch,reg.down))
        bins=c(0,500,1e3,1e4,1e5,Inf)
        levels=c("<500","500-1000","1e3-1e4","1e4-1e5",">1e5")
        dists$dist = cut(x=dists$dist,breaks=bins,labels=levels)
        dists$dist[is.na(dists$dist)] = levels[5]

        ## which into what
        ovl = findOverlaps(dat.gr,reg.gr)
        dat.ovl = bind_cols(dat[queryHits(ovl),],
                            id = reg.gr$id[subjectHits(ovl)],
                            updown = reg.gr$updown[subjectHits(ovl)],
                            svcell = reg.gr$score[subjectHits(ovl)],
                            reg = reg.gr$reg[subjectHits(ovl)])
        ## per region methylation
        dat.sum = dat.ovl %>%
            group_by(cell,calltype,id,updown,svcell,reg) %>%
            summarize(meth=sum(meth),
                      unmeth=sum(unmeth),
                      mfreq = meth/(meth+unmeth))
        ## delta between upstream and downstream
        dat.del = dat.sum %>% ungroup() %>%
            dplyr::select(-meth,-unmeth) %>%
            spread(updown,mfreq) %>% na.omit() %>%
            mutate(del=`2`-`1`) %>%
            filter(svcell==combos$twoname[i]) # only those occurring in cancer
        ## match sv by cell
        dat.cell = dat.del %>%
            dplyr::select(-`1`,-`2`) %>%
            spread(cell,del) %>% na.omit()

        ## add distance
        idx = match(paste(dat.cell$id,dat.cell$svcell),paste(dists$id,dists$cell))
        dat.cell$dist = dists$dist[idx]

        ## plot
        scatter <- function(data){
            ## linear regression
            lm.tb = data %>%
                group_by(calltype,reg) %>%
                do(model = lm(sampletwo ~ sampleone,data=.))
            lm.tb = lm.tb %>%
                mutate(rsq = summary(model)$r.squared,
                       slope = summary(model)$coefficients[,"Estimate"]["sampleone"],
                       intercept = summary(model)$coefficients[,"Estimate"]["(Intercept)"],
                       regrsq = paste(reg,":",round(rsq,2)),
                       svcell = data$svcell[1])
            data$regrsq = lm.tb$regrsq[match(data$reg,lm.tb$reg)]
            ggplot(data,aes(x=sampleone,y=sampletwo,color=dist))+
                facet_wrap(regrsq~svcell) +
                geom_point(alpha=0.5)+
                geom_abline(mapping=aes(slope=slope,intercept=intercept),
                            data=lm.tb,linetype=2)+
                coord_fixed()
        }
        boxplot <- function(data){
            ggplot(data,aes(x=reg,y=del))+
                facet_wrap(cell~svcell)+
                geom_boxplot()
        }
        scatter.cpg = scatter(dat.cell[which(dat.cell$calltype=="cpg"),]) +
            lims(x=c(-1,1),y=c(-1,1))+ggtitle("cpg")
        scatter.gpc = scatter(dat.cell[which(dat.cell$calltype=="gpc"),]) +
            lims(x=c(-.5,.5),y=c(-.5,0.5))+ggtitle("gpc")
        box.cpg = boxplot(dat.del[which(dat.del$calltype=="cpg"),]) +
            ggtitle("cpg")
        box.gpc = boxplot(dat.del[which(dat.del$calltype=="gpc"),]) +
            ggtitle("gpc")

        pdf(plotpath,width=4,height=3)
        print(scatter.cpg)
        print(scatter.gpc)
        print(box.cpg)
        print(box.gpc)
        dev.off()
    }

}
