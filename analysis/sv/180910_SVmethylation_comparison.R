#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../plot/methylation_plot_utils.R")
cores=10
parser="../../script/parseMethylFreq.py"
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/mfreq_all")
svroot=file.path(root,"pooled/sv")
plotdir=file.path(root,"plots/sv")
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
samples = tibble(cell=cells,
            cpg=file.path(datroot,paste(cell,"cpg.methfreq.txt.gz",sep=".")),
            gpc=file.path(datroot,paste(cell,"gpc.methfreq.txt.gz",sep=".")))
bedpaths = tibble(cell=cells,
                 fpath = file.path(svroot,paste0(cell,".svregion.250.bed")))
pd.list = lapply(seq_along(bedpaths$fpath),function(i){
    samples %>% mutate(bedcell=bedpaths$cell[i],bedpath=bedpaths$fpath[i])
})
pd = do.call(rbind,pd.list)%>%
    gather(calltype,filepath,-cell,-bedpath,-bedcell)

comp.path = file.path(svroot,"SVcomparison.vcf")

# read in regs
extracols = c("svtype","part","location")
reg.list = lapply(bedpaths$fpath,function(x){
    load_db(x,extracols)
})
# read in comparison
comp = read_tsv(comp.path,col_names=F)
comp.filt = comp[which(comp$X11 %in% c("xxxo") & comp$X5 == "<TRA>"),]
# read in avg meth
coms = lapply(seq(dim(pd)[1]),function(i){
    paste0("head -n100 ",pd$bedpath[i]," | python -u ",parser," by-region -v -t 12 -i ",pd$filepath[i])
})
coms = lapply(seq(dim(pd)[1]),function(i){
    paste0("python -u ",parser," by-region -v -t 12 -i ",pd$filepath[i]," -r ",pd$bedpath[i])
})

cnames = c("chrom","start","end","id","score","strand",extracols,"freq","cov","numsites")
dat.raw = lapply(coms,function(x){
    print(x)
    y = system(x,intern=T)
    out = as.tibble(do.call(rbind,strsplit(y,"\t")))
    names(out) = cnames;out
})

dat.named = lapply(seq_along(dat.raw),function(i){
    dat.raw[[i]]%>%type_convert()%>%
        mutate(cell=pd$cell[i],
               calltype=pd$calltype[i],
               svcell=pd$bedcell[i])
})

dat.all = do.call(rbind,dat.named)
dat.spread = dat.all%>%
    filter(cov>5 & location != "breakpoint" & svcell=="MDAMB231")%>%
    select(-cov,-numsites)%>%
    spread(cell,freq)
length(unique(paste0(dat.spread$id,dat.spread$svcell)))

dat.comp = dat.spread[which(dat.spread$id %in% comp.filt$X3),] %>%
    mutate(del = MDAMB231-MCF10A,lab = paste(part,location,sep=" : ")) %>%
    arrange(desc(del))
reg.insert = as.tibble(as.data.frame(reg.list[[4]]))%>%filter(location=="region")
dat.comp$svlen = reg.insert$width[match(dat.comp$id,reg.insert$id)]

dat.comp %>% group_by(calltype,part,location)%>%
    summarize(n())
length(unique(paste0(dat.comp$id,dat.comp$svcell)))
# spread to get cpg vs gpc
dat.chrom = dat.comp %>%
    select(id,part,location,calltype,del,lab)%>%
    spread(calltype,del)

# plotting
g = ggplot(dat.comp)+ theme_bw()
g.gpc = ggplot(dat.comp[dat.comp$calltype=="gpc",])+theme_bw()

g.scatter = g + geom_point(aes(x=MCF10A,y=MDAMB231),inherit.aes=F)+
    geom_abline(slope=1,linetype="dotted")+
    facet_grid(lab~calltype)

g.scatter.gpc = g.gpc + geom_point(aes(x=MCF10A,y=MDAMB231),inherit.aes=F)+
    geom_abline(slope=1,linetype="dotted")+
    facet_grid(lab~calltype)

g.box = g + geom_boxplot(aes(x=paste(lab,calltype),y=del,group=paste(lab,calltype),color=part))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(y="MDA231-MCF10A")
g.freqpoly = g + geom_freqpoly(aes(x=del),inherit.aes=F)+
    facet_grid(lab~calltype)+
    labs(x="MDA231-MCF10A")
g.chrom = ggplot(dat.chrom,aes(x=cpg,y=gpc))+theme_bw()+
    geom_point()+
    facet_grid(lab~.)+
    ggtitle("MDAMB231-MCF10A gpc vs cpg")

plotpath = file.path(plotdir,"180912_svcomparison_mda231TRA_mda231minusmcf10a_250bpwindow.pdf")
pdf(plotpath,height=10,width=4,useDingbats=FALSE)
print(g.scatter.gpc)
print(g.chrom)
print(g.scatter)
print(g.box)
print(g.freqpoly)
dev.off()
