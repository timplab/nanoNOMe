#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
source("../../script/methylation_plot_utils.R")
cores=10
parser="../../script/SVmethylation.py"
root="/dilithium/Data/Nanopore/projects/nomeseq/analysis"
datroot=file.path(root,"pooled/methylation/methbyread_all")
plotdir=file.path(root,"plots/sv")
cells=c("GM12878","MCF10A","MCF7","MDAMB231")
pd = tibble(cell=cells,
            cpgpath=file.path(datroot,paste(cell,"cpg.pooled.meth.bed.gz",sep=".")),
            gpcpath=file.path(datroot,paste(cell,"gpc.pooled.meth.bed.gz",sep=".")),
            bampath=file.path(root,"pooled/bam",paste0(cell,".pooled.bam")),
            svpath=file.path(root,paste0("pooled/sv/",cell,".sniffles.vcf")))
pd.gather = pd %>%
    gather(calltype,path,cpgpath,gpcpath)%>%
    mutate(calltype=gsub("path","",calltype))
    
svcomp=file.path(root,paste0("pooled/sv/SVcomparison.vcf"))
# read in average methylation
coms = lapply(seq(dim(pd)[1]),function(i){
    paste("grep -E 'xxoo|xxxo'",svcomp,"| python",parser,"-v -t 10 -b",pd$bampath[i],
          "-c",pd$cpgpath[i],"-g",pd$gpcpath[i])
})
raw.list = lapply(coms,function(x){
    system(x,intern=T)})
cnames = c("chrom1","start1","end1","chrom2","start2","end2","readname","svid","strand1","strand2","sv","tag","cpgcnt","gpccnt","cpgcov","gpccov")
dat.list = mclapply(mc.cores=cores,seq_along(raw.list),function(i){
    x = as.tibble(do.call(rbind,strsplit(raw.list[[i]],"\t")))
    colnames(x) = cnames
    x %>% mutate(cell=pd$cell[i])
})
dat.all = do.call(rbind,dat.list) %>% type_convert()

# pick ones with decent methylation
dat.sum = dat.all%>%group_by(cell,svid,sv,tag)%>%
    summarize(numreads=n(),
              totcpgcov=sum(cpgcov),
              cpg = sum(cpgcnt)/totcpgcov,
              totgpccov=sum(gpccov),
              gpc = sum(gpccnt)/totgpccov)%>%
    filter(totcpgcov>20)

dat.del = dat.sum %>%ungroup()%>%select(cell,svid,sv,tag,cpg)%>%
    filter(cell=="MCF10A" | cell=="MDAMB231")%>%
    spread(cell,cpg)%>%na.omit()%>%
    mutate(del=MDAMB231-MCF10A)%>%
    arrange(desc(abs(del)))

id.order = unique(dat.del$svid)

# make granges for reading methylation
sv.tb = dat.all %>% group_by(chrom1,start1,end1,chrom2,start2,end2,svid)%>%
    summarize(numreads=n())%>%
    filter(svid %in% id.order)
sv1 = sv.tb[,c(1:3,7)]%>% mutate(context="target")%>%
    setNames(c("chrom","start","end","id","context"))
sv2 = sv.tb[,c(4:7)]%>%mutate(context="insert")%>%
    setNames(names(sv1))
regs=rbind(sv1,sv2)
regs.gr=GRanges(regs)
win=2000
start(regs.gr) = start(regs.gr)-win
end(regs.gr) = end(regs.gr)+win

# reading methylation data
numplot = 100
regs.sub = regs.gr[which(regs.gr$id %in% id.order[1:numplot])]

meth.svreg = mclapply(mc.cores=cores,seq(dim(pd.gather)[1]),function(i){
    tabix_mbed(pd.gather$path[i],dbpath=regs.sub,by="call")%>%
        mutate(cell=pd.gather$cell[i],calltype=pd.gather$calltype[i])
})
meth.all = do.call(rbind,meth.svreg) %>%
    na.omit()
# determine context and plot for each sv
meth.gr = GRanges(meth.all)
plotwidth = 2000
plotpath = file.path(plotdir,"SVcomparison_readlevel.pdf")

pdf(plotpath,useDingbats=F,width=12,height=6)
for ( i in seq(1,100)){
    id = id.order[i]
    queries = dat.all[which(dat.all$svid==id),]
    reg = resize(regs.gr[regs.gr$id==id.order[i]],width=plotwidth,fix="center")
    ovlidx=overlapsAny(meth.gr,reg)
    contextidx = subjectHits(findOverlaps(meth.gr[ovlidx],reg))
    svcontext = reg$context[contextidx]
    center=start(reg)[contextidx]+plotwidth/2
    meth = meth.all[ovlidx,]%>%mutate(tag=svcontext,
                                      dist=start-center)

    dat.match = bind_cols(meth,queries[match(meth$qname,queries$readname),"sv"]) %>%
        na.omit()

    dat.avg = dat.match %>% group_by(chrom,dist,cell,sv,tag,calltype)%>%
        summarize(cov=n(),
                  freq=sum(mcall)/cov) %>%
        na.omit()%>%
        mutate(lab=paste(cell,sv,tag,calltype,sep="_"))
    dat.sum = dat.avg %>% group_by(sv,tag,calltype)%>%
        summarize(cov=mean(cov))
    if(sum(dat.sum$cov[which(dat.sum$sv=="noSV")]<5)!=0){next}
    # plot
    print(i)
    title=gsub(", "," > ",toString(resize(reg,width=1,fix="center")))
    g.avg = ggplot(dat.avg,aes(dist,freq,group=calltype,color=calltype))+theme_bw()+
        facet_grid(cell~sv+tag)+
        lims(x=c(-plotwidth/2,plotwidth/2),y=c(0,1))+
        geom_point(size=0.5,alpha=0.7)+
        geom_smooth(se=F,size=0.5)+
        ggtitle(title)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"))
    
    print(g.avg)
}
dev.off()
