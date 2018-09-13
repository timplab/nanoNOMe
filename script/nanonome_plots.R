#!/usr/bin/Rscript
# nanoNOMe plots

parse_mfreq = function(argsin){
    arglist=list(
        make_option(c("-c","--cpg"),dest="cpgpath",type="character",default=NULL,
                    help="cpg methylation frequency file path",metavar="/path/to/cpg"),
        make_option(c("-g","--gpc"),dest="gpcpath",type="character",default=NULL,
                    help="gpc methylation frequency file path",metavar="/path/to/gpc"),
        make_option(c("-r","--region"),dest="dbpath",type="character",default=NULL,
                    help="bed file of region(s) to plot"),
        make_option(c("-o","--output"),dest="plotpath",type="character",default=NULL,
                    help="output plot pdf file"),
        make_option(c("-w","--window"),type="numeric",default=50,
                    help="rolling mean window")
    )
    argparser=OptionParser(option_list=arglist)
    args=parse_args(argparser,args=argsin)
    args
}

# plotting
#plotpath  = args$plotpath
#win = args$window    
aggregateByDistance <- function(cpg.tb,gpc.tb,db.gr,plotpath,win){
    db.center = resize(shift(db.gr,shift=width(db.gr)/2),width=1,ignore.strand=T)
    dat.list=list(cpg.tb,gpc.tb)
    dat.ag=lapply(dat.list,function(x){
        x.filt=x%>%filter(cov>0 & trinuc!="GCG" & trinuc != "CCG")
        x.dist=x.filt%>%
            bind_cols(getDistance(x.filt,db.center))
        aggregate_methylation(x.dist,win)
    })

    dat.ag[[1]]$lab="CpG Methylation"
    dat.ag[[2]]$lab="GpC Accessibility"

    dat.plt=do.call(rbind,dat.ag)
    
    # plot
    g=ggplot(dat.plt,aes(x=dist,y=freq,color=lab,group=lab))+
        theme_bw()+geom_line()+
        lims(y=c(0,1))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"))+
        labs(x="Binned distance", y="Average methylation")

    pdf(plotpath,width=6,height=4,useDingbats=F)
    print(g)
    dev.off()
}
heatmapByDistance <- function(cpg.tb,gpc.tb,db.gr,plotpath,win){
    library(parallel)
    db.center = resize(shift(db.gr,shift=width(db.gr)/2),width=1,ignore.strand=T)
    dat.list=list(cpg.tb,gpc.tb)
    dat.dist=lapply(dat.list,function(x){
        x.noGCG=x%>%filter(trinuc!="GCG")
        x.dist=x.noGCG%>%mutate(freq=meth/(meth+unmeth))%>%
            bind_cols(getDistance(x.noGCG,db.center))
    })
    by=win/2
    cat("rollmean\n")
    dat.rollmean=lapply(dat.dist,function(x){
        rollrange=seq(from=min(x$dist)+win/2,to=max(x$dist)-win/2,by=by)
        rollmean.list=mclapply(seq_along(db.center),mc.cores=8,FUN=function(i){
            x.reg=x[which(x$index.feature==i),]
            calculate_rollmean(x.reg,rollrange,win=win/2)
        })
        out=as.tibble(do.call(rbind,rollmean.list))
        colnames(out)=rollrange
        out
    })

    cat("tidy data\n")
    dat.cpg=dat.rollmean[[1]]
    dat.gpc=dat.rollmean[[2]]
    order.tb=tibble(idx=order(rowMeans(dat.gpc[,-1],na.rm=T)),
                    order=seq(length(db.center)))%>%arrange(idx)
    dat.cpg$ind=dat.gpc$ind=order.tb$order
    cpg.plt=dat.cpg%>%
        gather(key=x,value=freq,-ind)%>%
        na.omit()%>%
        mutate(x=as.numeric(x))

    gpc.plt=dat.gpc%>%
        gather(key=x,value=freq,-ind)%>%
        na.omit()%>%
        mutate(x=as.numeric(x))
    # plot
    cat("plot\n")
    g.cpg=ggplot(cpg.plt)+
        geom_rect(aes(xmin=x-(win/2),xmax=x+(win/2),
                      ymin=ind,ymax=ind+1,fill=freq))+
        theme_bw()+
        labs(x="Binned distance", y="Average methylation")
    g.gpc=ggplot(gpc.plt)+
        geom_rect(aes(xmin=x-(win/2),xmax=x+(win/2),
                      ymin=ind,ymax=ind+1,fill=freq))+
        theme_bw()+
        labs(x="Binned distance", y="Average methylation")


    pdf(plotpath,width=5,height=5,useDingbats=F)
    print(g.cpg)
    print(g.gpc)
    dev.off()
}

test="aggregateByDistance,-c,/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation/mfreq_all/GM12878.cpg.methfreq.txt.gz,-g,/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/methylation/mfreq_all/GM12878.gpc.methfreq.txt.gz,-r,/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/ctcf/GM12878_CTCF.2kb.bed,-o,/dilithium/Data/Nanopore/projects/nomeseq/analysis/plots/aggregate/GM12878.ctcf.aggregate.pdf"
argsin=strsplit(test,",")[[1]]
if (! interactive()){
    modules=c("aggregateByDistance","heatmapByDistance")
    argsin=commandArgs(TRUE)
#    print(paste(argsin,collapse=","))
#    quit()
    command=argsin[1]
    if (command == "-h" | command == "--h" |
        ! command %in% modules){
        cat(paste("Possible modules :",paste(modules,collapse=", "),"\n"))
        quit()
    }
    library(optparse)
    library(getopt)
    source(file.path(dirname(get_Rscript_filename()),"methylation_plot_utils.R"))
    library(tidyverse)
    library(cluster)
    library(GenomicRanges)

    if( command=="aggregateByDistance"){
        args=parse_mfreq(argsin[-1])
        db.gr=load_db(args$dbpath)
        cpg.tb=tabix_mfreq(args$cpgpath,args$dbpath)
        gpc.tb=tabix_mfreq(args$gpcpath,args$dbpath)
        cat("plotting\n")
        aggregateByDistance(cpg.tb,gpc.tb,db.gr,args$plotpath,args$window)
    } else if( command == "heatmapByDistance"){
        args=parse_mfreq(argsin[-1])
        db.gr=load_db(args$dbpath)
        cpg.tb=tabix_mfreq(args$cpgpath,args$dbpath)
        gpc.tb=tabix_mfreq(args$gpcpath,args$dbpath)
        heatmapByDistance(cpg.tb,gpc.tb,db.gr,args$plotpath,args$window)

    } else if( command == "readlevelByDistance"){
        arglist=args_meth()
        args=parseArgs(arglist,argsin[-1])
        data.tb=readIntersectReadlevel(args$input)
        readlevelByDistance(data.tb,args$output,args$coverage)
   }
}
