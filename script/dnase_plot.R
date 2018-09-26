#!/usr/bin/Rscript
parse_dnase = function(argsin){
    arglist=list(
        make_option(c("-i","--input"),dest="input",type="character",default=NULL,
                    help="signal bedgraph.gz file"),
        make_option(c("-o","--output"),dest="plotpath",type="character",default=NULL,
                    help="output plot pdf file")
    )
    argparser=OptionParser(option_list=arglist)
    args=parse_args(argparser,args=argsin)
    args
}

aggregateByDistance <- function(inpath,plotpath){
    script=file.path(dirname(get_Rscript_filename()),"parseDNAse.py")
#    script=("./parseDNAse.py")
    com = paste0("python ",script," by-distance -v -i ",
                 inpath)
#    com = paste0("python ./parseDNAse.py by-distance -v -i ",inpath)
    dat.raw = system(com,intern=T)
    dat.tb = as.tibble(do.call(rbind,strsplit(dat.raw,"\t"))) %>%
        type_convert()
    names(dat.tb) = c("dist","freq")
    win = 50
    roll.range=seq(from=min(dat.tb$dist)+win/3,
                   to=max(dat.tb$dist)-win/2,by=1)
    dat.roll = tibble(Distance = roll.range,
                      Signal = calculate_rollmean(dat.tb,roll.range,win))
    names(dat.tb) = c("Distance","Signal")
    # plot
    g=ggplot(dat.roll,aes(x=Distance,y=Signal))+
        lims(y=c(0,2))+
        theme_bw()+geom_line()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"))+
        labs(x="Distance to center", y="Average Signal")

    pdf(plotpath,width=6,height=4,useDingbats=F)
    print(g)
    dev.off()
}
test="aggregateByDistance,-i,/dilithium/Data/Nanopore/projects/nomeseq/analysis/distancebed/GM12878_dnase_distance_ctcf.bedGraph,-r,/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/ctcf/GM12878_ctcf.center.bed,-o,/dilithium/Data/Nanopore/projects/nomeseq/analysis/plots/aggregate/GM12878_dnase_aggregate_ctcf.pdf"
argsin = strsplit(test,",")[[1]]
if (! interactive()){
    modules=c("aggregateByDistance")
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
#    source("methylation_plot_utils.R")
    library(tidyverse)
    library(cluster)
    library(GenomicRanges)

 
    if( command=="aggregateByDistance"){
        args=parse_dnase(argsin[-1])
        cat("plotting\n")
        aggregateByDistance(args$input,args$plotpath)
   }
}
