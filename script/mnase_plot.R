#!/usr/bin/Rscript
parse_mnase = function(argsin){
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

metaplotByDistance <- function(inpath,plotpath){
    # read data
    cat("reading data\n")
    script = file.path(dirname(get_Rscript_filename()),"parseMNase.py")
    com = paste0("python ",script," by-distance -v -i ",inpath)
    dat.raw=system(com,intern=T)
    dat.tb = as.tibble(do.call(rbind,strsplit(dat.raw,"\t")))%>%
        type_convert()
    names(dat.tb) = c("dist","freq")
    ymax = max(dat.tb$freq)
    # aggregate
    cat("aggregating\n")
    win = 50
    roll.range=seq(from=min(dat.tb$dist)+win/3,
                   to=max(dat.tb$dist)-win/2,by=1)
    dat.roll = tibble(Distance = roll.range,
                      Signal = calculate_rollmean(dat.tb,roll.range,win))
    # plot
    cat("plotting\n")
    g=ggplot(dat.roll,aes(x=Distance,y=Signal))+
        theme_bw()+geom_line()+
        lims(y=c(0,ymax))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"))+
        labs(x="Distance to center", y="Average Signal")
    pdf(plotpath,width=6,height=4,useDingbats=F)
    print(g)
    dev.off()
}
if (! interactive()){
    modules=c("metaplotByDistance")
    argsin=commandArgs(TRUE)
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
    library(GenomicRanges)

    if( command=="metaplotByDistance"){
        args=parse_mnase(argsin[-1])
        metaplotByDistance(args$input,args$plotpath)
   }
}
