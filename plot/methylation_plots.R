#!/usr/bin/Rscript
# script is for plotting current distribution using subsetting eventalign files

parseArgs = function(argsin){
    library(optparse)
    arglist=list(
        make_option(c("-i","--input"),type="character",default=NULL,
                    help="input file path",metavar="/path/to/input"),
        make_option(c("-o","--output"),type="character",default=NULL,
                    help="output file path",metavar="/path/to/out")
        )
    
    argparser=OptionParser(option_list=arglist)
    args=parse_args(argparser,args=argsin)
    ##Check if there are any variables, and if not, show help
    if (is.null(args$input)||is.null(args$output)){
        print_help(argparser)
        stop("input and output must be provided",call.=FALSE)
    }
    args
}

methByDistance = function(infile,outfile){
    library(tidyverse)
    cnames=c("distance","meth","unmeth")
    distmeth=read_tsv(file=infile,col_names=cnames)
    b=20
    meth.sum=mutate(distmeth,bindist=round(distance/b)*b)%>%
        group_by(bindist)%>%
        summarize(cov=n(),
                  meth=sum(meth),
                  unmeth=sum(unmeth),
                  freq=meth/(meth+unmeth))
    # plot
    g=ggplot(meth.sum,aes(x=bindist,y=freq))+theme_bw()+
        geom_point(size=0.5)+lims(y=c(0,1))+
        labs(title="Methylation by distance",
             x="Distance", y="Average Methylation")
    gsmooth=g+geom_smooth(method="loess",span=0.1,se=F)
    gline=g+geom_line()
    gcov=ggplot(meth.sum,aes(x=bindist,y=cov))+theme_bw()+
        geom_point(size=0.5)+labs(title="Coverage by distance",
                                  x="Distance",
                                  y="Coverage")
    pdf(outfile,width=10,height=4,useDingbats=F)
    print(gline)
    print(gsmooth)
    print(gcov)
    dev.off()
}

if (! interactive()){
    argsin=commandArgs(TRUE)
    command=argsin[1]
    args=parseArgs(argsin[-1])
    if( command=="methByDistance"){
        methByDistance(args$input,args$output)
    } else if( command == "null"){
        print("null")
    }
    
}

  
  

