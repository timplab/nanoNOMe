#!/usr/bin/Rscript
# methylation utils

# data reading
read_data <- function(infp){
    if (infp == "stdin"){
        infp <- file(infp)
    }
    read_tsv(infp,col_names=F)
}

load_db <- function(fpath,extracols=c("regstart","regend")){
    bedcnames=c("chrom","start","end","id","score","strand")
    cnames=c(bedcnames,extracols)
    clen=count_fields(fpath,tokenizer_tsv())[1]
    bed.tb=read_tsv(fpath,col_names=cnames[1:clen])%>%
        mutate(start=start+1)
    bed.gr=GRanges(bed.tb)
    bed.gr
}
tabix <- function(querypath,dbpath,col_names=NULL,verbose=TRUE){
        if ("GRanges" %in% class(dbpath)){
            # input region is a GRanges object
            if (verbose) cat(paste0("reading regions defined by GRanges object",
                                    " in ",querypath,"\n"))
            corenum=detectCores()-2
            if (verbose) cat(paste0("using ",corenum," cores\n"))
            raw.list = mclapply(mc.cores=corenum,dbpath,function(reg.gr){
                if (verbose) cat(paste0(as.character(reg.gr),"\n"))
                reg.char=paste0(as.character(seqnames(reg.gr)),":",
                                as.character(start(reg.gr)),"-",
                                as.character(end(reg.gr)))
                command=paste("tabix",querypath,reg.char,sep=" ")
                region.raw = system(command,intern=TRUE)
                if (verbose) cat("converting to tibble\n")
                region=do.call(rbind,strsplit(region.raw,"\t"))
                region.tb=as.tibble(region)
                if (!is.null(col_names)) colnames(region.tb)=col_names
                region.tb
            })
            if (verbose) cat("joining\n")
            region.tb = do.call(rbind,raw.list)
        } else {
            if (verbose) cat(paste0("reading regions defined by ",
                                    dbpath," in ",querypath,"\n"))
            command=paste("tabix",querypath,"-R",dbpath,sep=" ")
            region.raw=system(command,intern=TRUE)
            if (verbose) cat("converting to tibble\n")
            region=do.call(rbind,strsplit(region.raw,"\t"))
            region.tb=as.tibble(region)
            if (!is.null(col_names)) colnames(region.tb)=col_names
        }
        region.tb %>% type_convert()
}

mbedByCall <- function(mbed,smooth=FALSE,ns=10,h=50,verbose=T){
    if (verbose) cat("parsing data into single call per line\n")
    corenum=(detectCores()-2)/2
    if (verbose) cat(paste0("using ",corenum," cores\n"))
    out.list=mclapply(mc.cores=corenum,seq(dim(mbed)[1]),function(i){
        if (verbose) cat(paste0(i,"\n"))
        read=mbed[i,]
        mstring=read$mstring    
        call=str_extract_all(mstring,"[a-z]")[[1]]
        pos=cumsum(as.numeric(strsplit(mstring,"[a-z]")[[1]]))
        out=tibble(chrom=read$chrom,
               start=read$start+pos+1,
               end=start,
               qname=read$readname,
               mcall=call,
               score=as.numeric(strsplit(read$scores,",")[[1]]),
               context=strsplit(read$context,",")[[1]]
               )
        out$mcall[which(out$mcall=="m")]=1
        out$mcall[which(out$mcall=="u")]=0
        out$mcall=as.numeric(out$mcall)
        if ( dim(out)[1] == 0 ) return(out)
        if ( smooth == TRUE ){
            require(locfit)
            smooth = locfit(mcall ~ lp(start, nn = ns, h = h), data = out, maxk =10000)
            pp = preplot(smooth,where="data",band="local")
            out$smooth = pp$trans(pp$fit)
        }
        if (verbose) cat(paste0("done with ",i,"\n"))
        out
    })
    calls = do.call(rbind,out.list)
    calls
}
    
tabix_mbed <- function(querypath,dbpath=NULL,by=c("read","call"),verbose=TRUE){
    mbedcnames=c("chrom","start","end","readname","mstring","scores","context")
    if (!is.null(dbpath)) {
        region.tb=tabix(querypath,dbpath,mbedcnames,verbose=verbose)
        region.tb$start=as.numeric(region.tb$start)
        region.tb$end=as.numeric(region.tb$end)
        if (verbose) cat("removing redundant loci\n")
        out.tb=unique(region.tb)
    }else{
        if (verbose) cat(paste0("reading the entire data of ",querypath,"\n"))
        out.tb=read_tsv(querypath,col_names=mfreqcnames)
    }
    parsembed=function(mbed.tb,what,verbose=TRUE){
        if (what == "call"){
            mbedByCall(out.tb,verbose=verbose)
        } else if (what == "read"){
            mbed.tb
        }
    }
    if (length(by)>1){
        lapply(by,function(x){parsembed(out.tb,x,verbose)})
    }else{ parsembed(out.tb,by) }

}
    
tabix_mfreq <- function(querypath,dbpath=NULL,cov=2,trinuc_exclude="GCG",verbose=TRUE){
    mfreqcnames=c("chrom","start","strand","meth","unmeth","dinuc","trinuc")
    if (!is.null(dbpath)) {
        region.tb=tabix(querypath,dbpath,mfreqcnames,verbose=verbose)
        region.tb$meth=as.numeric(region.tb$meth)
        region.tb$unmeth=as.numeric(region.tb$unmeth)
        region.tb$start=as.numeric(region.tb$start)
        if (verbose) cat("removing redundant loci\n")
        out.tb=unique(region.tb)
    }else{
        if (verbose) cat(paste0("reading the entire data of ",querypath,"\n"))
        out.tb=read_tsv(querypath,col_names=mfreqcnames)
    }
    if (!is.null(trinuc_exclude)){
        if (verbose) cat(paste0("removing ",trinuc_exclude,"\n"))
        out.tb=out.tb[which(out.tb$trinuc!=trinuc_exclude),]
    }
    if (verbose) cat("calculating coverage and frequency\n")
    out.tb$cov=out.tb$meth+out.tb$unmeth
    if (!is.null(cov)){
        if (verbose) cat(paste0("removing sites with less than ",cov,"x coverage\n"))        
        out.tb=out.tb[which(out.tb$cov>=cov),]
    }
    out.tb$freq=out.tb$meth/out.tb$cov
    out.tb$end=out.tb$start
    na.omit(out.tb)
}
# functions
getCenter <- function(db.gr){
    resize(shift(db.gr,shift=width(db.gr)/2),width=1,ignore.strand=T)
}
getRegionMeth <- function(query,subject,thr=2,verbose=TRUE){
    if (verbose) cat("converting data to GRanges if not supplied as them\n")
    if (class(query)[1] != "GRanges") query=GRanges(query)
    if (class(subject)[1] != "GRanges") subject=GRanges(subject)
    if (verbose) cat("finding overlaps\n")
    ovl=findOverlaps(query,subject)
    freq.tb=tibble(freq=query$freq[queryHits(ovl)],
                   cov=query$cov[queryHits(ovl)],
                   feature.index=subjectHits(ovl))
    if (verbose) cat("calculating methylation by region\n")
    freq.tb[which(freq.tb$cov>=thr),] %>% group_by(feature.index)%>%
        summarize(totcov=sum(cov),
                  numsites=n(),
                  freq=mean(freq))
}

    
getDistance <- function(query,subject){
    # input subject to this must be center of features (done with getCenter)
    cat("getting distances\n")
    if (class(query)[1] != "GRanges") query=GRanges(query)
    if (class(subject)[1] != "GRanges") subject=Granges(subject)
    nearest.idx=nearest(query,subject,ignore.strand=T)
    sub.hit=subject[nearest.idx]
    dist.tb=tibble(index.feature=nearest.idx,
                   chrom.feature=as.character(seqnames(sub.hit)),
                   start.feature=start(sub.hit),
                   end.feature=end(sub.hit),
                   strand.feature=as.character(strand(sub.hit)),
                   dist=ifelse(strand.feature=="-",
                               end.feature-start(query),
                               start(query)-start.feature))
    dist.tb
}
calculate_rollmean <- function(dat.tb,rollrange,win){
    rollmeans=numeric()
    for ( center in rollrange ){
        dat.win=dat.tb[which(dat.tb$dist >= center-win &
                             dat.tb$dist <= center+win),]
        rollmean=mean(dat.win$freq)
#        if (is.nan(rollmean)) rollmean=rollmeans[length(rollmeans)]
#        rollmean=sum(dat.win$totmeth)/
#            (sum(dat.win$totmeth)+sum(dat.win$totunmeth))
        rollmeans=c(rollmeans,rollmean)
    }
    rollmeans
}
aggregate_methylation <- function(dat.dist,win=50){
    cat("aggregating and calculating rolling mean\n")
    if (! "freq" %in% names(dat.dist)) {
        dat.dist = dat.dist %>% 
          mutate(freq=meth/(meth+unmeth)) %>%
          na.omit()
    }
    dat.ag=dat.dist%>%
        group_by(dist)%>%
            summarize(freq=mean(freq)) %>%
        na.omit()
    roll.range=seq(from=min(dat.ag$dist)+win/2,
                   to=max(dat.ag$dist)-win/2,by=1)
    rollmeans=calculate_rollmean(dat.ag,rollrange=roll.range,win=win/2)
    tibble(dist=roll.range,freq=rollmeans)
}
