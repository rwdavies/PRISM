#' @title Convert fasta into RData object files
#' @param fasta path to fasta file
#' @param chrList vector of chromosome names to prep
#' @param verbose whether to print messages. verbose 1 = Yes, 0 = No
#' @param nCores number of cores to use
#' @export
prepFasta=function(fasta,chrList,verbose=1,nCores=1)
{
  prefix=""
  if(substr(fasta,nchar(fasta)-2,nchar(fasta))==".gz") prefix="z"
  # description     
  # 1 - get location of different chromosomes
  # 2 - rip into one text file per chromosome
  # 3 - turn into RData file
  #
  # 1 - get location
  #
  if(verbose==1) print(paste("Getting position of chromosomes in fasta, ",date(),sep=""))
  a=system(paste(prefix,"cat ",fasta," | grep -n dna ",sep=""),intern=TRUE)
  e=system(paste(prefix,"cat ",fasta," | wc -l ",sep=""),intern=TRUE)
  # now - reformat into something useable
  m=t(sapply(strsplit(a,":"),function(x) x))
  m2=cbind(m[,1],t(sapply(strsplit(m[,2]," "),function(x) substr(x,2,nchar(x))))[,1])
  m3=rbind(m2,c(as.integer(e)+1,NA))
  #
  # 2 - rip into one file per chromosome
  #
  out=mclapply(chrList,mc.cores=nCores,function(chr)
  {
    if(verbose==1)   print(paste("Converting fasta to RData object - chromosome ",chr,", ",date(),sep=""))
    x=match(chr,m2[,2])
    start=as.integer(m[x,1])+1
    end=as.integer(m[x+1,1])-1
    refX=system(paste(prefix,"cat ",fasta," | awk ",shQuote(paste("{if((NR>=",start,") && (NR<=",end,")) {print $0}}",sep="")),"",sep=""),intern=TRUE)
    # turn into -1=N, 0=A, 1=C, 2=G, 3=T format
    refX2=unlist(strsplit(refX,""))
    ref=match(refX2,c("N","A","C","G","T"))-2
    s=sum(is.na(ref))
    if(s>0) print(paste("While making the fasta for chromosome ",chr,", ",s," unrecognized bases were converted to N",sep=""))
    ref[is.na(ref)==TRUE]=-1
    save(ref,file=paste(fasta,".",chr,".RData",sep=""))
  })
  if(verbose==1)   print(paste("Done converting fasta to RData objects, ",date(),sep=""))
  return(NULL)
}



#' @title Prep repeat masker file
#' @param rmsk Path to repeat masker file
#' @param simpleRepeats Path to simple repeats file
#' @param chrList vector of chromosome names to prep
#' @param verbose whether to print messages. verbose 1 = Yes, 0 = No
#' @param nCores number of cores to us
#' @param useChr Whether to add chr to repeat masker genoName list
#' @export
prepRepeat=function(rmsk,simpleRepeats,chrList,verbose=1,nCores=1,useChr=TRUE)
{
  # description     
  # 1 - load all results
  # 2 - rip into one RData file per chromosome
  #
  # 1 - load in
  #
  if(verbose==1) print(paste("Loading files, ",date(),sep=""))
  data=read.table(rmsk,header=TRUE,comment.char="@")
  dataSR=read.table(simpleRepeats,header=TRUE,comment.char="@")  
  #
  # 2 - rip into one file per chromosome
  #
  if(verbose==1) print(paste("Save repeat files per chromosome, ",date(),sep=""))
  out=mclapply(chrList,mc.cores=nCores,function(chr)
  {
    ### repeat mask
    if(useChr==TRUE) mask=data[data[,"genoName"]==paste("chr",chr,sep=""),] 
    if(useChr==FALSE) mask=data[data[,"genoName"]==chr,]
    save(mask,file=paste(rmsk,".",chr,".RData",sep=""))
    ### simple repeat
    if(useChr==TRUE) rep=dataSR[dataSR[,"chrom"]==paste("chr",chr,sep=""),] 
    if(useChr==FALSE) rep=dataSR[dataSR[,"chrom"]==chr,]
    save(rep,file=paste(simpleRepeats,".",chr,".RData",sep=""))
  })
  if(verbose==1)   print(paste("Done, ",date(),sep=""))
  return(NULL)
}






#' @title Call hotspots from rate files
#' @param prefix for rate files, prefix, ie path to file without chromosome  
#' @param suffix for rate files, suffix, ie path after chromosome
#' @param chrList chromosome list to find hotspots on
#' @param minRate minimum rate to declare hotspot
#' @param minFrac when a hotspot has been called, the continuous region with rate at this fraction of the maximum rate
#' @param minLength remove identified hotspots with a width smaller than this
#' @param maxLength remove hotspots with width smaller than this
#' @param hotspotSmoothingLength smooth recombination rates over this interval
#' @param positionColumn column name for physical position in recombination rate file
#' @param rateColumn column name for local (i.e. not cumulative) recombination rate in recombination rate file
#' @param verbose whether to print messages. 1 = Yes, 0 = No
#' @param nCores number of cores to use
#' @export
callHotspots=function(prefix,suffix,chrList,minRate=10,minFrac=0.8,minLength=1000,maxLength=8000,maxTrim=NA,hotspotSmoothingLength=2000,positionColumn="Position.bp.",rateColumn="Rate.cM.Mb.",verbose=1,nCores=1)
{
  #
  # check column names exist
  #
  data=read.table(paste(prefix,chrList[1],suffix,sep=""),header=TRUE,nrow=10)
  if(is.na(match(positionColumn,colnames(data))))  STOP("position column name misspecified")
  if(is.na(match(rateColumn,colnames(data))))  STOP("rate column name misspecified")  
  #
  # loop over each chromosome
  #
  out=mclapply(chrList,mc.cores=nCores,function(chr) {
    if(verbose>=2) print(paste("Call hotspots for chromosome ",chr,", ",date(),sep=""))
    #
    # load rate 
    #
    data=read.table(paste(prefix,chr,suffix,sep=""),header=TRUE)
    #
    # turn into one vector
    #
    # position in bp, rate in intervals in cM /Mbp
    rL=cbind(as.integer(data[,positionColumn]),as.numeric(data[,rateColumn]))
    #
    # general strategy
    #
    # 1) make rates chromosome wide
    # 2) make smoothed per-bp rates over 1kb on either side
    # 3) find regions of local maxima
    # 4) find hotspot bounds for each
    # 5) for overlapping hotspots, remove wider one
    #
    n=nrow(rL)
    chrLength=max(rL[,1])-1 # maximum defined rates
    x=rL[-1,1]-rL[-n,1]
    y=rL[-n,2]
    z=c(rep(0,rL[1,1]-1),unlist(sapply(1:(n-1),function(i) rep(y[i],x[i]))))
    # smooth 1kbp either side
    z2=cumsum(c(0,z))
    w=hotspotSmoothingLength
    z3= z2[-(1:(w*2+1))] - z2[-(chrLength+1 + (-(2*w):0))]
    z4=array(0,chrLength)
    z4[-c( 1:(w) , (chrLength + (-(w-1):0)) )]=z3
    z4=z4/(2*w+1)
    rS=z4 # smoothed rate - to 4 decimal places
    ### final local maxima above a certain rate
    # three options - increasing, decreasing, flat
    # want to find regions where it goes from increasing to decreasing
    # or, where it does so, broken up by regions of 0s
    # so want both sides lower and above rate
    ## make bed file from differences
    z1=1*as.integer(diff(rS)>0) + (-1)*as.integer(diff(rS)<0)
    z2=z1[-length(z1)]!=z1[-1]
    c=(1:chrLength)[c(TRUE,z2,FALSE)] # first position of change
    changeMat=cbind(c,c(c[-1],chrLength),rS[c],z1[c])
    colnames(changeMat)=c("start","end","rateAtStart","changeType")
    # each entry specifies within that range,
    # either increasing, decreasing, or staying the same
    # range is inclusive of all bases in change
    changeMatOK=changeMat[changeMat[,"rateAtStart"]>minRate,]
    x=changeMatOK
    # merge adjacent regions
    # peak value by definition has to be in that interval
    # go to interval, find maximum rate
    names=c("peakSpot","start","end","maxRate","avRate") # fake hot matrix if things fail
    hot=matrix(NA,ncol=5,nrow=2)
    colnames(hot)=names
    if(is.null(nrow(changeMatOK))==FALSE) # dammit indian population
    {
      toSwitch=changeMatOK[-1,1]==changeMatOK[-nrow(changeMatOK),2]
      toSwitchL=(1:length(toSwitch))[toSwitch]
      for(j in toSwitchL[length(toSwitchL):1])
        changeMatOK[j,2]=changeMatOK[j+1,2]
      changeMatOK=changeMatOK[-(toSwitchL+1),]
      if(sum((changeMatOK[-1,1]-changeMatOK[-nrow(changeMatOK),2])<=0)>0)
        print(paste("WARNING - no hotspots on chromosome ",chr,sep=""))
      # check discontinuous
      changeMatOK=changeMatOK[,-c(3:4)] # not informative anymore
    }
    if(length(changeMatOK)==4) changeMatOK=matrix(changeMatOK[1:2],nrow=1)      
    ### OK - now fine map within these regions
    # choose peak value within each region
    if(nrow(changeMatOK)==0) return(hot)
    hot=t(apply(changeMatOK,1,function(x) {
      v=x[1]:x[2]
      w=rS[v]==max(rS[v])
      w2=floor(median(v[w])) # take one of the top values within the hotspot
      if(is.na(maxTrim)==TRUE)
      {
        # take top value, find out where less than fraction
        y=seq(w2,min(w2+5000,chrLength))
        y2=(rS[y]) > (minFrac * rS[w2[1]])
        w3=which.min(y2) # first one under
        if(w3==1) w3=5000+1
        # do left side
        y=seq(w2,max(1,w2-5000))
        y2=(rS[y]) > (minFrac * rS[w2])
        w4=which.min(y2) # first one under
        if(w4==1) w4=5000+1
        w3=w3-2; w4=w4-2
      } else {
        w4=round(maxTrim/2); w3=round(maxTrim/2)
      }
      # take only values close to center until a
      return(c(w2,w2-w4,w2+w3,rS[w2],mean(rS[seq(w2-w4,w2+w3)])))
    }))
    colnames(hot)=c("peakSpot","start","end","maxRate","avRate")
    #
    # remove overlapping hotspots - keep smaller ones
    #
    if(verbose>=2)     print(paste("Done calling hotspots for chromosome ",chr,", ",date(),sep=""))    
    return(hot)
    # long!
  })
  names(out)=chrList
  #
  #
  # end of each chromosome - now assemble whole thing
  #
  hot=array(0,c(0,6))
  colnames(hot)=c("chr","peakSpot","start","end","maxRate","avRate")
  for(chr in chrList)   hot=rbind(hot,  cbind(rep(chr,nrow(out[[as.character(chr)]])),out[[as.character(chr)]]))
  hot=hot[-1,]
  hot=hot[is.na(hot[,2])==FALSE,]
  # re-order
  hot=hot[,c(1,3,4,5,6,2)]
  colnames(hot)=c("Chr","Pos_start_bp","Pos_end_bp","maxRate","avRate","peakSpot")
  #
  # remove if too big or too small
  #
  if(verbose>=1)   print(paste("Number of hotspots before removing due to size: ",nrow(hot),sep=""))
  width=(as.integer(hot[,3])-as.integer(hot[,2]))
  hot=hot[width>=minLength & width<=maxLength,]
  hot=hot[order(hot[,1],as.integer(hot[,2])),] 
  if(verbose>=1)  print(paste("Number of hotspots after removing due to size: ",nrow(hot),sep=""))
  if(verbose>=1)   print("Quantile distribution of hotspot sizes")
  if(verbose>=1)   print(quantile(width,probs=c(0,0.05,0.1,0.25,0.5,0.75,0.9,0.95,1)))
  #print(paste("hotspotSmoothingLength = ",hotspotSmoothingLength,sep=""))  
  #print(paste("minRate = ",minRate,sep=""))
  #print(paste("minFrac = ",minFrac,sep=""))
  #print(paste("minLength = ",minLength,sep=""))
  #print(paste("maxLength = ",maxLength,sep=""))
  #
  # done! 
  #
  return(hot)
}



check_and_shrink_hot <- function(hotC, ref, chr) {
    keep <- rep(TRUE, nrow(hotC))
    for(i in 1:nrow(hotC)) {
        if (hotC[i, 1] > length(ref)) {
            keep[i] <- FALSE
        }
    }
    if (sum(keep) > 0) {
        warning(paste0("chromosome ", chr, " has ", sum(!keep), " out of ", length(keep), " hotspots outside of the reference genome. This suggests an incorrect use of reference genome for the given hotspots"))
        hotC <- hotC[keep, ]
    }
    hotC
}



#' @title Get enrichment for a set of hotspots
#' @param chrList List of chromosomes to run over
#' @param hot A set of hotspots
#' @param hotExclude Hotspots to exclude
#' @param hotCenterDist Distance from hotspot centre each way to look for enrichment
#' @param K k-mer length to seed over
#' @param nCores How many computer cores to use
#' @param filterSequences Whether to filter hotspots for repeats and simple repeats
#' @param fasta Path to fasta file
#' @param rmsk Path to repeat masker file
#' @param simpleRepeats Path to simple repeats file
#' @param verbose whether to print messages. 2 = frequent, 1 = moderate, 0 = No
#' @export
getEnrichment <- function(
    chrList,
    hot,
    hotExclude = NULL,
    hotspotCenterDist = 200,
    K = 8,
    nCores = 1,
    filterSequences = TRUE,
    fasta = NULL,
    rmsk = NULL,
    simpleRepeats = NULL,
    verbose = 1
) {

    if(verbose>=1)
        print(paste("Perform searching for motif enrichment, ",date(),sep=""))

    ## check columns of hotspots
    if(sum(is.na(match(c("Chr","Pos_start_bp","Pos_end_bp"),colnames(hot))))>0)
        stop("hotspot matrix column names don't include necessary columns")

    if(verbose>=1)
        print(paste("Loading sequences of chromosomes and motifs, ",date(),sep=""))

    ## get sequence of the motifs
    out <- mclapply(
        chrList,
        mc.cores = nCores,
        function(chr) {

        if(verbose>=2)
            print(paste("Load chromosome ",chr,", ",date(),sep=""))
        
        ## load ref
        file <- paste(fasta,".",chr,".RData",sep="")
        if (!file.exists(file)) {
            stop(paste0("Cannot find file:", file))
        }
        load(file)
        ref[ref==-1]=4  ## OK so 4 is bad
        ## mask out repeats
        if(filterSequences==TRUE) {
            load(paste(rmsk,".",chr,".RData",sep=""))
            for(i in 1:nrow(mask)) ref[ (mask[i,"genoStart"]+1):mask[i,"genoEnd"]]=4
            ## mask out simple repeats
            load(paste(simpleRepeats,".",chr,".RData",sep=""))
            for(i in 1:nrow(rep)) ref[ (rep[i,"chromStart"]+1):rep[i,"chromEnd"]]=4
        }
        
        ## mask our previous hotspots
        if(length(hotExclude)>0) {
            coldL=hotExclude[hotExclude[,"Chr"]==chr,]
            if(length(coldL)>0) {
                for(i in 1:nrow(coldL))
                    ref[coldL[i,2]:coldL[i,3]]=4
            }
        }

        
        ## get hotspots
        hotL <- hot[hot[,"Chr"]==chr,]
        hotL2 <- cbind(as.integer(hotL[,2]),as.integer(hotL[,3]))
        hotL2 <- check_and_shrink_hot(hotL2, ref, chr)
        
        seqs=lapply(1:nrow(hotL),function(i) {
            x=hotL[i,]
            return(ref[as.integer(x["Pos_start_bp"]):as.integer(x["Pos_end_bp"])])
        })
        
        ## in hotspot, in middle 200 bases
        hotC=t(apply(hotL[,c("Pos_start_bp","Pos_end_bp")],1,function(x) {
            x=as.integer(x)
            a1=round(mean(x)-hotspotCenterDist)
            a2=round(mean(x)+hotspotCenterDist)
            return(c(max(a1,x[1]),min(a2,x[2])))
        }))
        hotC <- check_and_shrink_hot(hotC, ref, chr)
        
        
        ## get positions of those within or outside 
        inhotCenter <- unique(unlist(apply(hotC,1,function(x) seq(x[1],x[2]))))
        inhot <- unique(unlist(apply(hotL2,1,function(x) seq(x[1],x[2]))))
        inhotEdges <- setdiff(inhot,inhotCenter)
        ## get the sequences involved!
        hotCenterSeq <- ref[inhotCenter]
        hotEdgesSeq <- ref[inhotEdges]        
        coldSeq <- ref[-inhot] ## rest of genome!
        
        ## hash
        hotCenterSeqH <- simonHash2(hotCenterSeq, K, check = TRUE)
        hotEdgesSeqH <- simonHash2(hotEdgesSeq, K, check = TRUE) 
        coldSeqH <- simonHash2(coldSeq, K, check = TRUE)
        
        ## count - note - simonHash2 is 1-based, so remove first entry
        hotCenterSeqC <- increment(y=as.numeric(hotCenterSeqH),yT=as.integer(length(hotCenterSeqH)),xT=as.integer(4^K))[-1]
        hotEdgesSeqC <- increment(y=as.numeric(hotEdgesSeqH),yT=as.integer(length(hotEdgesSeqH)),xT=as.integer(4^K))[-1]
        coldSeqC <- increment(y=as.numeric(coldSeqH),yT=as.integer(length(coldSeqH)),xT=as.integer(4^K))[-1]
        
        ## return counts
        if(verbose>=2) print(paste("Done chromosome ",chr,", ",date(),sep=""))    
        return(
            list(
                hotCenterSeqC = hotCenterSeqC,
                hotEdgesSeqC = hotEdgesSeqC,
                coldSeqC = coldSeqC,
                seqs = seqs
            )
        )
    })
    names(out) <- chrList

    ##
    ## check for errors
    ##
    check_mclapply_OK(out, "an error occured during getEnrichment")
    
    ##
    ## get sequences
    ##
    seqs <- as.list(sum(unlist(lapply(out,function(x) length(x$seqs)))))
    c=1
    for(j in 1:length(out)) {
        n=length(out[[j]]$seqs)
        seqs[c:(c+n-1)]=out[[j]]$seqs
        c=c+n
    }
    ##
    ## sum counts
    ##
    if(verbose>=1)   print(paste("Add results, build matrices, ",date(),sep=""))
    ##
    ##
    hotCenterSeq=array(0,4**K)
    for(i in chrList) hotCenterSeq=hotCenterSeq + out[[as.character(i)]][[1]]
    hotEdgesSeq=array(0,4**K)
    for(i in chrList) hotEdgesSeq=hotEdgesSeq + out[[as.character(i)]][[2]]
    coldSeq=array(0,4**K)
    for(i in chrList) coldSeq=coldSeq + out[[as.character(i)]][[3]]
    ##
    ## get OR and p!
    ##
    if(verbose>=1)   print(paste("Get ORs and p-values, ",date(),sep=""))
    hashLink=buildHashLink(K)  
    ##
    ##
    ## original - hotspots versus rest of genome
    ##
    x=apply(hashLink[,(ncol(hashLink)+1-K):ncol(hashLink)],1,hashF,K=K) + 1 ##re-order?
    compO1=makeCompO2(r=cbind(hotCenterSeq[x]+ hotEdgesSeq[x],coldSeq[x]),which=array(TRUE,4**K),hashLink=hashLink,K=K,num=2)
    out1=getORAndP(compO=compO1,K=K,num=2,nCores=nCores,verbose=FALSE,data=NULL,gc.content=NULL,CG=NULL,version="robbie",comp2=NULL,outFile=NULL,cgte=0,pvalmet="chisq",rgte=0)
    ##
    ## also - hotspots vs edges
    ##
    compO2=makeCompO2(r=cbind(hotCenterSeq[x],hotEdgesSeq[x]),which=array(TRUE,4**K),hashLink=hashLink,K=K,num=2)
    out2=getORAndP(compO=compO2,K=K,num=2,nCores=nCores,verbose=FALSE,data=NULL,gc.content=NULL,CG=NULL,version="robbie",comp2=NULL,outFile=NULL,cgte=0,pvalmet="chisq",rgte=0)
    ##
    ## also - hotspots in two parts vs rest of genome
    ##
    compO3=makeCompO2(r=cbind(hotCenterSeq[x],hotEdgesSeq[x],coldSeq[x]),which=array(TRUE,4**K),hashLink=hashLink,K=K,num=3)
    out3=getORAndP(compO=compO3,K=K,num=3,nCores=nCores,verbose=FALSE,data=NULL,gc.content=NULL,CG=NULL,version="robbie",comp2=NULL,outFile=NULL,cgte=0,pvalmet="chisq",rgte=0)
    ##
    ## return! both matrices
    ##
    compX1=cbind(out1$p[,1],out1$or[,1],compO1)
    compX1=compX1[order(compX1[,1]),]
    compX2=cbind(out2$p[,1],out2$or[,1],compO2)
    compX2=compX2[order(compX2[,1]),]
    compX3=cbind(out3$p,out3$or,compO3)
    compX3=compX3[order(compX3[,1]),]
    if(verbose>=1) print(paste("Done searching for motif enrichment, ",date(),sep=""))
    return(
        list(
            compX1 = compX1,
            compX2 = compX2,
            compX3 = compX3,
            seqs = seqs
        )
    )
}







#' @title Get enrichment for a set of hotspots
#' @param seqsH hotspot sequences
#' @param seqsC coldspot sequences
#' @param hotCenterDist Distance from hotspot centre each way to look for enrichment
#' @param K k-mer length to seed over
#' @param nCores How many computer cores to use
#' @param filterSequences Whether to filter hotspots for repeats and simple repeats
#' @param rmsk Path to repeat masker file
#' @param simpleRepeats Path to simple repeats file
#' @param verbose whether to print messages. 2 = frequent, 1 = moderate, 0 = No
#' @export
getEnrichmentFromSequences <- function(
    seqsH,
    seqsC,
    hotspotCenterDist = 200,
    K = 8,
    verbose = 1,
    nCores = 1
) {

    if(verbose>=1)
        print_message("Perform searching for motif enrichment")


    hotCenterSeq <- array(0,4**K)
    hotEdgesSeq <- array(0,4**K)

    for(seq in seqsH) {
        
        middle <- round(nchar(seq) / 2)
        a1 <- round(middle - hotspotCenterDist)
        a2 <- round(middle + hotspotCenterDist)
        
        hotCenterSeqX <- match(strsplit(substr(seq, a1, a2), "")[[1]], c("A", "C", "G", "T", "N"))
        hotCenterSeqH <- simonHash2(hotCenterSeqX - 1,K) ## needs 0-based, is 1-based
        hotCenterSeq <- hotCenterSeq + increment(
            y = as.numeric(hotCenterSeqH),
            yT = as.integer(length(hotCenterSeqH)),
            xT = as.integer(4^K)
        )[-1]
        
        hotEdgesSeqX <- match(
            c(
                strsplit(substr(seq, 1, a1 - 1), "")[[1]],
                strsplit(substr(seq, a2 + 1, nchar(seq)), "")[[1]]
            ),
            c("A", "C", "G", "T", "N")
        )
        hotEdgesSeqH <- simonHash2(hotEdgesSeqX - 1,K)
        hotEdgesSeq <- hotEdgesSeq + increment(
            y = as.numeric(hotEdgesSeqH),
            yT = as.integer(length(hotEdgesSeqH)),
            xT = as.integer(4^K)
        )[-1]
        
    }

    coldSeq <- array(0, 4 ** K)
    
    for(seq in seqsC) {
        
        middle <- round(nchar(seq) / 2)
        a1 <- round(middle - hotspotCenterDist)
        a2 <- round(middle + hotspotCenterDist)
        
        coldSeqX <- match(strsplit(seq, "")[[1]], c("A", "C", "G", "T", "N"))
        coldSeqH <- simonHash2(coldSeqX - 1,K)
        coldSeq <- coldSeq + increment(
            y = as.numeric(coldSeqH),
            yT = as.integer(length(coldSeqH)),
            xT = as.integer(4^K)
        )[-1]
        
    }
    
    if(verbose>=1)
        print_message("Combine results")

    hashLink <- buildHashLink(K)

    
    ## original - hotspots versus rest of genome
    x <- apply(
        hashLink[,(ncol(hashLink)+1-K):ncol(hashLink)],
        1,
        hashF,
        K=K
    ) + 1 ## re-order (none performed)
    compO1 <- makeCompO2(
        r = cbind(hotCenterSeq[x]+ hotEdgesSeq[x],coldSeq[x]),
        which = array(TRUE,4**K),
        hashLink = hashLink,
        K = K,
        num = 2
    )
    
    out1 <- getORAndP(
        compO = compO1,
        K = K,
        num = 2,
        nCores = nCores,
        verbose = FALSE,
        data = NULL,
        gc.content = NULL,
        CG = NULL,
        version = "robbie",
        comp2 = NULL,
        outFile = NULL,
        cgte = 0,
        pvalmet = "chisq",
        rgte = 0
    )
    compX1 <- cbind(out1$p[,1],out1$or[,1],compO1)
    compX1 <- compX1[order(compX1[,1]),]

    
    ##
    ## also - hotspots vs edges
    ##
    compO2=makeCompO2(r=cbind(hotCenterSeq[x],hotEdgesSeq[x]),which=array(TRUE,4**K),hashLink=hashLink,K=K,num=2)
    out2=getORAndP(compO=compO2,K=K,num=2,nCores=nCores,verbose=FALSE,data=NULL,gc.content=NULL,CG=NULL,version="robbie",comp2=NULL,outFile=NULL,cgte=0,pvalmet="chisq",rgte=0)
    compX2=cbind(out2$p[,1],out2$or[,1],compO2)
    compX2=compX2[order(compX2[,1]),]
    
    ##
    ## also - hotspots in two parts vs rest of genome
    ##
    compO3=makeCompO2(r=cbind(hotCenterSeq[x],hotEdgesSeq[x],coldSeq[x]),which=array(TRUE,4**K),hashLink=hashLink,K=K,num=3)
    out3=getORAndP(compO=compO3,K=K,num=3,nCores=nCores,verbose=FALSE,data=NULL,gc.content=NULL,CG=NULL,version="robbie",comp2=NULL,outFile=NULL,cgte=0,pvalmet="chisq",rgte=0)
    compX3=cbind(out3$p,out3$or,compO3)
    compX3=compX3[order(compX3[,1]),]
    
    
    if(verbose>=1)
        print_message("Done searching for motif enrichment")

    ## seqs here 1-4 numbers
    seqs <- lapply(seqsH, function(seq) -1 + match(strsplit(seq, "")[[1]], c("A", "C", "G", "T")))
    
    return(
        list(
            compX1 = compX1,
            compX2 = compX2,
            compX3 = compX3,
            seqs = seqs
        )
    )
}


































#' @title After calculating enrichment, initialize a PWM
#' @param motifEnrichment Specific input from getEnrichment
#' @param K k-mer length to consider
#' @param initializationMethod Either "simple" or "advanced". Simple just takes top k-mer, advanced builds simple PWM
#' @param eligible Which k-mers are eligible to be considered
#' @param verbose Verbosity. 1 = Full
#' @export
initializePWM=function(motifEnrichment,K=8,initializationMethod="simple",eligible=NULL,verbose=1,pwmValue=0.9)
{
  if(verbose>=1) print(paste("Initialize PWM, ",date(),sep=""))
  # these aren't sorted
  compX2=motifEnrichment$compX2
  compX3=motifEnrichment$compX3
  # sort together
  compX2=compX2[order(compX2[,"F"]),]
  compX3=compX3[order(compX3[,"F"]),]
  # build nucleotide diversity on the fly
  a=sapply(0:3,function(i) {rowSums(compX2[,ncol(compX2)+1+(-K):(-1)]==i)})
  nd=rowSums(a>0)
  # check
  if(sum(compX2[,"F"]!=compX3[,"F"])>0) print("WARNING - FAILED ALIGNMENT") 
  orm=compX2[compX2[,2]>1 & is.na(compX2[,2])==FALSE & nd!=1 ,"F"]
  compX3G=compX3[match(orm,compX3[,"F"]),]
  compX3G=compX3G[order(compX3G[,1]),]
  pwm=compX3G[1,ncol(compX3G)+(-K+1):0]
  dimvec=length(pwm)
  pwmX=array((1-pwmValue)/3,c(length(pwm),4))
  pwmX[cbind(1:length(pwm),pwm+1)]=pwmValue
  #c(0,2,0,3,2,1,3,0,1,3)
  scorematset=log(pwmX)
  findWordsOut=NULL
  if(initializationMethod=="advanced")
  {
    # consider a further seeder here
    pThresh=0.05/((K-1)**2+K*4)
    pwmMinimum=0.1
    plineage=compX3[,1] # actually needed
    orlineage=compX3[,4] # only for reporting
    plineage[is.na(plineage)]=1
    orlineage[is.na(orlineage)]=1
    compOX=compX3 # used for forward, reverse
    # use compX2 and compX3 from above
    eligibleX=compX2[,2]>1 & compX3[,4]> 1 & is.na(compX2[,2])==FALSE & nd & plineage<pThresh & nd>1
    if(is.null(eligible)==TRUE) eligible=eligibleX
    if(is.null(eligible)==FALSE) eligible=eligibleX & eligible
    # set p and OR
    whichPWM=which.min(compX3[eligible,1])
    whichPWM=match(compX3[eligible,"F"][whichPWM],compX3[,"F"])
    pwm=compX3[whichPWM,ncol(compX3)+(-K+1):0]
    findWordsOut=findWords(motifSeed=pwm,eligible=eligible,plineage=plineage,orlineage=orlineage,K=K,compOX=compOX,pThresh=pThresh,pwmMinimum=pwmMinimum)
    scorematset=log(t(findWordsOut$sum))
    dimvec=nrow(scorematset)
    # remove it from eligibility of course
    eligible[whichPWM]=FALSE
    pwm=NULL
 }
  if(verbose>=1) print(paste("Done initializing PWM, ",date(),sep=""))
  return(list(pwm=pwm,dimvec=dimvec,scorematset=scorematset,findWordsOut=findWordsOut))
}






#' @title Get sequences for a set of hotspots
#' @param chrList List of chromosomes to run over
#' @param hot A set of hotspots
#' @param nCores How many computer cores to use
#' @param filterSequences Whether to filter hotspots for repeats and simple repeats
#' @param rmsk Path to repeat masker file
#' @param simpleRepeats Path to simple repeats file
#' @param verbose whether to print messages. 2 = frequent, 1 = moderate, 0 = No
#' @export
getSequences=function(chrList,hot,nCores=1,filterSequences=TRUE,rmsk=NULL,simpleRepeats=NULL,verbose=1)
{
  #
  if(verbose>=1) print(paste("Get original sequences, ",date(),sep=""))
  #
  # check columns of hotspots
  if(sum(is.na(match(c("Chr","Pos_start_bp","Pos_end_bp"),colnames(hot))))>0)
    stop("hotspot matrix column names don't include necessary columns")
  #
  #
  if(verbose>=1) print(paste("Loading sequences of chromosomes and motifs, ",date(),sep=""))
  #
  # get sequence of the motifs
  out=mclapply(chrList,mc.cores=nCores,function(chr) {
    #
    if(verbose>=2) print(paste("Load chromosome ",chr,", ",date(),sep=""))
    # load ref
    load(paste(fasta,".",chr,".RData",sep=""))
    ref[ref==-1]=4 
    # mask out repeats
    if(filterSequences==TRUE)
    {
      load(paste(rmsk,".",chr,".RData",sep=""))
      for(i in 1:nrow(mask)) ref[ (mask[i,"genoStart"]+1):mask[i,"genoEnd"]]=4
      # mask out simple repeats
      load(paste(simpleRepeats,".",chr,".RData",sep=""))
      for(i in 1:nrow(rep)) ref[ (rep[i,"chromStart"]+1):rep[i,"chromEnd"]]=4
    }
    #
    # get hotspots
    #
    hotL=hot[hot[,"Chr"]==chr,]
    seqs=lapply(1:nrow(hotL),function(i) {
      x=hotL[i,]
      return(ref[as.integer(x["Pos_start_bp"]):as.integer(x["Pos_end_bp"])])
    })
    return(list(seqs=seqs,x=NA))
  })
  a=sapply(out,function(x) length(x))
  if(sum(a!=2)>0) {
    print("WARNING - getSequences failed")
    print(out[[which.min(a<2)]])
  }
  names(out)=chrList
  #
  # get sequences
  #
  seqs=as.list(sum(unlist(lapply(out,function(x) length(x$seqs)))))
  c=1
  for(j in 1:length(out))
  {
    n=length(out[[j]]$seqs)
    seqs[c:(c+n-1)]=out[[j]]$seqs
    c=c+n
  }
  #
  # return! both matrices
  #
  if(verbose>=1) print(paste("Done getting original sequences, ",date(),sep=""))
  return(list(seqs=seqs))
}



