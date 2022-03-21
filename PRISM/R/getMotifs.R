#' @title Refine motif in hotspots from seed sequence
#' @param hot A set of hotspots
#' @param seqsH A list of hotspots (use initializationMethod advanced, dataType sequence)
#' @param seqsC A list of coldspots (use initializationMethod advanced, dataType sequence)
#' @param hotExclude Hotspots to exclude (for instance, if these are known to contain matches to another motif)
#' @param eligible k-mers already used as seeds (or NULL)
#' @param nCores Number of computer cores to use
#' @param K Length of k-mers to consider for seeding
#' @param rmsk Path to repeat masker file (assumes preparation already performed)
#' @param simpleRepeats Path to simple repeats file (assumes preparation already performed)
#' @param initializationMethod Whether to initialize using a single top k-mer (simple), by constructing a rudimentary PWM (advanced), or to use pre-existing motif probabilities (existing)
#' @param dataType Whether input data is sequences (sequence) or chromosome and fasta based 
#' @param pwmValue If initializing simply, the probability the PWM contains the sequence of the k-mer
#' @param alpha To complete
#' @param incprob To complete
#' @param maxits To complete
#' @param plen To complete
#' @param updatemot To complete
#' @param updatealpha To complete
#' @param ourprior To complete
#' @param updateprior To complete
#' @param bg To complete
#' @param plotting Whether to print an output PDF ("pdf"), to screen ("screen"), or not (anything else)
#' @param plotPrefix If outputting PDFs, what to call them
#' @param verbose Verbosity. 1 = Standard, 2 = Detailed, otherwise no printing
#' @param scorematset To complete
#' @param dimvec To complete
#' @param pwm When using your own initialization, 
#' @export
getMotifs <- function(
    hot,
    seqsH,
    seqsC,
    hotExclude=NULL,
    eligible=NULL,
    nCores=1,
    K=8,
    hotspotCenterDist=200,
    filterSequences=FALSE,
    rmsk=NULL,
    simpleRepeats=NULL,
    discovery=TRUE,
    initializationMethod="advanced",
    dataType="normal",
    pwmValue=0.9,
    alpha=0.5,
    incprob=0.99999,
    maxits=30,
    plen=0.05,
    updatemot=1,
    updatealpha=1,
    ourprior=NULL,
    updateprior=1,
    bg=-1,
    plotting="pdf",
    plotPrefix="getMotifs.plot",
    verbose=1,
    scorematset=NULL,
    dimvec=NULL,
    pwm=NULL
) {

    ## scorematset=NULL,pwm=NULL ,seqs,dimvec=NULL
  
    ## backmethod - which background method to use
    ## 1 = original tri method
    ## 2 = do based on scores

    ## NOTE - assume sequences already sanity checked for maximum and minimum lengths
    if(verbose>=1) {
        print(paste(" ------------------------------------------------------- ",sep=""))
        print(paste(" ---- Program Initializing ",date()," -----",sep=""))
        print(paste(" ------------------------------------------------------- ",sep=""))
    }


    

    ## Whether to initialize using a single top k-mer (simple), by constructing a rudimentary PWM (advanced), or to use pre-existing motif probabilities (existing)
    if(initializationMethod!="existing") {
        ##
        ## step -1 - do motif enrichment
        ##
        if (dataType != "sequence") {
            motifEnrichment <- getEnrichment(
                chrList = chrList,
                hot = hot,
                hotExclude = hotExclude,
                nCores = nCores,
                K = K,
                hotspotCenterDist = hotspotCenterDist,
                filterSequences = TRUE,
                rmsk = rmsk,
                simpleRepeats = simpleRepeats,
                verbose = verbose,
                dataType = dataType
            )
        } else {
            motifEnrichment <- getEnrichmentFromSequences(
                seqsH = seqsH,
                seqsC = seqsC,
                K = K,
                hotspotCenterDist = hotspotCenterDist,
                verbose = verbose,
                nCores = nCores
            )
        }
        ##
        ## step 0 - do motif enrichment
        ##
        initialData <- initializePWM(
            motifEnrichment = motifEnrichment,
            K = K,
            initializationMethod = initializationMethod,
            verbose = verbose,
            eligible = eligible
        )
        pwm <- initialData$pwm
        dimvec <- initialData$dimvec
        scorematset <- initialData$scorematset
        seqsX <- motifEnrichment$seqs
        ##
        ##
        ##
    } else if (initializationMethod == "existing") {
        ##
        ##
        ##
        seqsX <- getSequences(
            chrList = chrList,
            hot = hot,
            nCores = nCores,
            filterSequences = filterSequences,
            rmsk = rmsk,
            simpleRepeats = simpleRepeats,
            verbose = verbose
        )$seqs
        ## if a PWM is supplied, build scorematset
        if(is.null(pwm)==FALSE & is.null(scorematset)==TRUE) {
            pwmX=array((1-pwmValue)/3,c(length(pwm),4))
            pwmX[cbind(1:length(pwm),pwm+1)]=pwmValue
            ##c(0,2,0,3,2,1,3,0,1,3)
            scorematset=log(pwmX)
            dimvec=length(pwm)
        }
        ## otherwise, go with scorematset, which is already defined
    } else {
        stop("Bad initialization approach")
    }



    ##
    ## reconfigure some inputs
    ##
    ##if(input=="seq") seqsL=lapply(seqs,function(x) unlist(strsplit(x,"")))
    ##if(input=="list") {
    seqs=lapply(seqsX,function(x) c("A","C","G","T","N")[x+1])
    seqsL=seqs; seqs=sapply(seqs,function(x) paste(x,collapse=""))
    
    origseqs=seqs
    maxwidth=max(sapply(seqs,length))
    ##if(is.na(minLength)==FALSE) sequencesToKeep=sapply(seqsL,length)>max(len,minLength)
    ##seqsL=seqsL[sequencesToKeep]
    ##seqs=seqs[sequencesToKeep]


    ##
    ## null inputs
    ##
    bgstore=NULL
    qvecstore=NULL
    bg2=NULL


  


    ##
    ## new variables
    ##
    nSequences=length(seqsL)
    originalSequenceLength=sapply(seqsL,length)
    maxwidth=max(originalSequenceLength)



    ##
    ## initialize and some more new variables
    ##
    alphas=matrix(nrow=0,ncol=length(dimvec))
    if(length(alpha)!=length(dimvec)){
	print("Initialising fractions as uniform")
	alpha=rep(1/length(dimvec),length(dimvec))/2
    }
    logfact=cumsum(log(1:(length(seqs)+4)))
    logfact=c(0,logfact)


    ##
    ## simon stuff
    ##
    qvec=vector(length=0)
    fullseqs=seqs
    seqs=gsub("A","1",seqs)
    seqs=gsub("C","2",seqs)
    seqs=gsub("G","3",seqs)
    seqs=gsub("T","4",seqs)
    seqs=gsub("N","5",seqs)
    temp=rep("5",maxwidth)
    temp=paste(temp,collapse="")
    for(i in 1:length(seqs)) if(originalSequenceLength[i]<maxwidth) seqs[i]=paste(seqs[i],substring(temp,1,(maxwidth-originalSequenceLength[i])),sep="")

    ##
    ##index each sequence by its triplet
    ##
    newmat=array(5,c(nSequences,maxwidth))
    for(i in 1:nSequences) { 
        x=match(seqsL[[i]],c("A","C","G","T","N"))
        newmat[i,1:length(x)]=x
    }
    newmatOriginal=newmat
    index=array(1,dim(newmat)) #newmat # first two cols dont make sense
    index[,3:ncol(newmat)]=16*(newmat[,1:(ncol(newmat)-2)]-1)+4*(newmat[,2:(ncol(newmat)-1)]-1)+newmat[,3:ncol(newmat)] # OK
                                        # if it is a 5 - block all 3 possible entries
    index[newmat==5]=65
    index[,2:ncol(newmat)][newmat[,1:(ncol(newmat)-1)]==5]=65
    index[,3:ncol(newmat)][newmat[,1:(ncol(newmat)-2)]==5]=65

    ##
    ## prior probabilities
    ##
    if(is.null(ourprior))
    {
        prior=rep(0.1,10) # 10 regions
    } else {
        prior=ourprior
        prior=prior/sum(prior)
    }






    ##
    ## begin iterations
    ##


    its=0;
    results=as.list(1:maxits)

    while(its<=maxits) {
        
        its=its+1
        if(verbose>=1)   print(paste("-------- Start it=",its," -------- " ,date(),sep=""))



        ##
        ## setup! 
        ##
        if(verbose>=2)  print(paste("Setting up - ",date()))
                                        #
        ## get start and end positions of each motif within the score matrix.
        ## note need to do this each time, as size of things can change
                                        #
        starts=c(1,cumsum(dimvec)+1)
        starts=starts[1:length(dimvec)]
        ends=cumsum(dimvec)
        motifLengths=ends-starts+1
                                        #print(paste("starts - ",starts,collapse=","))
                                        #print(paste("ends - ",ends,collapse=","))
        if(verbose>=2)   print(paste("dimvec - ",dimvec,collapse=","))    
        ##score
        bestpos=matrix(0,nrow=length(seqs),ncol=length(starts))
        beststrand=bestpos
        newmat=matrix(5,ncol=maxwidth,nrow=length(seqs))
        scores=newmat
        scores2=newmat
        if(its==1 | updatemot==1)
        { 
            overallscores=matrix(0,nrow=length(seqs)*length(starts),ncol=maxwidth)
            overallscores2=overallscores
            background=overallscores
            priormat=overallscores
        }
        posmat=matrix(0,nrow=length(seqs)*length(starts),ncol=maxwidth)
        for(i in 1:maxwidth) posmat[,i]=i
        newmat=newmatOriginal
        
        
        ##scoremat now just prob of each base given binding
        ##we'll get probs of a sequence without binding in the matrix "background"
#####get sets of probabilities


        ##
        ## start some probability initialization
        ##
        if(verbose>=2)   print(paste("Beginning scoring - ",date()))
        for(j in 1:length(starts)) {
            if(verbose>=2)     print(paste("Scoring motif ",j," - ",date())) 
                                        #
            scoremat=scorematset[starts[j]:ends[j],]
            scoremat=cbind(scoremat,rep(-10000,nrow(scoremat))) #nick removed pointless maxwidth=200 setting
                                        #
            compmat=scoremat[,c(4:1,5)]
            compmat=compmat[nrow(compmat):1,]
            if(updatemot==1 | its==1)
            {
                startrange=seq(1,maxwidth-nrow(scoremat)+1)
                maxes=max(startrange)-1 # this is a value?
                                        # get scores
                scores=addMatMat(newmat-1,scoremat)
                scores2=addMatMat(newmat-1,compmat)
                                        #
                                        # original calculations
                                        #
                overall=scores
                overall[scores2>scores]=scores2[scores2>scores]
                for(i in 1:nrow(overall))
                {
                    bestpos[i,j]=sample(c(which(overall[i,]==max(overall[i,])),which(overall[i,]==max(overall[i,]))),1)
                    score1=scores[i,bestpos[i,j]]
                    score2=scores2[i,bestpos[i,j]]
                    beststrand[i,j]=1
                    if(score2>score1) beststrand[i,j]=0 # I'm not sure if this is used? 
                }
                                        #
                                        # compare against genomic background
                                        #
                if(1==0) # other background method - tried - can re-instane
                {
                    print(paste("Background scoring - ",date()))
                                        # note - no need to compare to background later
                    out=getProbabilityUnderGenomeNull(chrs=18:19,tempdir=tempdir,nCores=nCores,scoremat=scoremat,compmat=compmat,nAlternative=nAlternative)
                    oA=out$oA # distribution under alternative
                    oB=out$oB # distribution under background
                    fileName=paste(outputdir,"plots/",plotPrefix,".dist.",its,".forward.png",sep="")
                    scores=log(getNormalize(s=scores,oA=oA,oB=oB,naVal=1e-9,plot=TRUE,fileName=fileName,its=its))
                                        #fileName=paste(outputdir,"plots/",plotPrefix,".dist.",its,".reverse.png",sep="")
                    scores2=log(getNormalize(s=scores2,oA=oA,oB=oB,naVal=1e-9,plot=FALSE,fileName=fileName,its=its))
                                        # note - be careful - to switch binning
                                        # see both getNormalize and normalize
                }
                overallscores[(nrow(newmat)*(j-1)+1):(nrow(newmat)*j),]=-1e9
                overallscores[(nrow(newmat)*(j-1)+1):(nrow(newmat)*j),1:ncol(scores)]=scores
                overallscores2[(nrow(newmat)*(j-1)+1):(nrow(newmat)*j),]=-1e9
                overallscores2[(nrow(newmat)*(j-1)+1):(nrow(newmat)*j),1:ncol(scores)]=scores2
            } # 
        } # closes for(j in 1:length(starts)
                                        #if(updatemot!=1 & its!=1 & backMethod==1)
                                        #{
        ##after first iteration, if not changing motif then unsubtract background from last iteration
        overallscores=overallscores+background
        overallscores2=overallscores2+background




        ##
        ## set up background
        ##
        if(verbose>=2)   print(paste("Set up background - ",date(),sep=""))
        out=setUpBackground(bg=bg,qvec=qvec,originalSequenceLength=originalSequenceLength,newmat=newmat,bgstore=bgstore,qvecstore=qvecstore,bg2=bg2)
        bg=out$bg
        bgstore=out$bgstore
        qvec=out$qvec
        qvecstore=out$qvec
        bg2=out$bg2
        

        ##
        ## make a background matrix for the different motifs, based on newqvec
        ##
                                        # backMethod==1 - original background method
        if(verbose>=2)   print(paste("Putting in background information - ",date(),sep=""))
                                        #out=makeBackground2(starts=starts,ends=ends,qvec=qvec,newqvec=newqvec,index=index,nSequences=nSequences,maxwidth=maxwidth)
        out=makeBackground2(starts=starts,ends=ends,qvec=qvec,newqvec=newqvec,index=index,nSequences=nSequences,maxwidth=maxwidth)
        background=out$background
        overallscores=overallscores-background
        overallscores2=overallscores2-background


        
        
        





        ##
        ## get prior mat
        ##
                                        # priormat is a matrix of nrow = nSequences, 1 entry per site, equal to the prior probability of it being in that region
        if(verbose>=2)   print(paste("Get prior matrix - ",date(),sep=""))
        out=getPriorMat(nSequences=nSequences,starts=starts,originalSequenceLength=originalSequenceLength,motifLengths=motifLengths,posmat=posmat,prior=prior,alpha=alpha)
        endpos=out$endpos
        priormat=out$priormat


        

        


        
        ##
        ## sampling probabilities for motifs
        ##
        if(verbose>=2) print(paste("Make sampling matrix - ",date(),sep=""))
                                        #
                                        # samplemat has rows equal to number of sequences, 1 set of columns for each 
                                        # motif. later, use cumulative sums to determine 
                                        #
        scoremat=scorematset
        postforward=exp(overallscores) * priormat
        postbackward=exp(overallscores2) * priormat
        samplemat=matrix(nrow=nrow(newmat),ncol=maxwidth*length(starts)*2)
        s=maxwidth
        for(j in 1:length(starts))
        {
            samplemat[,(s*(j-1)*2+1):(s*(2*j-1))]=postforward[(nSequences*(j-1)+1):(nSequences*j),];
            samplemat[,(s*(2*j-1)+1):(s*(2*j))]=postbackward[(nSequences*(j-1)+1):(nSequences*j),];
        } 
        samplemat=samplemat/(rowSums(samplemat)+(1-sum(alpha)))

        ##
        ## sampling! 
        ##
        if(verbose>=2)   print(paste("Determine state probabilities - ",date(),sep=""))
        ##need to record which motif if any
        ##which strand
        ##best position,  best strand
        ##start point 
        regprobs=matrix(0,nrow=nSequences,ncol=(length(starts)))
        s=ncol(postforward)
        for(j in 1:length(starts))
            regprobs[,j]=rowSums(samplemat[,(s*(j-1)*2+1):(s*(2*j))])
                                        # testing of the method - not needed
        regprob=rowSums(regprobs)

        

        


        

        
        ##
        ## sampling! 
        ##
        if(verbose>=2)   print(paste("Sample motif locations and whether bound - ",date(),sep=""))
        mot=array(0,nSequences)
        testmat=samplemat
        for(j in 2:ncol(samplemat))
            testmat[,j]=testmat[,(j-1)]+testmat[,j]
        whichcol=array(0,nSequences)
        q=runif(nSequences)
        for(j in ncol(samplemat):1)    whichcol[q<=testmat[,j]]=j # wow, interesting
        mot[whichcol!=0]=1
        whichmot=floor((whichcol-1)/2/ncol(overallscores))+1
        whichstrand=(floor((whichcol-1)/ncol(overallscores))+1)%%2
        whichpos=whichcol-(whichmot-1)*ncol(overallscores)*2-(1-whichstrand)*ncol(overallscores)
        whichmot=whichmot[mot==1]
        whichpos=whichpos[mot==1]
        whichstrand=whichstrand[mot==1]
        ##



        ##pdf("~/tmp.pdf",height=5,width=20)
        ##plot(samplemat[30,])
        ##dev.off()
        ##system("rsync -av ~/tmp.pdf dense:~/tmp.pdf")
        ##seqsL[[30]][  2001-whichpos[30] +-10:10]
        

        
#####get a prior on positions
####update prior using sampled positions
        totals=1:length(starts)
        for(i in 1:length(totals)) totals[i]=sum(whichmot==i)
        totals=c(totals,sum(mot==0))
###make sure prior prob of motif is 1/2 in dirichlet
        totals[length(totals)]=totals[length(totals)]+length(starts)-1
####need mcmc pack for rdirichlet(n,alpha)
        alphanew=rdirichlet(1,alpha=1+totals)
        ##alphanew=rbeta(1,shape1=sum(mot==1)+1,shape2=sum(mot==0)+1)
####sample start pos
        alphanew=alphanew[1:length(starts)]
        if(verbose>=2)   print("Alpha values sampled:")
        if(verbose>=2)   print(c(alphanew,sum(alphanew)))
        whichregs=which(mot==1)
        v=hist((whichpos-1)/(originalSequenceLength[whichregs]-dimvec[whichmot]),breaks=seq(0,1,0.1),plot=FALSE)
        saveForHist=(whichpos-1)/(originalSequenceLength[whichregs]-dimvec[whichmot])
        if(updateprior==1)
        {
            prior=v$counts+5
            prior=prior/sum(prior)
        }
                                        # fancier saveForHist - for multiple samples
        saveForHistL=lapply(1:length(starts),function(iMotif) {
            return((whichpos[whichmot==iMotif]-1)/(originalSequenceLength[whichregs[whichmot==iMotif]]-dimvec[iMotif]))
        })



        ##
        ##
        ##
        strand=whichstrand
        if(verbose>=2)   print(paste("Sampling sequences - ",date(),sep=""))
        bg=bgstore
        qvec=qvecstore
        ourcounts=matrix(nrow=65,ncol=0)
        oursummary=matrix(nrow=5,ncol=0)
###new, stronger background model
        newbackground=matrix(0,nrow=length(starts),ncol=64)
        for(j in 1:length(starts)){
            scoremat=scorematset[starts[j]:ends[j],]
            tempregs=whichregs[whichmot==j]
            newbackground[j,]=colSums(bg2[whichregs[whichmot==j],])
            ##  debug
            tempstarts=whichpos[whichmot==j]
            tempends=tempstarts+nrow(scoremat)-1
            tempstrand=whichstrand[whichmot==j]
            ##
            ourseqs=matrix(nrow=length(tempregs),ncol=nrow(scoremat)+50)
            sampleseqs=seqs[tempregs]
            v=substring(sampleseqs,tempstarts-25,tempends+25)
            for(i in 1:length(v)) if(tempstarts[i]<=25){
                                      v[i]=paste(c(rep("5",25-tempstarts[i]+1),v[i]),collapse="")
                                  }
                                        # robbie - small fix
            xxx=nchar(sampleseqs)
            for(i in 1:length(v)) if(tempends[i]+25>xxx[i]){
                                      v[i]=paste(c(v[i],rep("5",tempends[i]+25-nchar(sampleseqs)[i])),collapse="")
                                  }
            ##  have subsequences
            for(k in 1:(nrow(scoremat)+50))
                ourseqs[,k]=as.double(substring(v,k,k))
            ## generate the complement
            ourseqs2=ourseqs[,ncol(ourseqs):1]
            ourseqs2=5-ourseqs2
            ourseqs2[ourseqs2==0]=5
            ## end of generating complement
            ourseqs[tempstrand==0,]=ourseqs2[tempstrand==0,] # flip! 
            summary=matrix(nrow=4,ncol=ncol(ourseqs))
            for(i in 1:4) summary[i,]=colSums(ourseqs==i)
            oursummary=cbind(oursummary,rbind(summary,j))
            ## this gives us counts - now need to get likelihood under a background model (more work!)
            ## enables us to sample a new motif in a relatively, though not completely, "principled" manner
            newtemp=16*(ourseqs[,1:(ncol(ourseqs)-2)]-1)+4*(ourseqs[,2:(ncol(ourseqs)-1)]-1)+(ourseqs[,3:(ncol(ourseqs)-0)]-1)
            newtemp[ourseqs[,1:(ncol(ourseqs)-2)]==5]=-1
            newtemp[ourseqs[,2:(ncol(ourseqs)-1)]==5]=-1
            newtemp[ourseqs[,3:(ncol(ourseqs)-0)]==5]=-1
            newtemp=newtemp+1
            motcounts=matrix(ncol=ncol(newtemp),nrow=64)
            for(i in 1:64) motcounts[i,]=colSums(newtemp==i)
            ##newbackground[j,]=newbackground[j,]-rowSums(motcounts)
            bg=bg-rowSums(motcounts)
            ourcounts=cbind(ourcounts,rbind(motcounts,j))	
###note can still "fold over" later
        }

        
        ## manually inspect alignment
        ##pdf(paste(plotPrefix,"ourseqs.",its,".pdf",sep=""),height=20,width=5)
        ##n1=nrow(ourseqs2)
        ##n2=ncol(ourseqs2)
        ##plot(x=0,y=0,ylim=c(0,n1),xlim=c(0,n2))
        ##colStore=c("blue","red","green","black","grey")
        ##for(j in 1:n2)
        ##  rect(ybottom=0:(n1-1),ytop=1:(n1),xleft=j-1,xright=j,col=colStore[ourseqs[,j]],border=NA)
        ##dev.off()
        ##system("rsync -av ~/tmp.pdf dense:~/tmp.pdf")
        

        ##
        if (verbose >= 2) {
            print_message("More background stuff")
        }
        e=c("A","C","G","T")
        for(i in 1:2) e=c(paste("A",e,sep=""),paste("C",e,sep=""),paste("G",e,sep=""),paste("T",e,sep=""))
        d=c("T","G","C","A")
        for(i in 1:2) d=c(paste(d,"T",sep=""),paste(d,"G",sep=""),paste(d,"C",sep=""),paste(d,"A",sep=""))
        names(bg)=e
        bg=bg+bg[d]
        ##add one for prior
        bg=bg+1
        colnames(newbackground)=e
        newbackground=newbackground+newbackground[,d]
        qvec=as.vector(bg)/sum(bg)
        newqvec=newbackground/rowSums(newbackground)


        
        ##
        ## get some vectors
        ##
        if(verbose>=2) print(paste("Get prediction vectors - ",date(),sep=""))    
        out=getPredFrac(
            starts = starts,
            qvec = qvec,
            newqvec = newqvec
        )
        predfrac=out$predfrac
        newpredfrac=out$newpredfrac  



        

        ##
        ## 
        ##
        if(verbose>=2)   print(paste("Summarizing - ",date(),sep=""))      
        newnewdim=vector(length=0)
        newmatset=matrix(nrow=0,ncol=4)
        bindmatset=matrix(nrow=0,ncol=4)
        for(j in 1:length(starts))
        {
            motcounts=ourcounts[1:64,ourcounts[65,]==j]
            motcounts=t(motcounts)
            summary=oursummary[1:4,oursummary[5,]==j]
            expcounts=matrix(nrow=nrow(motcounts),ncol=4)
            ##add one for robustness
            for(i in 1:4) expcounts[,i]=motcounts %*% newpredfrac[j,,i]+1
            expcounts=expcounts/rowSums(expcounts)
            summary=summary[,2:(ncol(summary)-1)]
            expsummary=t(expcounts)
            summary2=summary/expsummary*rowSums(expcounts)
            noninclogprob=colSums(summary*t(log(expcounts)))
            ## have a uniform prior for bases included, four bases
            ## need likelihood for an included base
            ## uniform dirichlet prior leads to following posterior after integrating out frequencies
            inclogprob=log(6)+logfact[summary[1,]+1]+logfact[summary[2,]+1]+logfact[summary[3,]+1]+logfact[summary[4,]+1]-logfact[colSums(summary)+4]
            increl=log(incprob*exp(inclogprob-noninclogprob)+(1-incprob))
            increl[is.infinite(increl)]=(inclogprob-noninclogprob)[is.infinite(increl)]
            newrel=c(0,cumsum((increl)))
            ## get lhood for each possible start and end position
            lhood=matrix(0,nrow=ncol(summary),ncol=ncol(summary))
            for(start in 1:ncol(summary))
            {
                for(end in start:ncol(summary))
                    lhood[start,end]=log(plen)*(end-start)+log(1-plen)+newrel[end+1]-newrel[start]
            }
####these are relative log-probs
            lhood2=exp(lhood-max(lhood[lhood!=0]))
            lhood2[lhood==0]=0
            lhood2=lhood2/sum(lhood2[lhood!=0])
            motpos=sample(length(lhood2),1,prob=as.double(lhood2))
            start=motpos %% nrow(lhood2)
            if(start==0) start=nrow(lhood2)
            end=(motpos-start)/nrow(lhood2)+1
            ## have sampled new motif start and end positions
            ## sample new parameters - should be dirichlet but use expectations instead?
            temp=summary[,start:end]+1
            temp=matrix(temp,nrow=4)
            temp=t(t(temp)/colSums(temp))
            temp=matrix(log(temp),nrow=4)
#####for binding, get rid of background
            temp2=summary2[,start:end]+1
            temp2=matrix(temp2,nrow=4)
            ##temp2=t(t(temp2)/colSums(temp2))
            temp2=matrix(log(temp2),nrow=4)
            newmat=matrix(nrow=end-start+1,ncol=4)
            newmat[,1:4]=t(temp)
            newmat2=t(temp2)
###make a new combined matrix for looking at...
            newmatset=rbind(newmatset,newmat)
            bindmatset=rbind(bindmatset,newmat2)
            newnewdim=c(newnewdim,nrow(newmat))
        }



        ##
        ## set new params
        ##
        if(verbose>=2)   print(paste("Set new params - ",date(),sep=""))
        scorematsetold=scorematset
        if(updatealpha==1) alpha=alphanew
        if(updatemot==1)
        {
            scorematset=newmatset[,1:4]
            dimvec=newnewdim
        }  else {
            scorematset=scorematset[,1:4]
        }
        alphas=rbind(alphas,alpha)
        ##remove motifs if not viable
        if(min(length(fullseqs)*alpha)<=10)
        {
            print("Some motifs have <=10 expected copies, removing")
            print(which(length(fullseqs)*alpha<=10))
            newmat=matrix(nrow=0,ncol=4)
            newmat2=matrix(nrow=0,ncol=4)
            newstarts=c(1,cumsum(dimvec)+1)
            newends=cumsum(dimvec)
            for(i in 1:length(dimvec))
            {
                if(length(fullseqs)*alpha[i]>10)
                {
                    newmat=rbind(newmat,scorematset[newstarts[i]:newends[i],])
                    newmat2=rbind(newmat2,bindmatset[newstarts[i]:newends[i],])
                }
            }
            scorematset=newmat
            bindmatset=newmat2
            dimvec=dimvec[length(fullseqs)*alpha>10]
            ##remove offending motif
            alphas=alphas[,length(fullseqs)*alpha>10]
            alpha=alpha[length(fullseqs)*alpha>10]
        }
###remove motifs if not long enough
        if(min(dimvec)<=3)
        {
            print("Some motifs have  length <=3, removing")
            print(which(dimvec<=3))
            newmat=matrix(nrow=0,ncol=4)
            newmat2=matrix(nrow=0,ncol=4)
            newstarts=c(1,cumsum(dimvec)+1)
            newends=cumsum(dimvec)
            for(i in 1:length(dimvec))
            {
                if(dimvec[i]>3)
                {
                    newmat=rbind(newmat,scorematset[newstarts[i]:newends[i],])
                    newmat2=rbind(newmat2,bindmatset[newstarts[i]:newends[i],])
                }
            }
            scorematset=newmat
            bindmatset=newmat2
            scorematset=newmat
###remove offending motif
            alphas=alphas[,dimvec>3]
            alpha=alpha[dimvec>3]
            dimvec=dimvec[dimvec>3]
        }


        ##
        ## plot here
        ##
        if(verbose>=2)   print(paste("Plotting - ",date(),sep=""))  
        if(plotting=="pdf" | plotting=="screen")  plotIteration2(plotting=plotting,fileName=paste(plotPrefix,".it.",its,".pdf",sep=""),dimvec=dimvec,scorematset=scorematset,height=4,width=16,saveForHistL=saveForHistL,whichmot=whichmot,nSequences=nSequences,regprobs=regprobs,its=its)


        ##
        ## save things 
        ##
        results[[its]]=list(scoremat=scoremat)
        
    } # end of iterations



    if(verbose>=1) print(paste(" --------------------------------------------- ",sep=""))
    if(verbose>=1)  print(paste(" ----- Program Done ",date()," -----", sep=""))
    if(verbose>=1)  print(paste(" --------------------------------------------- ",sep=""))

    ##
    ##  end of iterations
    ##



    return(list(seqs=origseqs,alphas=alphas,beststrand=beststrand, trimmedseqs=fullseqs,prior=prior,alpha=alpha,bindmat=bindmatset,scoremat=scorematset,scorematdim=dimvec,regprob=regprob,regprobs=regprobs,bestmatch=bestpos,whichregs=whichregs,whichpos=whichpos,background=qvec,whichmot=whichmot, whichstrand=strand,pwm=pwm,motifEnrichment=motifEnrichment,initialData=initialData))

} # end of getMotifs
