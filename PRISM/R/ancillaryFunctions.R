print_message <- function(x, include_mem = FALSE) {
    if (include_mem) {
        mem <- system("ps auxww | grep 'scripts/profile.R' | grep slave | grep -v 'grep' | awk -v OFS='\t' '$1=$1' | cut -f6", intern = TRUE)
        if (length(mem) > 0) {
            mem <- paste0(paste0(round(as.integer(mem) / 2 ** 20, 3), collapse = ", "), " - ")
        } else {
            mem <- ""
        }
    } else {
        mem <- ""
    }
    message(
        paste0(
            "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", mem, x
        )
    )
}


check_mclapply_OK <- function(out, stop_message = "An error occured during QUILT. The first such error is above") {
    te <- (sapply(out, class) == "try-error") | sapply(out, is.null)
    if (sum(te) > 0) {
        print_message(out[[which(te)[1]]]) # print first error
        stop(stop_message)
    }
    return(NULL)
}




## some hashing functions
hashF=function(x,K) sum(x*4^(0:(K-1)))
hashFC=function(x,K) sum((3-x[K:1])*4^(0:(K-1)))
hashM=function(x,K)  min(hashF(x,K),hashFC(x,K))
conv=function(x,K) 3-x[K:1]

simonHash2 <- function(seqv,K, check = FALSE) {
    ## fixed from simonHash - was backwards!!!
    nonrep <- array(0,length(seqv)-K+1)
    for(i in 1:K) {
        z=seqv[i:(i+length(nonrep)-1)]
        nonrep=nonrep+4^(i-1)*z # flipped K-i to i-K
        nonrep[z>3]=-Inf
    }
    nonrep2 <- nonrep[!is.infinite(nonrep)]+1
    if (check) {
        if (sum(is.na(nonrep2)) > 0) {
            stop("Something went wrong, simonHash2 has NA values")
        } else if (min(nonrep2) < 1) {
            stop("Something went wrong, simonHash2 wants to return values < 1")
        } else if (max(nonrep2) > 4 ** K) {
            stop("Something went wrong, simonHash2 wants to return values > 4 ** K")
        }
    }
    return(nonrep2)
}
  



# make compO matrix
makeCompO2=function(r,which,hashLink,K,num)
{
  # no point using 0 entries at all!
  toKeep=which
  #toKeep=rowSums(r>0)>0
  if(sum(which)>1)
    compO=cbind(r[toKeep,],hashLink[toKeep,])
  if(sum(which)==1)
    compO=matrix(c(r[toKeep,],hashLink[toKeep,]),nrow=1)
  if(sum(which)==0)
    return(compO=NULL)
  colnames(compO)=c(rep("",ncol(r)),colnames(hashLink))
  # if only exists in F - OK
  # if it exists in R - add to current value, remove R
  a=match(compO[,"F"],compO[,"R"])
  a[compO[,"F"]==compO[,"R"]]=NA
  # in the rare case A is only NA - no need to do next steps
  if(sum(is.na(a)==FALSE)>0)
  {
    b=cbind((1:length(a))[is.na(a)==FALSE],a[is.na(a)==FALSE])
    b[b[,2]<b[,1],c(1,2)]=b[b[,2]<b[,1],c(2,1)]
    which=match(unique(b[,1]),b[,1])
    compO[b[which,1],1:num]=compO[b[which,1],1:num]  + compO[b[which,2],1:num]
    # now remove the second entry
    compO=compO[-b[which,2],]
  }
  # at this point apply filters as well
  #complex=compO[,"mr"]<=mrle & compO[,"nd"]>=ndge
  #compO=compO[complex==TRUE,]
  if(is.null(dim(compO)))
    compO=NULL # no point moving forward with one dimension
  return(compO)
}



# generate p-values for some counts! 
getORAndP=function(compO=NULL,K=NULL,num,nCores=1,verbose=FALSE,data=NULL,gc.content=NULL,CG=NULL,version="robbie",comp2=NULL,outFile=NULL,cgte=0,pvalmet="fishers",rgte=0)
{
  #
  # modifications - jan 2015 - do filtering after p-value calculation
  #  
  # further modifications - better use of type (numeric) - modification september 9, 2014
  # 
  if(version=="afi")
  {
    # afi version
    type=as.integer(gc.content*K + CG)
    typeList=names(table(type))
    # get first data
    x1=data
  } else if(version=="robbie") {
    #
    ### define type
    #
    # stratify by AT/GC content and CpGs
    type=as.integer(compO[,"at"]*K + compO[,"cpgs"])
   #type=as.integer(compO[,"mr"]*K+compO[,"at"]*6*K + compO[,"cpgs"])
    ### new code for more refined granular comparison
    x=compO[,(num+4+1):(num+4+K)]
    if(1==0) {
      nucleo=cbind(rowSums(x==0),rowSums(x==1),rowSums(x==2),rowSums(x==3))
      nucleoM=nucleo
      nucleoM[nucleoM[,4]<nucleoM[,1],c(1,4)]=nucleoM[nucleoM[,4]<nucleoM[,1],c(4,1)]
      nucleoM[nucleoM[,3]<nucleoM[,2],c(2,3)]=nucleoM[nucleoM[,3]<nucleoM[,2],c(3,2)]
      type=as.integer(nucleoM[,1] + 10*nucleoM[,2] + 100*nucleoM[,3] + 1000*nucleoM[,4] + 10000*compO[,"cpgs"]+5*10000*compO[,"mr"])
    }
    ###
    typeN=as.numeric(type)
    typeList=names(table(type))
    # we need 4 things or make the OR
    # first, we need each entry for every SNP under consideration
    x1=compO[,1:num]
  } else {
    return(NULL) # only two version supported
  }
  #
  ### NEW - SEPT 9 2014 - allow second way to calculate p-values
  #
  if(version=="robbie" & is.null(comp2))
  {
    ### second, for every SNP type, want the sum against every other group
    x2=rowSums(x1)-x1
    x3b=apply(x1,2,function(xxx) increment2N(y=as.numeric(xxx),z=typeN,yT=as.integer(dim(x1)[1]),xT=as.integer(range(type)[2])))
    x3=apply(x3b,2,function(y) y[type+1]) - x1
    x4b=increment2N(y=as.numeric(rowSums(x1)),z=typeN,yT=as.integer(dim(x1)[1]),xT=as.integer(range(type)[2]))
    x4=x4b[type+1] - x3 - x2 - x1
  } else {
    x2=comp2[,1:num]
    x3b=apply(x1,2,function(xxx) increment2N(y=as.numeric(xxx),z=typeN,yT=as.integer(dim(x1)[1]),xT=as.integer(range(type)[2])))
    x3=apply(x3b,2,function(y) y[type+1]) - x1
    # x4 is now in a very similar way
    x4b=apply(x2,2,function(xxx) increment2N(y=as.numeric(xxx),z=typeN,yT=as.integer(dim(x1)[1]),xT=as.integer(range(type)[2])))
    # OVERFLOW! ARGH
    x4=apply(x4b,2,function(y) y[type+1]) - x2
     #i=10;j=8 # uncomment to try manual check
     #c(x1[i,j],x2[i,j],x3[i,j],x4[i,j])
     #c(x1[i,j],x2[i,j],sum(x1[type[i]==type,j])-x1[i,j],sum(x2[type[i]==type,j])-x2[i,j])
  }
  # for those with only a single entry - this won't work. blank later
  a=table(typeN)
  toBlank=is.na(match(typeN,names(a)[a==1]))==FALSE
  x1[toBlank,]=0 # they will get filtered out later
  x2[toBlank,]=5 # they will get filtered out later
  x3[toBlank,]=10 # they will get filtered out later
  x4[toBlank,]=20 # they will get filtered out later
  #
  ### now, calculate p-values!
  # make some ORs
  or=x1/x3*x4/x2 # 
  se=sqrt(1/x1+1/x2+1/x3+1/x4) # OR for log OR
  #
  # calculate p-values here
  #
  p=mclapply(1:ncol(x1),pvalmet=pvalmet,function(i,pvalmet)
  {
   #which=1:nrow(x1)
   #which=typeN==590
    a=x1[,i]
    b=x2[,i]
    c=x3[,i]
    d=x4[,i]
    n=a+b+c+d
    p=array(NA,length(a)) # set to NA unless at least 5 entries
    which=a>=cgte
    # entries are m, n, k, lo, hi, x from interior of Fisher's Exact test code see tempFISHER.R
    m=cbind(a+c,b+d,a+b,0,a+b,a)
    t1=m[,3]-m[,2];    m[t1>0,4]=t1[t1>0]
    t1=m[,1]>m[,3];    m[t1,5]= m[t1,3]
    m2=cbind(a,b,c,d)
    p=apply(m2,1,fastChisq)
    if(verbose==TRUE) print(paste(i,",sum(which),",sum(which),",",date(),sep=""))
    if(pvalmet=="chisq" & sum(which)>1) ps=apply(m2[which,],1,fastChisq)
    if(pvalmet=="chisq" & sum(which)==1) ps=fastChisq(m2[which,])
    if(pvalmet=="chisq" & sum(which)==0) ps=1
    if(pvalmet=="fishers" & sum(which)>1) ps=apply(m[which,],1,fastFishers)
    if(pvalmet=="fishers" & sum(which)==1) ps=fastFishers(m[which,])
    if(pvalmet=="fishers" & sum(which)==0) ps=1
    if(verbose==TRUE) print(paste(i,",done,",date(),sep=""))
    p[which]=ps
    # now do fisher's exact if sum is too small
    p[p=="NaN"]=NA
    return(p)
  },mc.cores=nCores)
  # fill in properly
  p=sapply(p,function(x) x)
  or[is.na(p)==TRUE]=NA
  se[is.na(p)==TRUE]=NA
  p[toBlank,]=NA # those from single entries
  or[toBlank,]=NA
  se[toBlank,]=NA
  #
  ### now - make things smaller
  #
  # now - remove if the row is too small
  which=rowSums(compO[,1:num])>=rgte
  compO=makeSubMatrix(compO,which=which,num+4+K) # all we need for p-value 1
  if(is.null(comp2)==FALSE)
    comp2=makeSubMatrix(comp2,which=which,num+4+K) # all we need for p-value 1  
  p=makeSubMatrix(p,which=which,num)
  or=makeSubMatrix(or,which=which,num)
  se=makeSubMatrix(se,which=which,num)
  # make sure each one is a matrix
  if(is.null(outFile))
  {
    return(list(p=p,or=or,se=se))
  } else {
    # only save if dimension greater than 0
    if(nrow(compO)>0)
      save(compO,comp2,p,or,se,file=outFile) 
  }
}



makeSubMatrix=function(x,which,w)
      {
        y=x[which,]
        if((length(y)/(w))==1) {
          z=matrix(y,nrow=1)
          colnames(z)=names(y)
          y=z
        }
        return(y)
}
      


fastChisq=function(x) {
      x=matrix(x,ncol=2)
        nr <- as.integer(nrow(x))
        nc <- as.integer(ncol(x))
        sr <- rowSums(x)
        sc <- colSums(x)
      n=sum(x)
        E <- outer(sr, sc, "*")/n
        YATES <- 0
        STATISTIC <- sum((abs(x - E) - YATES)^2/E)
        PARAMETER <- (nr - 1L) * (nc - 1L)
        PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
      return(PVAL)
    }


fastFishers=function(y)
{
        logdc <- dhyper(y[4]:y[5], y[1], y[2], y[3], log = TRUE) # range
        d <- logdc + log(1) * y[4]:y[5]
        d <- exp(d - max(d))
        d = d/sum(d)
        return(sum(d[d <= d[y[6] - y[4] + 1] * (1 + 10^(-7))]))
}

buildHashLink=function(K)
{
    # build up has link here
    hashLink=array(0,c(4^K,K+1))
    hashLink[,1]=0:(4^K-1)
    for(k in 1:(K-1))  hashLink[,k+1]=  rep(sort(rep(0:3,4^(k-1))),4^(K-k-1))
    hashLink[,K+1]=  sort(rep(0:3,4^(K-1)))
    # add in hash of converse
    a=apply(hashLink[,-1],1,function(x) hashFC(x,K))
    hashLink=cbind(hashLink[,1],a,hashLink[,-1])
    colnames(hashLink)=c("F","R",rep("",dim(hashLink)[2]-2))
    ### for each entry in the hash, calculate number of CpGs and GC content
    x=  hashLink[,-c(1:2)]
    at=rowSums(x==0 | x==3) # AT number
    cpgs=rowSums(x[,-1]==2 & x[,-K]==1 )
    ### also flag low complexity ones
    # maximum run length - for each succesive run, ask if it exists, if yes, then store
    hashLink=cbind(at,cpgs,hashLink)
    colnames(hashLink)[1:2]=c("at","cpgs")
    ### save this variable for now
    # add rownames as human readable
    rownames(hashLink)=apply(hashLink[,-c(1:4)],1,function(x) paste(c("A","C","G","T")[x+1],collapse="",sep=""))
    return(hashLink)
} # end of hashLink building function



# for a motif seed, find close words
findWords=function(motifSeed,eligible,plineage,orlineage,K,compOX,pThresh,pwmMinimum)
{
  #
  ### iterate
  #
  # first iteration
  matches=as.list(1:10)
  for(i in 1:10) matches[[i]]=list(0,0,0,0)
  cct=pThresh # cct = current cluster threshold
  matches[[1]]=getCloseMatches(motifSeed,eligible,plineage=plineage,K,compOX,first=TRUE,iteration=1,cct=cct)
  eligible=matches[[1]]$eligible
  # if there are NO entries here - we're done
  if(length(matches[[1]][[1]])>K)
  {
    # start with a seed
    iteration=1
    while(iteration<100)
    {
      iteration=iteration+1
      # go through each entry, get matches
      matchesP=matches[[iteration-1]][[1]]
      whereP=matches[[iteration-1]][[2]]
      # get new cluster threshold
      # basically, take old one, and add in new expected tests
      cct=0.05 / (0.05/cct + length(matches[[iteration-1]][[2]])*328)
      #cct=0.05 / (0.05/cct + sum(unlist(lapply(matches[[iteration-1]],function(x) return(x$numberTested)))))
      ### now loop  
      n=length(matchesP)/K
      matchesL=as.list(1:n)
      for(i in 1:n)
      {
        matchesL[[i]]=getCloseMatches(motifSeed=matrix(matchesP,ncol=K)[i,],eligible,plineage=plineage,K,compOX,iteration=iteration,cct=cct)
        eligible=matchesL[[i]]$eligible # update which ones are eligible
        # adjust distance based on where it came from
        matchesL[[i]][[2]]=matchesL[[i]][[2]]+whereP[i]
      }
      ### reformat 
      matchesN=matrix(unlist(lapply(matchesL,function(x) t(x[[1]]))),ncol=K,byrow=TRUE)
      matchesN=matrix(matchesN,ncol=K)
      whereN=unlist(lapply(matchesL,function(x) return(x[[2]])))
      namesN=unlist(lapply(matchesL,function(x) return(x[[3]])))
      whereItMatchesToN=unlist(lapply(matchesL,function(x) return(x[[4]])))
      fromWhereN=unlist(lapply(matchesL,function(x) return(x[[5]])))      
      ### remove those we've seen before
      remove=match(whereItMatchesToN,unlist(lapply(matches[1:(iteration-1)],function(x) x[[4]])))
      matchesN=matrix(matchesN[is.na(remove)==TRUE,],ncol=K)
      whereN=whereN[is.na(remove)==TRUE]
      namesN=namesN[is.na(remove)==TRUE]
      whereItMatchesToN=whereItMatchesToN[is.na(remove)==TRUE]
      fromWhereN=fromWhereN[is.na(remove)==TRUE]      
      ### also remove duplicates
      keep=match(unique(whereItMatchesToN),whereItMatchesToN)
      matches[[iteration]]=list(matchesN[keep,],whereN[keep],namesN[keep],whereItMatchesToN[keep],fromWhereN[keep])
      if(length(namesN)==0)
      {
        matches=matches[1:(iteration-1)]
        iteration=100
      }
    } # end of iteration loop
    ###
    ### reformat output
    ###
    matchesN=matrix(unlist(lapply(matches,function(x) t(x[[1]]))),ncol=K,byrow=TRUE)
    whereN=unlist(lapply(matches,function(x) return(x[[2]])))
    namesN=unlist(lapply(matches,function(x) return(x[[3]])))
    whereItMatchesToN=unlist(lapply(matches,function(x) return(x[[4]]))) # reorder
    fromWhereN=unlist(lapply(matches,function(x) return(x[[5]]))) # reorder    
    a=order(whereN)
    matchesN=matchesN[a,]
    whereN=whereN[a]
    namesN=namesN[a]
    fromWhereN=fromWhereN[a]
    whereItMatchesToN=whereItMatchesToN[a]        
  } else {
    # this is only if there are no matches
    matchesN=matrix(motifSeed,ncol=K)
    whereN=0
    x=hashF(motifSeed,K)
      t1=match(x,compOX[,"F"])
      t2=match(x,compOX[,"R"])
    t1[is.na(t1)]=t2[is.na(t1)] # t1 is where is matches to
    whereItMatchesToN=t1
    fromWhereN=0
    namesN=whereItMatchesToN
  }
  ###
  ### make a version we can display
  ###
  border=max(abs(range(whereN)))
  n=length(namesN)
  textOut=array(0,n)
  for(i in 1:n)
  {
    textOut[i]=paste(paste(rep(" ",border-whereN[i]),collapse=""),paste(c("A","C","G","T")[1+matchesN[i,]],collapse=""),paste(rep(" ",border+whereN[i]),collapse=""),sep="")
  }
  ###
  ### also count up every position
  ###
  x=lapply(textOut,function(x) strsplit(x,""))
  sum=array(0,c(5,K+2*border))
  orL=orlineage[whereItMatchesToN]   # or of motifs
  for(i in 1:length(textOut))
  {
    t1= match(unlist(x[[i]]),c("A","C","G","T"))
    for(j in 1:4) {    a=is.na(t1)==FALSE & t1==j;    sum[j,a]=sum[j,a]+orL[i]-1}
    # use OR-1 to determine summing
    sum[5,is.na(t1)==FALSE]=  sum[5,is.na(t1)==FALSE]+orL[i]-1
  }
  rownames(sum)=c("A","C","G","T","Sum")
  sum=sum[-5,]
  sumNum=sum
  sum[sum==0]=pwmMinimum
  sum=apply(sum,2,function(x) x/sum(x))
  ### return
  return(list(matchesN=matchesN,whereN=whereN,namesN=namesN,whereItMatchesToN=whereItMatchesToN,textOut=textOut,sum=sum,sumNum=sumNum,fromWhereN=fromWhereN,eligible=eligible))
}




### modified - send in p-values and if eligible, and a threshold
getCloseMatches=function(motifSeed,eligible,plineage,K,compOX,first=FALSE,iteration,cct)
{
  where=array(0,(K-1)*4*4*2+K*4)
  motifArray=array(0,c((K-1)*4*4*2+K*4,K))
  # get hashes for everything within 1, no moves
  for(i in 1:K)
  {    
    for(j in 0:3)
    {     
      motifArray[1+j+4*(i-1),]=motifSeed
      motifArray[1+j+4*(i-1),i]=j
      where[1+j+4*(i-1)]=0
    }
  }
  # now off by 1 + offset
  for(i in 1:(K-1))
  {    
    for(j1 in 0:3)
    {
      for(j2 in 0:3)
      {
        # new start base
        motifArray[1+4*K+4*4*(i-1)+4*j2+j1,]=c(j1,motifSeed[-K])
        motifArray[1+4*K+4*4*(i-1)+4*j2+j1,i+1]=j2
        where[1+4*K+4*4*(i-1)+4*j2+j1]=1
        # MODIFIED Sept 30, 2014 - flip definition of where
      }
    }
  }
  for(i in 1:(K-1))
  {    
    for(j1 in 0:3)
    {
      for(j2 in 0:3)
      {
        # new end base
        motifArray[1+4*K+4*4*(K-1)+4*4*(i-1)+4*j2+j1,]=c(motifSeed[-1],j1)
        motifArray[1+4*K+4*4*(K-1)+4*4*(i-1)+4*j2+j1,i]=j2
        where[1+4*K+4*4*(K-1)+4*4*(i-1)+4*j2+j1]=-1
      }
    }
  }
  # now get hash functions for these
  x=apply(motifArray,1,hashF,K=K)
  rownames(motifArray)=x
  # remove duplicates
  motifArray=motifArray[match(unique(x),x),]
  where=where[match(unique(x),x)]  
  x=unique(x)
  #### modified feb 2015 -  ALSO - remove 1 from a pair of converses
  y=apply(motifArray,1,hashFC,K=K)
  z=match(x,y)
  if(sum(is.na(z))>0)
  {
    M=cbind(x[is.na(match(x,y))==FALSE],y[is.na(match(x,y))==FALSE])
    a=M[,1]
    a[M[,1]<M[,2]]=M[M[,1]<M[,2],2]
    # from a pair - keep the smaller one
    # BUT - if the motif seed is there, keep it!!!
    mh=hashF(motifSeed,K=K)
    if(is.na(match(mh,a))==TRUE)
      a[match(mh,M[,1])]=mh
    z2=unique(a) # these we want to keep
    # end of make sure motif seed is in there
    e=array(FALSE,length(x))
    e[match(z2,x)]=TRUE
    keep = is.na(match(y,x))==TRUE |  e
    # now apply - with this, we remove one instance of converses
    motifArray=motifArray[keep,]
    where=where[keep]  
    x=unique(x[keep])
  }
  ### end of modification feb 8, 2015
  fromWhere=rep(iteration,length(x))
  if(first==TRUE) fromWhere[colSums(apply(motifArray,1,function(x) x==motifSeed))==K]=0
  # remove the motif itself, except the first time
  if(first==FALSE)
  {
    e=hashF(motifSeed,K=K)
    y=match(e,x)
    motifArray=motifArray[-y,]
    where=where[-y]
    x=x[-y]
    fromWhere=fromWhere[-y]
  }
  ### now find in compOL
  t1=match(x,compOX[,"F"])
  t2=match(x,compOX[,"R"])
  t1[is.na(t1)]=t2[is.na(t1)] # t1 is where is matches to
  ### remove if neither can be matched
  remove=is.na(t1)
  motifArray=matrix(motifArray[remove==FALSE,],ncol=K)
  where=where[remove==FALSE]
  fromWhere=fromWhere[remove==FALSE]  
  x=x[remove==FALSE]
  t1=t1[remove==FALSE]
  # keep the eligibles
  which1=eligible[t1] & plineage[t1]<cct
  # remove these from being eligible in the future
  eligible[t1][plineage[t1]<cct]=FALSE 
  # return the eligibles
  return(list(a=matrix(motifArray[which1,],ncol=K),w=where[which1],x=x[which1],t1[which1],fromWhere=fromWhere[which1],eligible=eligible))
}


# can plot multiple PWMs as once!
plotIteration2=function(plotting,fileName,dimvec,scorematset,height=6,width=16,saveForHistL,whichmot,nSequences,regprobs,its)
{
  nMotifs=length(dimvec)
  height=height * (nMotifs) # upsize!
  #
  if(plotting=="pdf") pdf(fileName,height=height,width=width)
  newstarts=c(1,cumsum(dimvec)+1)
  newends=cumsum(dimvec)
  for(i in 1:length(dimvec))
  {
    iB=nMotifs+1-i
    #
    # add probability hotspot contains motif
    #
    new=TRUE
    if(i==1) new=FALSE
    par(fig=c((0)/4,1/4, (iB-1)/(nMotifs),(iB)/(nMotifs) ),new=new)
    v=hist(regprobs[,i],breaks=seq(0,1,0.05),plot=TRUE,main=paste("# sampled = ",sum(whichmot==i)," / ",nSequences,"\n# >50% = ",sum(regprobs[,i]>0.50),sep=""),xlab="Probability hotspot contains motif")
    #
    # add position within hotspot
    #
    par(fig=c((1)/4,2/4, (iB-1)/(nMotifs),(iB)/(nMotifs) ),new=TRUE)
    v=hist(saveForHistL[[i]],breaks=seq(0,1,0.05),plot=TRUE,main=paste("motif ",i,", iteration ",its,sep=""),xlab="Position within hotspot")
    #
    # motifs themselves
    #
    for(k in 1:2)
    {
      par(fig=c((k+1)/4,(k+2)/4, (iB-1)/(nMotifs),(iB)/(nMotifs) ),new=TRUE)
      # x1, x2, y1, y2
      testmat=scorematset[newstarts[i]:newends[i],]
      compmat=testmat[,c(4:1)]
      compmat=compmat[nrow(compmat):1,]
      if(k==2) testmat=compmat
      # sequence
      newSeq(sum=t(exp(testmat)),n=-1,outputPlot=FALSE,pdfName=NA,main=c("Forward","Reverse")[k])
    }
  }
  if(plotting=="pdf") dev.off()
}



##tidy up the endpoints using 0.25
setUpBackground=function(bg,qvec,originalSequenceLength,newmat,bgstore,qvecstore,bg2)
{
  if(length(bg)==1)
  {
    if(length(qvec)==0)
    {
      bg=array(0,64)
      bg2=matrix(0,nrow=length(originalSequenceLength),ncol=length(bg))
      for(i in 1:length(originalSequenceLength))
      {
        #if(!i%%1000){ print(i);print(bg)}
        temp=newmat[i,1:originalSequenceLength[i]]
        newtemp=1:(length(temp)-2)
        newtemp=16*(temp[1:(length(temp)-2)]-1)+4*(temp[2:(length(temp)-1)]-1)+(temp[3:(length(temp)-0)]-1)
        # mask out those which are 5s
        newtemp[temp[1:(length(temp)-2)]==5]=-1
        newtemp[temp[2:(length(temp)-1)]==5]=-1
        newtemp[temp[3:(length(temp)-0)]==5]=-1
        ##!!! 
        for(k in 1:64)  bg2[i,k]=bg2[i,k]+sum(newtemp==(k-1))
        bg=bg+bg2[i,]
      } 
      ####allow either strand because under null don't expect this to matter
      ####just getting triplet probs
      e=c("A","C","G","T")
      for(k in 1:2) e=c(paste("A",e,sep=""),paste("C",e,sep=""),paste("G",e,sep=""),paste("T",e,sep=""))
      d=c("T","G","C","A")
      for(k in 1:2) d=c(paste(d,"T",sep=""),paste(d,"G",sep=""),paste(d,"C",sep=""),paste(d,"A",sep=""))
      names(bg)=e
      bg=bg+bg[d]
      #####add one on for the prior, and to make robust to lack of data
      bg=bg+1
      bgstore=bg
      qvec=as.vector(bg)/sum(bg)
      qvecstore=qvec
    }
  }
  return(list(bg=bg,bgstore=bgstore,qvec=qvec,qvecstore=qvecstore,bg2=bg2))
}




# make a background matrix 
makeBackground2=function(starts,ends,qvec,newqvec,index,nSequences,maxwidth)
{
  ###now use qvec to get scores for all our sequences under the null
  ###first, get conditional probabilities for all triplets
  ###adjust below!!
  letter1=floor(((1:64)-1)/16)+1
  letter2=floor(((1:64)-(letter1-1)*16-1)/4)+1
  letter3=floor(((1:64)-(letter1-1)*16-(letter2-1)*4-1))+1
  newqvec=array(0,64)
  for(i in 1:64)  newqvec[i]=qvec[i]/sum(qvec[letter1==letter1[i] & letter2==letter2[i]])
  #####add an extra term for missing bases
  newqvec=c(newqvec,0.25)
  newqvec=log(newqvec)
  ####now  apply to the matrix storing our triplet lookup information
  ###take logs
  ###now add up to give the background equivalent to the length of the motif
  background=matrix(0,nrow=nSequences*length(starts),ncol=maxwidth)  
  background[is.na(background)]=0
  #
  # make background, add to it
  #
  background=matrix(0,nrow=nSequences*length(starts),ncol=maxwidth)  
  for(j in 1:length(starts)) # do for multiple motifs
  {
    w1=(nSequences*(j-1)+1):(nSequences*j)
    motifLength=ends[j]-starts[j]+1
    ### calculate new background
    xxx=addMatMat(index-1,matrix(rep(newqvec,motifLength),nrow=motifLength,byrow=TRUE))
    ### pad up the rest with 0s
    cols=maxwidth-motifLength+1
    background[w1,1:cols] =background[w1,1:cols] + xxx
  }
  return(list(background=background))
}


getPriorMat=function(nSequences,starts,originalSequenceLength,motifLengths,posmat,prior,alpha)
{
  rm(priormat,endpos)
  endpos=array(0,nSequences * length(starts))
  for(j in 1:length(starts))
    endpos[(nSequences * (j-1) + 1): (nSequences*j)]=originalSequenceLength-motifLengths[j]+1
  endpos=c(endpos)
  priormat=posmat/endpos
  priormat[priormat>1]=0 # get rid of the other side
  priormat=priormat-1/endpos
  priormat=floor(priormat*10)+1
  for(i in 1:10) priormat[priormat==i]=prior[i]
  priormat[priormat<=0]=0
  priormat[priormat>10]=0
  priormat=priormat/rowSums(priormat)/2
  # october 6 2015 - added into here - cleaner
  for(j in 1:length(starts))
    priormat[(nSequences*(j-1)+1):(nSequences*j),]=priormat[(nSequences*(j-1)+1):(nSequences*j),]*alpha[j];
  return(list(endpos=endpos,priormat=priormat))
}


getPredFrac=function(starts,qvec,newqvec)
{
  letter1=floor(((1:64)-1)/16)+1
  letter2=floor(((1:64)-(letter1-1)*16-1)/4)+1
  letter3=floor(((1:64)-(letter1-1)*16-(letter2-1)*4-1))+1
  predfrac=matrix(nrow=64,ncol=4)
  for(i in 1:64){
  predfrac[i,1]=qvec[letter1==letter1[i] & letter3==letter3[i] & letter2==1]
  predfrac[i,2]=qvec[letter1==letter1[i] & letter3==letter3[i] & letter2==2]
  predfrac[i,3]=qvec[letter1==letter1[i] & letter3==letter3[i] & letter2==3]
  predfrac[i,4]=qvec[letter1==letter1[i] & letter3==letter3[i] & letter2==4]
  }
  predfrac=predfrac/rowSums(predfrac)
  newpredfrac=array(0,dim=c(length(starts),64,4))
  for(i in 1:64){
  newpredfrac[,i,1]=newqvec[,letter1==letter1[i] & letter3==letter3[i] & letter2==1]
  newpredfrac[,i,2]=newqvec[,letter1==letter1[i] & letter3==letter3[i] & letter2==2]
  newpredfrac[,i,3]=newqvec[,letter1==letter1[i] & letter3==letter3[i] & letter2==3]
  newpredfrac[,i,4]=newqvec[,letter1==letter1[i] & letter3==letter3[i] & letter2==4]
  }
  for(j in 1:length(starts))
    newpredfrac[j,,]=newpredfrac[j,,]/rowSums(newpredfrac[j,,])
  return(list(predfrac=predfrac,newpredfrac=newpredfrac))
}



###
### plotting functions - from previous functions
###


### make better sequence logo plot
newSeq=function(sum,n,outputPlot=FALSE,pdfName=NULL,main=NULL,height=4,width=8)
{
  #
  ### make plottable PWM from this
  #
  # turn into ic
  #en=(1/log(2))*(4-1)/(2*n)
  Hi=-colSums(( sum * log2(sum)))
  # Ri=2-(Hi+en)
  Ri=2-(Hi)  
  # relative height - easy
  relativeHeight=t(Ri*t(sum)) # height of each base
  m=ncol(sum)
  bases=matrix(rep(c("A","C","G","T"),m),ncol=m)
  order=apply(relativeHeight,2,order)
  # now reorder and make cumulative sum
  for(i in 1:m)
  {
    bases[,i]=bases[,i][order[,i]]
    relativeHeight[,i]=relativeHeight[,i][order[,i]]
  }
  toPlotHeight=rbind(0,apply(relativeHeight,2,cumsum))
  #
  ### can now plot
  #
  if(outputPlot==TRUE) pdf(pdfName,height=height,width=width)
  # open a plot window
  plot(0,col="white",xlab="",ylab="",xlim=c(0,m+1),ylim=c(0,2),axes=FALSE,main=main)
  # add some letters
  for(iB in 1:m)
  {
    for(iL in 1:4)
    {
      # try making as rectangle
      a=letterR(x.pos=iB,y.pos=toPlotHeight[iL,iB],ht=relativeHeight[iL,iB],wt=1,letter=bases[iL,iB])
      for(j in unique(a$id))
        polygon(a$x[a$id==j],a$y[a$id==j],col=a$fill[j],border=NA)
    } # end of letter
  } # end of base
  # add axes
  mtext("IC",side=2,at=1,line=3)
  axis(side=2,at=seq(0,2,0.5),labels=prettyNum(seq(0,2,0.5)),las=2)
  if(outputPlot==TRUE) dev.off()
}

letterR=function(x.pos,y.pos,ht,wt,letter)
{
  # source functions
  if(letter=="A") return(letterRA(x.pos,y.pos,ht,wt))
  if(letter=="C") return(letterRC(x.pos,y.pos,ht,wt))
  if(letter=="G") return(letterRG(x.pos,y.pos,ht,wt))
  if(letter=="T") return(letterRT(x.pos,y.pos,ht,wt))
}


letterRA <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
  y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- c(rep(1,9),rep(2,4))
  }else{
    id <- c(rep(id,9),rep(id+1,4))
  }
    
  fill <- c("green","white")
    
  list(x=x,y=y,id=id,fill=fill)
}

## T
letterRT <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,8)
  }else{
    id <- rep(id,8)
  }
  
  fill <- "red"
    
  list(x=x,y=y,id=id,fill=fill)
}

## C
letterRC <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }
  
  fill <- "blue"
    
  list(x=x,y=y,id=id,fill=fill)
}


## G
letterRG <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  h1 <- max(y.l1)
  r1 <- max(x.l1)
  
  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)

  

  if (is.null(id)){
    id <- c(rep(1,length(x)),rep(2,length(x.add)))
  }else{
    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
  }

  x <- c(rev(x),x.add)
  y <- c(rev(y),y.add)
    
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  
  fill <- c("orange","orange")
    
  list(x=x,y=y,id=id,fill=fill)

}


