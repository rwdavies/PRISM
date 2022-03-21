#' @param fracWithMotif what fraction contain motifs
#' @param nMotifs how many motifs
#' @param seq most likely characters in sequence
simulate_hotspots <- function(
    fracWithMotif = 0.4,
    nMotifs = 10000,
    seq = 1 + c(0,3,0,1,0,2,1,1,3,3)
) {
    ## simulate sequences with motif
    n <- length(seq)
    prob <- 0.5 + 0.5*runif(n) # probability base not chosen at random
    seqsC <- lapply(1:nMotifs,function(iMotif) {
        s <-sample(1:4,1000 + runif(1)*2000,replace=TRUE)
        if(runif(1)<fracWithMotif) {
            ## give 80% probability to being in central 200 + 200 bp
            e=array(0.2/(length(s)-n-401),length(s)-n) # use this to make prob 80% for central 200 bp
            e[(round(length(e)/2)-200):(round(length(e)/2)+200)]=0.8/401
            w=sample(length(s)-n,size=1,prob=e)+1:n
            s2=seq
            w2=runif(10)>prob
            s2[w2]=sample(1:4,sum(w2),replace=TRUE)
            s[w]=s2
        }
        ## turn into letters
        s <- paste(c("A","C","G","T")[s],collapse="")
        return(s)
    })
    return(seqsC)
}
    

    
    
