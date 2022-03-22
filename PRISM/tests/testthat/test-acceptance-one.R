if ( 1 == 0 ) {
    
    library("testthat")
    library("PRISM")
    dir <- "~/proj/PRISM/"
    setwd(paste0(dir, "/PRISM/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)


}

test_that("can run getMotifs using sequences", {

    ## this acceptance test runs the motif finder using 10,000 hotspots
    ## and 10,000 coldspots, where the hotspots are enriched for
    ## the motif "seq" below, and further enriched in the centre
    ## of the hotspots for that sequence
    ## the code below runs the motif finder,
    ## then check that the returned motif is a good match to the true one
    
    seq <- 1 + c(0,3,0,1,0,2,1,1,3,3) ## ATACAGCCTT or AAGGCTGTAT
    
    ## hotspots
    seqsH <- simulate_hotspots(
        fracWithMotif = 0.8,
        nMotifs = 10000,
        seq = seq
    )
    ## coldspots
    seqsC <- simulate_hotspots(
        fracWithMotif = 0,
        nMotifs = 10000,
        seq = seq
    )

    dir <- tempdir()
    dir.create(dir)
    setwd(dir)
    
    ## only do simple k-mer version
    out <- getMotifs(
        seqsH = seqsH,
        seqsC = seqsC,
        initializationMethod = "advanced",
        dataType = "sequence",
        K = 8,
        eligible = NULL,
        plotting = "pdf",
        plotPrefix = "plot",
        verbose = 2,
        nCores = 1,
        maxits = 2,
        pwmValue = 0.9,
        hotspotCenterDist = 200
    )

    ## check results here
    inferred_seq <- apply(t(out$scoremat), 2, function(x) which.max(x))

    ## checks forward vs backwards
    check_inferred_vs_true_seq(inferred_seq, seq)

})
