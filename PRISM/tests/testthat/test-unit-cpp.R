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

test_that("pass in a vector and a matrix and a length", {

    A <- matrix(as.numeric(sample(4,40,replace=TRUE)),ncol=8)-1 ## note 0-based
    b <- runif(4)

    result <- addMatVec(A=A,b=b)
    expect_equal(result[1, 1],  sum(b[1+A[1,1:4]]))
    expect_equal(result[2, 1],  sum(b[1+A[2,1:4]]))

})
                                        


test_that("pass in a vector and a matrix and a length matrix", {
    
    A <- matrix(as.numeric(sample(5,40,replace=TRUE)),ncol=8)-1 ## note 0-based
    B <- matrix(runif(15),ncol=5)

    ## note A is 0-based !
  
    result <- addMatMat(A=A,B=B)
    expect_equal(result[1, 1], sum(B[cbind(1:nrow(B),A[1,1:nrow(B)]+1)])) # entry 1,1
    expect_equal(result[2, 2], sum(B[cbind(1:nrow(B),A[2,1+1:nrow(B)]+1)])) # entry 2,2

})
  
