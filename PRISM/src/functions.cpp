// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <time.h>

using namespace Rcpp;





//
//
// c++ functions for PIMP
//
//


// EXAMPLE of how to export if desired
////' Reformat the reads for a single individual
////' 
////' @param x Describe parameters here.
////' @export
//// [[Rcpp::export]]



//' Add vector through a matrix
//' 
//' @param A - matrix ( explain what this is)
//' @param b - vector ( explain what this is)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix addMatVec(Rcpp::NumericMatrix A, Rcpp::NumericVector b) {
  //
  int c = b.size(); 
  // new variables
  int nrows = A.nrow();
  int ncols = A.ncol();
  int maxCols=ncols - c + 1; // 1 based - maximum possible columns
  Rcpp::NumericMatrix D(nrows,maxCols);
  // iterate over
  int i, j, k;
  for(i=0; i<nrows; i++)
  {
    for(j=0; j<maxCols; j++)
    {
      for(k=0; k<c; k++)
      {
        D(i,j)=D(i,j)+b(A(i,j+k));
      }
    }
  }
  return(D);
}





//' Add matrix through a matrix
//' 
//' @param A - matrix ( explain what this is)
//' @param B - vector ( explain what this is)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix addMatMat(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B)
{
  //Rcpp::NumericMatrix Ax = as<Rcpp::NumericMatrix>(A);
  //Rcpp::NumericMatrix Bx = as<Rcpp::NumericMatrix>(B);
  int c = B.nrow(); 
  // new variables
  int nrows = A.nrow();
  int ncols = A.ncol();
  int maxCols=ncols - c + 1; // 1 based - maximum possible columns
  Rcpp::NumericMatrix D(nrows,maxCols); // what is returned
  // iterate over
  int i, j, k;
  for(i=0; i<nrows; i++)
  {
    for(j=0; j<maxCols; j++)
    {
      for(k=0; k<c; k++)
      {
        D(i,j)=D(i,j)+B(k,A(i,j+k));
      }
    }
  }
  return(D);
}







//' For a vector y of length yT, return a (1-based) vector where entry i counts the number of times integer i-1 was seen in the data
//' 
//' @param y vector of integers
//' @param yT length of y
//' @param xT maxmimum integer to count to
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector increment(Rcpp::IntegerVector y, int yT, int xT)
{
  // new variables
  Rcpp::NumericVector x(xT+1);
  int t;
  for(t=0; t<=yT-1; t++)
    x[y[t]]++;
  return(x);
}




//' For a vector y of length yT, return a (1-based) vector where entry i counts the number of times integer i-1 was seen in the data
//' 
//' @param y vector of integers
//' @param yT length of y
//' @param xT maxmimum integer to count to
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector increment2N(Rcpp::NumericVector y, Rcpp::NumericVector z, int yT, int xT)
{
  //  int yTa = as<int>(yT);
  //int xTa = as<int>(xT);
  // Rcpp::NumericVector ya = as<Rcpp::NumericVector>(y);
  //Rcpp::NumericVector za = as<Rcpp::NumericVector>(z);
  // new variables
  //Rcpp::IntegerVector xa(xTa+1);
  Rcpp::NumericVector x(xT+1);
  int t;
  for(t=0; t<=yT-1; t++)
    x[z[t]]=x[z[t]]+y[t];
  return(x);
}
