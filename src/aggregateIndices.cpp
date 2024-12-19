#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix updColMaxBy(NumericMatrix mu, NumericMatrix fmax, IntegerVector fac, LogicalVector skip) {
  int n = mu.nrow();
  int m = mu.ncol();
  if (Rcpp::is_na(any(skip))) {stop("NA in skip.");}
  if (m != fmax.ncol()) {stop("Wrong dimension for fmax");}
  if (fac.length() != n || skip.length() != n) {stop("Wrong length for fac or skip.");}
  for (int i = 0; i < n; i++) {
    if (!skip(i)) {
      if (fac(i) < 1 || fac(i) > fmax.nrow()) {stop("updColMax: Index out of range for fac.");}
      for (int j = 0; j < m; j++) {
       if (fmax(fac(i) - 1, j) < mu(i, j)) {
          fmax(fac(i) - 1, j) = mu(i, j);
        }
      }
    }
  }
  for (int f = 0; f < fmax.nrow(); f++) {
    for (int j = 0; j < m; j++) {
      if (fmax(f, j) <= R_NegInf) {
        fmax(f, j) = 0.0; // If max is -Inf, replace with 0.
      }
    }
  }
  return fmax;
}



// [[Rcpp::export]]
NumericMatrix updBaseMat(NumericMatrix bmat, NumericMatrix mu, NumericMatrix fmax, IntegerVector fac, NumericVector weights, LogicalVector skip) {
  int n = mu.nrow();
  int m = mu.ncol();
  int nlev = fmax.nrow();
  NumericMatrix aggr (nlev, m);
  if (Rcpp::is_na(any(skip))) {stop("NA in skip.");}
  if (m != fmax.ncol()) {stop("Wrong dimension for fmax");}
  if (bmat.ncol() != fmax.ncol() || bmat.nrow() != fmax.nrow()) {stop("dimension mismatch bmat-fmax");}
  if (fac.length() != n || skip.length() != n || weights.length() != n) {stop("Wrong length for fac or skip or weights.");}
  std::fill(aggr.begin(), aggr.end(), R_NaReal);
  // Skip first loop if no aggregation needed?
  for (int i = 0; i < n; i++) {
    if (!skip(i)) {
      if (fac(i) < 1 || fac(i) > fmax.nrow()) {stop("updBaseMat: Index out of range for fac.");}
      for (int j = 0; j < m; j++) {
        if (R_IsNA(aggr(fac(i) - 1, j))) {
          aggr(fac(i) - 1, j) = weights(i) * exp(mu(i, j) - fmax(fac(i) - 1, j));
        } else {
          aggr(fac(i) - 1, j) += weights(i) * exp(mu(i, j) - fmax(fac(i) - 1, j));
        }
      }
    }
  }

  for (int f = 0; f < nlev; f++) {
    for (int j = 0; j < m; j++) {
      if (!R_IsNA(aggr(f, j))) {
        if (R_IsNA(bmat(f, j))) {
          bmat(f, j) = fmax(f, j) + log(aggr(f, j));
        } else {
          bmat(f, j) = fmax(f, j) + log(exp(bmat(f, j) - fmax(f, j)) + aggr(f, j));
        }
      }
    }
  }
 return bmat;
}

