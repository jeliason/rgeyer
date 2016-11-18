#include <Rcpp.h>
#include <vector>

#include "Image.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rtest_im_c(List im, NumericMatrix x) {
  NumericVector out(x.nrow());
  Image z(im);


  for(int i=0; i<x.nrow(); i++) {
    out(i) =z.getValue(x(i,0), x(i,1));
  }

  return out;
}
