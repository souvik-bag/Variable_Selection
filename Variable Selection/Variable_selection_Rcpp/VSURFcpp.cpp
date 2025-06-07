#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List VSURFcpp(NumericMatrix x, NumericVector y, int ntree, bool verbose, bool parallel) {
  NumericVector imp;
  List vsurf;
  
  // Run VSURF algorithm
  vsurf =   R::VSURF(x, y, ntree, verbose, parallel);
  imp = R::as.vector(vsurf$varselect.pred)
  // Return VSURF results as R list
  return imp;
}
