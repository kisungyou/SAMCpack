// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"

using namespace Rcpp;
using namespace arma;

/////////////////////////////////////////////////////////////////////////////////////////////////
// Export Functions
// [[Rcpp::export]]
Rcpp::List exec_SAMCoptim(Function func, const int nv, arma::vec& energy, arma::mat& domain,
                     const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize,
                     arma::mat& trange,  arma::vec& init){
  
  Rcpp::List hey;
  return(hey);
}