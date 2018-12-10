#ifndef _SAMCpack_CPP_AUXILIARY_H
#define _SAMCpack_CPP_AUXILIARY_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Auxiliary 1 : Sampling from random-walk proposal
arma::vec sampling_rw(arma::vec xold, arma::mat domain, const double stepsize);
arma::vec sampling_rwvec(arma::vec xold, arma::mat domain, arma::vec stepsize);
// Auxiliary 2 : find a suitable location for energy vector
int find_location(double x, arma::vec y);
// Auxiliary 3 : adjust weight update
arma::vec adjust_weights(arma::vec& weight,arma::mat& trange);
// Auxiliary 4 : find min, max
arma::vec find_minmax(arma::vec x);
arma::vec find_rowminmax(arma::rowvec x);
// Auxiliary 5 : rescale stuffs
arma::colvec rescale_colvec(arma::colvec x);
arma::rowvec rescale_rowvec(arma::rowvec x);
arma::mat rescale_vert2(arma::mat A);
arma::mat rescale_hori2(arma::mat A);
// Auxiliary 6 : get grid vector, whow!
arma::vec get_gridvec(const double pstart, const double pend, const int N);
// Auxiliary 7 : evaludate CSAMC2 g(lambda(x))
double evalcsamc_g2(arma::mat g, arma::vec gvec1, arma::vec gvec2, arma::vec lbdx);
// Auxiliary 8 : kernel density update for 2dimensional case
arma::mat evalcsamc_ks2(arma::vec vec1, arma::vec vec2, arma::mat ysamples, arma::mat H);
// Auxiliary 9 : initialize theta to be a vector of [-2,2] uniform
arma::vec init_theta(int m);
// Auxiliary 10 : sample k random integers in 0:(n-1)
arma::uvec sample_int(int n, int k);
#endif