// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"

using namespace Rcpp;
using namespace arma;

/////////////////////////////////////////////////////////////////////////////////////////////////
// Export Functions
// [[Rcpp::export]]
Rcpp::List exec_SAMCoptim(Function func, const int nv, arma::vec& energy, arma::mat& domain,
                     const int niter, arma::vec& vecpi, double t0, const double xi, arma::vec stepsize,
                     arma::vec& init, double tau_h, double tau_star, double tau_t0,
                     double trunc_init, double trunc_r){
  
  // 1-1. setting
  arma::mat samples(niter,nv,fill::zeros); // (n-by-p) convention
  const int m = vecpi.n_elem;              // number of energy domains
  arma::vec visits(m,fill::zeros);         // visiting frequency
  arma::vec thetas = init_theta(m);        // weights to be undated

  // 1-2. variable settings
  arma::vec xold = init;
  double Hy,Hxold;   // energy
  int    Jy,Jxold;   // location
  double ttxold,tty; // according theta values
  double r;          // acceptance
  
  // 1-3. xold-related computations
  Hxold  = sum(as<NumericVector>(func(xold))); // need to wrap it for the pointer part.
  if ((Hxold < energy.min())||(Hxold > energy.max())){
    stop("* SAMC backend : Hxold is not within the partition range.");
  }
  Jxold  = find_location(Hxold,energy);
  ttxold = thetas(Jxold);
  
  // 1-4. main iteration
  double gfactor = 0.0;
  arma::vec Ijxold(m,fill::zeros);
  double Mbound = trunc_init;
  double temp_now = 0.0;
  for (int i=0;i<niter;i++){
    // additional ; computing temperature now
    temp_now = tau_h*(std::sqrt(static_cast<float>(tau_t0/std::max(static_cast<double>(i+1),tau_t0)))) + tau_star;
    
    // A-1. sample generation
    arma::vec y = sampling_rwvec(xold,domain,stepsize);
    // A-2. compute ratio
    Hy  = sum(as<NumericVector>(func(y)));
    Jy  = find_location(Hy,energy);
    tty = thetas(Jy);
    r = static_cast<double>(std::exp(static_cast<float>((-(Hy-Hxold)/temp_now)+(ttxold-tty))));
    // A-3. accept or reject
    if (R::runif(0,1)<=std::min(1.0,r)){ // Accepted
      samples.row(i) = y.t();
      visits(Jy) += 1;
      xold = y;
      Hxold = Hy;
      Jxold = Jy;
      ttxold = tty;
    } else {
      samples.row(i) = xold.t();
      visits(Jxold) += 1;
    }
    // B-1. temporary update
    arma:vec thalf(m,fill::zeros);
    gfactor = std::pow(t0/(std::max(t0,static_cast<double>(i+1))),xi);
    
    Ijxold.fill(0.0);
    Ijxold(Jxold) = 1;
    for (int j=0;j<m;j++){
      thalf(j) = thetas(j) + gfactor*(Ijxold(j)-vecpi(j));
    }
    // B-2. Accept the update of weights or not
    if (arma::norm(thalf,2)<Mbound){
      thetas = thalf;
    } else {
      Mbound = Mbound*trunc_r;
      thetas = init_theta(m);
    }
  }
  
  // 1-5. return output
  Rcpp::List output;
  output["samples"] = samples;
  output["frequency"] = visits;
  output["theta"] = thetas;
  
  return(output);
}