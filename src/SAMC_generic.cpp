// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"

using namespace Rcpp;
using namespace arma;


/////////////////////////////////////////////////////////////////////////////////////////////////
// DENSITY type
template <typename T>
List core_samc_density(T func,const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                       const int niter, arma::vec& vecpi, const double t0, const double xi,
                       arma::vec stepsize, arma::mat& trange, arma::vec& init) {
  // 1-1. setting
  arma::mat samples(niter,nv,fill::zeros); // (n-by-p) convention
  const int m = vecpi.n_elem;              // number of energy domains
  arma::vec visits(m,fill::zeros);         // visiting frequency
  arma::vec thetas(m,fill::zeros);         // weights to be undated

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
  for (int i=0;i<niter;i++){
    
    // A-1. sample generation
    arma::vec y = sampling_rwvec(xold,domain,stepsize);
    // A-2. compute ratio
    Hy  = sum(as<NumericVector>(func(y)));
    Jy  = find_location(Hy,energy);
    tty = thetas(Jy);
    r = exp((-(Hy-Hxold)/tau)+(ttxold-tty));
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
    double gfactor = t0/(std::max(t0,pow((double)(i),xi)));
    arma::vec Ijxold(m,fill::zeros);
    for (int j=0;j<m;j++){
      Ijxold(j) = 0;
    }
    Ijxold(Jxold) = 1;
    for (int j=0;j<m;j++){
      thalf(j) = thetas(j) + gfactor*(Ijxold(j)-vecpi(j));
    }
    // B-2. Accept the update of weights or not
    if ((all(thalf<=trange.col(1)))&&(all(trange.col(0)<=thalf))){
      thetas = thalf;
    } else {
      thetas = adjust_weights(thalf,trange);
    }
  }

  // 1-5. return output
  Rcpp::List output;
  output["samples"] = samples;
  output["frequency"] = visits;
  output["theta"] = thetas;
  
  return(output);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Export Functions
// [[Rcpp::export]]
Rcpp::List exec_SAMC(Function func, const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                               const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize,
                               arma::mat& trange,  arma::vec& init){
  return core_samc_density<Function>(func,nv,energy,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init);
}
