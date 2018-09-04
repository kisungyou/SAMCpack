// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"

using namespace Rcpp;
using namespace arma;

/*
 * Type0 : core_samefast_type0 : No Data
 * Type1 : core_samefast_type1 : arma::vec data
 * Type2 : core_samefast_type2 : arma::mat data
 * Type3 : core_samefast_type3 : Rcpp::List data
 */

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Type0 : core_samefast_type0 : No Data
template <typename T>
List core_samcfast_type0(T func,const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                         const int niter, arma::vec& vecpi, const double t0, const double xi,
                         const double stepsize, arma::mat& trange, arma::vec& init) {
  // 2-1. setting
  arma::mat samples(niter,nv,fill::zeros); // (n-by-p) convention
  const int m = vecpi.n_elem;              // number of energy domains
  arma::vec visits(m,fill::zeros);         // visiting frequency
  arma::vec thetas(m,fill::zeros);         // weights to be undated
  
  // 2-2. variable settings
  arma::vec xold = init;
  double Hy,Hxold;   // energy
  int    Jy,Jxold;   // location
  double ttxold,tty; // according theta values
  double r;          // acceptance
  
  // 2-3. xold-related computations
  Hxold  = sum(as<NumericVector>(func(xold))); // need to wrap it for the pointer part.
  Jxold  = find_location(Hxold,energy);
  ttxold = thetas(Jxold);
  
  // 2-4. main iteration
  for (int i=0;i<niter;i++){
    // A-1. sample generation
    arma::vec y = sampling_rw(xold,domain,stepsize);
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
  
  // 2-5. return output
  Rcpp::List output;
  output["samples"] = samples;
  output["frequency"] = visits;
  output["theta"] = thetas;
  
  return(output);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
// Type1 : core_samefast_type1 : arma::vec data
template <typename T>
List core_samcfast_type1(T func,const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                         const int niter, arma::vec& vecpi, const double t0, const double xi,
                         const double stepsize, arma::mat& trange, arma::vec& init, arma::vec data) {
  // 2-1. setting
  arma::mat samples(niter,nv,fill::zeros); // (n-by-p) convention
  const int m = vecpi.n_elem;              // number of energy domains
  arma::vec visits(m,fill::zeros);         // visiting frequency
  arma::vec thetas(m,fill::zeros);         // weights to be undated
  
  // 2-2. variable settings
  arma::vec xold = init;
  double Hy,Hxold;   // energy
  int    Jy,Jxold;   // location
  double ttxold,tty; // according theta values
  double r;          // acceptance
  
  // 2-3. xold-related computations
  Hxold  = sum(as<NumericVector>(func(xold,data))); // need to wrap it for the pointer part.
  Jxold  = find_location(Hxold,energy);
  ttxold = thetas(Jxold);
  
  // 2-4. main iteration
  for (int i=0;i<niter;i++){
    // A-1. sample generation
    arma::vec y = sampling_rw(xold,domain,stepsize);
    // A-2. compute ratio
    Hy  = sum(as<NumericVector>(func(y,data)));
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
  
  // 2-5. return output
  Rcpp::List output;
  output["samples"] = samples;
  output["frequency"] = visits;
  output["theta"] = thetas;
  
  return(output);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
// Type2 : core_samefast_type2 : arma::mat data
template <typename T>
List core_samcfast_type2(T func,const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                         const int niter, arma::vec& vecpi, const double t0, const double xi,
                         const double stepsize, arma::mat& trange, arma::vec& init, arma::mat data) {
  // 2-1. setting
  arma::mat samples(niter,nv,fill::zeros); // (n-by-p) convention
  const int m = vecpi.n_elem;              // number of energy domains
  arma::vec visits(m,fill::zeros);         // visiting frequency
  arma::vec thetas(m,fill::zeros);         // weights to be undated
  
  // 2-2. variable settings
  arma::vec xold = init;
  double Hy,Hxold;   // energy
  int    Jy,Jxold;   // location
  double ttxold,tty; // according theta values
  double r;          // acceptance
  
  // 2-3. xold-related computations
  Hxold  = sum(as<NumericVector>(func(xold,data))); // need to wrap it for the pointer part.
  Jxold  = find_location(Hxold,energy);
  ttxold = thetas(Jxold);
  
  // 2-4. main iteration
  for (int i=0;i<niter;i++){
    // A-1. sample generation
    arma::vec y = sampling_rw(xold,domain,stepsize);
    // A-2. compute ratio
    Hy  = sum(as<NumericVector>(func(y,data)));
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
  
  // 2-5. return output
  Rcpp::List output;
  output["samples"] = samples;
  output["frequency"] = visits;
  output["theta"] = thetas;
  
  return(output);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
// Type3 : core_samefast_type3 : Rcpp::List data
template <typename T>
List core_samcfast_type3(T func,const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                         const int niter, arma::vec& vecpi, const double t0, const double xi,
                         const double stepsize, arma::mat& trange, arma::vec& init, Rcpp::List data) {
  // 2-1. setting
  arma::mat samples(niter,nv,fill::zeros); // (n-by-p) convention
  const int m = vecpi.n_elem;              // number of energy domains
  arma::vec visits(m,fill::zeros);         // visiting frequency
  arma::vec thetas(m,fill::zeros);         // weights to be undated
  
  // 2-2. variable settings
  arma::vec xold = init;
  double Hy,Hxold;   // energy
  int    Jy,Jxold;   // location
  double ttxold,tty; // according theta values
  double r;          // acceptance
  
  // 2-3. xold-related computations
  Hxold  = sum(as<NumericVector>(func(xold,data))); // need to wrap it for the pointer part.
  Jxold  = find_location(Hxold,energy);
  ttxold = thetas(Jxold);
  
  // 2-4. main iteration
  for (int i=0;i<niter;i++){
    // A-1. sample generation
    arma::vec y = sampling_rw(xold,domain,stepsize);
    // A-2. compute ratio
    Hy  = sum(as<NumericVector>(func(y,data)));
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
  
  // 2-5. return output
  Rcpp::List output;
  output["samples"] = samples;
  output["frequency"] = visits;
  output["theta"] = thetas;
  
  return(output);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
// Export Functions
// Pointer Definition
typedef SEXP (*ptrtype0)(arma::vec);
typedef SEXP (*ptrtype1)(arma::vec, arma::vec);
typedef SEXP (*ptrtype2)(arma::vec, arma::mat);
typedef SEXP (*ptrtype3)(arma::vec, Rcpp::List);

// [[Rcpp::export]]
Rcpp::List exec_samcfast_type0(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                                   const int niter, arma::vec& vecpi, const double t0, const double xi, const double stepsize,
                                   arma::mat& trange,  arma::vec init){
  ptrtype0 func = *XPtr<ptrtype0>(func_);
  return core_samcfast_type0<ptrtype0>(func,nv,energy,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init);
}
// [[Rcpp::export]]
Rcpp::List exec_samcfast_type1(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                               const int niter, arma::vec& vecpi, const double t0, const double xi, const double stepsize,
                               arma::mat& trange,  arma::vec init, arma::vec data){
  ptrtype1 func = *XPtr<ptrtype1>(func_);
  return core_samcfast_type1<ptrtype1>(func,nv,energy,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data);
}
// [[Rcpp::export]]
Rcpp::List exec_samcfast_type2(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                               const int niter, arma::vec& vecpi, const double t0, const double xi, const double stepsize,
                               arma::mat& trange,  arma::vec init, arma::mat data){
  ptrtype2 func = *XPtr<ptrtype2>(func_);
  return core_samcfast_type2<ptrtype2>(func,nv,energy,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data);
}
// [[Rcpp::export]]
Rcpp::List exec_samcfast_type3(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau,
                               const int niter, arma::vec& vecpi, const double t0, const double xi, const double stepsize,
                               arma::mat& trange,  arma::vec init, Rcpp::List data){
  ptrtype3 func = *XPtr<ptrtype3>(func_);
  return core_samcfast_type3<ptrtype3>(func,nv,energy,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data);
}