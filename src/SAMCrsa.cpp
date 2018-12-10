// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"

using namespace Rcpp;
using namespace arma;

////////////////////////////////////////////////////////////////////////
double init_double(double min, double max){ // initialize with [min, max]
  double hey = arma::randu();
  hey = (max-min)*hey + min;
  return(hey);
}
double initlog_double(double min, double max){
  double hey = arma::randu();
  hey = (max-min)*hey + min;
  return(std::exp(hey));
}
arma::vec init_vec(int n, double min, double max){ // initialize a vector in [min, max]
  arma::vec hey(n,fill::zeros);
  for (int i=0;i<n;i++){
    hey(i) = (arma::randu()*(max-min)) + min;
  }
  return(hey);
}
arma::vec initlog_vec(int n, double min, double max){
  arma::vec hey(n,fill::zeros);
  for (int i=0;i<n;i++){
    hey(i) = std::exp((arma::randu()*(max-min)) + min);
  }
  return(hey);
}
arma::mat compute_Rz(arma::mat distmat, double phi){ // compute Rz(i,j)
  int n = distmat.n_rows;
  arma::mat output(n,n,fill::zeros);
  for (int i=0;i<n;i++){
    output(i,i) = 1.0;
  }
  double tmpval;
  for (int i=0;i<(n-1);i++){
    for (int j=(i+1);j<n;j++){
      tmpval = std::exp(-distmat(i,j)/phi);
      output(i,j) = tmpval;
      output(j,i) = tmpval;
    }
  }
  return(output);
}
arma::mat compute_dRzPhi(arma::mat distmat, double phi){ // compute dRz/dPhi
  int n = distmat.n_rows;
  arma::mat output(n,n,fill::zeros); // diagonal is now zero h_ij*exp(-h_ij) = 0*exp(-0) = 0*1 = 0
  double inval, outval, tmpval;
  for (int i=0;i<(n-1);i++){
    for (int j=(i+1);j<n;j++){
      outval = distmat(i,j)/(phi*phi);
      inval  = std::exp(-distmat(i,j)/phi);
      tmpval = outval*inval;
      
      output(i,j) = tmpval;
      output(j,i) = tmpval;
    }
  }
  
  return(output);
}
arma::vec compute_muz(double beta0, arma::vec betas, arma::mat X){
  arma::vec tmp = X*betas;
  for (int i=0;i<(tmp.n_elem);i++){
    tmp(i) = tmp(i) + beta0;
  }
  return(tmp);
}
double compute_increment(double beta0old, double beta0new, arma::vec betasold, arma::vec betasnew, double phiold, double phinew,
                         double sigma2old, double sigma2new, double tau2old, double tau2new){
  double inc1 = static_cast<double>(std::pow(static_cast<float>(beta0old-beta0new),2.0));
  double inc2 = static_cast<double>(std::pow(static_cast<float>(arma::norm(betasold-betasnew,2)),2.0));
  double inc3 = (phiold-phinew)*(phiold-phinew);
  double inc4 = (sigma2old-sigma2new)*(sigma2old-sigma2new);
  double inc5 = static_cast<double>(std::pow(static_cast<float>(tau2old-tau2new),2.0));
  double output = (static_cast<double>(std::sqrt(static_cast<float>(inc1+inc2+inc3+inc4+inc5))));
  return(output);
}
bool check_trunc(double beta0, arma::vec betas, double phi, double sig2, double tau2, int trunc){
  double dpi = static_cast<double>(trunc);
  bool cond1 = ((beta0<=(2.0+dpi))&&(beta0>=(-2.0-dpi)));
  bool cond2 = ((all(betas <= (2.0+dpi)))&&(all(betas>=(-2.0-dpi))));
  bool cond3 = ((phi<=std::exp(4.0+dpi))&&(phi>=std::exp(-dpi)));
  bool cond4 = ((sig2<=std::exp(2.0+dpi))&&(sig2>=std::exp(-2.0-dpi)));
  bool cond5 = ((tau2<=std::exp(2.0+dpi))&&(tau2>=std::exp(-2.0-dpi)));
  
  bool final;
  if (cond1&&cond2&&cond3&&cond4&&cond5){
    final = true;
  } else {
    final = false;
  }
  return(final);
}

/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List ksSAMCrsa(arma::mat coords, arma::vec Y, arma::mat X, int nsubset, int maxiter, double a0, double t0, double b0, double eta){
  // hooo
  double ee = 2.7182818284590452353602874;
  double ee2 = ee*ee;
  double ee4 = ee2*ee2;
  
  // preliminary : get parameters
  int n = coords.n_rows;
  int p = X.n_cols; 
  
  // preliminary : pairwise Euclidean distance among coordinates
  int i,j;
  arma::rowvec coord1(2,fill::zeros);
  arma::rowvec coord2(2,fill::zeros);
  arma::mat distcoords(n,n,fill::zeros);
  for (i=0;i<(n-1);i++){
    coord1 = coords.row(i);
    for (j=(i+1);j<n;j++){
      coord2 = coords.row(j);
      double distval = arma::norm(coord1-coord2,2);
      distcoords(i,j) = distval;
      distcoords(j,i) = distval;
    }
  }
  
  // initialize
  double beta0old = init_double(-2.0,2.0);
  double beta0new = 0.0;
  arma::vec betasold = init_vec(p, -2.0, 2.0);
  arma::vec betasnew(p,fill::zeros);
  double phiold = initlog_double(-1.0,4.0);
  double phinew = 0.0;
  double sig2old = initlog_double(-2.0,2.0);
  double sig2new = 0.0;
  double tau2old = initlog_double(-2.0,2.0);
  double tau2new = 0.0;
  
  // main computation
  int iter = 0;                             // iteration counter
  arma::uvec idpick(nsubset, fill::zeros);  // index for subsamples
  int truncation = 0;                       // pi_t : varying truncation counter
  
  arma::vec muz(nsubset, fill::zeros);              // \mu_z : partial empirical mean
  arma::vec subY(nsubset, fill::zeros);             // partial Y
  arma::mat subcoords(nsubset, 2, fill::zeros);     // partial coords
  arma::mat subX(nsubset, p, fill::zeros);          // partial X
  arma::mat subdist(nsubset, nsubset, fill::zeros); // partial distance matrix
  
  arma::mat Rz(nsubset, nsubset, fill::zeros);        // Rz : partial correlation
  arma::mat dRzPhi(nsubset, nsubset, fill::zeros);    // dRzPhi : dRz/dPhi : derivative
  arma::mat Sigmaz(nsubset, nsubset, fill::zeros);    // Sigmaz
  arma::mat invSigmaz(nsubset, nsubset, fill::zeros); // invSigmaz
  
  double a_t = 0.0; // a_t : gain factor
  double b_t = 0.0; // b_t : norm bound
  
  arma::vec record_beta0(maxiter,fill::zeros);
  arma::mat record_betas(maxiter,p,fill::zeros);
  arma::vec record_phi(maxiter,fill::zeros);
  arma::vec record_sig2(maxiter,fill::zeros);
  arma::vec record_tau2(maxiter,fill::zeros);
  
  while (iter < maxiter){
    // 0. record !
    record_beta0(iter) = beta0old;
    record_betas.row(iter) = betasold;
    record_phi(iter) = phiold;
    record_sig2(iter) = sig2old;
    record_tau2(iter) = tau2old;
    
    
    // 1. pick a random vector and choose subs
    idpick = sample_int(n, nsubset);
    subY = Y.rows(idpick);
    subX = X.rows(idpick);
    subcoords = coords.rows(idpick);
    subdist = distcoords.submat(idpick, idpick);
    muz  = compute_muz(beta0old, betasold, subX);
    
    // 2. compute Rz, dRzPhi, Sigmaz and inverse
    Rz     = compute_Rz(subdist, phiold);
    dRzPhi = compute_dRzPhi(subdist, phiold);
    Sigmaz = sig2old*Rz;
    for (int i=0;i<nsubset;i++){
      Sigmaz(i,i) = Sigmaz(i,i) + tau2old;
    }
    invSigmaz = arma::pinv(Sigmaz);
    
    // 3. update parameters
    a_t = a0*t0/static_cast<double>(std::max(static_cast<double>(iter),t0));
    beta0new = beta0old + (a_t*arma::as_scalar(arma::sum(invSigmaz*(subY-muz))));
    betasnew = betasnew + (a_t*(subX.t()*invSigmaz*(subY-muz)));
    phinew   = phiold   + a_t*((arma::sum(arma::diagvec(sig2old*invSigmaz*dRzPhi))/(-2.0))+((sig2old/2.0)*arma::as_scalar((subY.t()-muz.t())*invSigmaz*dRzPhi*invSigmaz*(subY-muz))));
    sig2new  = sig2old  + a_t*(((arma::sum(arma::diagvec(invSigmaz*Rz)))/(-2.0)) + (0.5* arma::as_scalar((subY.t()-muz.t())*invSigmaz*Rz*invSigmaz*(subY-muz))));
    tau2new  = tau2old  + a_t*((arma::sum(arma::diagvec(invSigmaz))/(-2.0)) + (0.5*arma::as_scalar((subY.t()-muz.t())*invSigmaz*invSigmaz*(subY-muz))));
    
    // 4. varying truncation
    b_t = b0*(static_cast<double>(std::pow(static_cast<float>(t0/std::max(static_cast<double>(iter),t0)),eta)));
    
    // //******************** SIMPLY ASSUME SO LARGE ***************************************//
    // beta0old = beta0new;
    // betasold = betasnew;
    // phiold   = phinew;
    // sig2old  = sig2new;
    // tau2old  = tau2new;
    
    if ((compute_increment(beta0old,beta0new,betasold,betasnew,phiold,phinew,sig2old,sig2new,tau2old,tau2new) < b_t)&&(check_trunc(beta0new,betasnew,phinew,sig2new,tau2new,truncation))){
    //if ((compute_increment(beta0old,beta0new,betasold,betasnew,phiold,phinew,sig2old,sig2new,tau2old,tau2new) < b_t)){
      // successfully accept and update
      beta0old = beta0new;
      betasold = betasnew;
      phiold   = phinew;
      sig2old  = sig2new;
      tau2old  = tau2new;
    } else {
      // varying truncation IN CHARGE !!!!!!!!!!!!!!!!!!!!!!!!!!! - yeah, i'm devastated.
      beta0old = init_double(-2.0,2.0);
      betasold = init_vec(p, -2.0, 2.0);
      phiold  = initlog_double(-1.0, 4.0);
      sig2old = initlog_double(-2.0, 2.0);
      tau2old = initlog_double(-2.0, 2.0);
      truncation += 1;
    }
    
    // update iteation
    iter += 1;
  }
  
  List hey;
  hey["beta0"] = beta0old;
  hey["betas"] = betasold;
  hey["phi"]  = phiold;
  hey["sig2"] = sig2old;
  hey["tau2"] = tau2old;
  
  hey["ntrunc"] = truncation;
  hey["rec_beta0"] = record_beta0;
  hey["rec_betas"] = record_betas;
  hey["rec_phi"]   = record_phi;
  hey["rec_sig2"]  = record_sig2;
  hey["rec_tau2"]  = record_tau2;
  return(hey);
}