// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"

using namespace Rcpp;
using namespace arma;

// auxiliary functions
Rcpp::List cpp_ttest_resample(arma::vec x, arma::vec y, int L){
  int nx = x.n_elem;
  int ny = y.n_elem;
  
  Rcpp::List uvecs_nx = two_perm_vec(nx, L);
  Rcpp::List uvecs_ny = two_perm_vec(ny, L);
  
  arma::uvec idx1 = as<arma::uvec>(uvecs_nx["vec1"]);
  arma::uvec idx2 = as<arma::uvec>(uvecs_nx["vec2"]);
  arma::uvec idy1 = as<arma::uvec>(uvecs_ny["vec1"]);
  arma::uvec idy2 = as<arma::uvec>(uvecs_ny["vec2"]);
  
  arma::vec xsample = arma::join_cols(x.rows(idx2), y.rows(idy1));
  arma::vec ysample = arma::join_cols(x.rows(idx1), y.rows(idy2));
  
  Rcpp::List output;
  output["xsample"] = xsample;
  output["ysample"] = ysample;
  return(output);
}
// [[Rcpp::export]]
double cpp_ttest_statistic(arma::vec x, arma::vec y, bool var_equal){
  double nx = static_cast<double>(x.n_elem);
  double ny = static_cast<double>(y.n_elem);
  double output = 0.0;
  if (var_equal==false){
    output = std::abs(static_cast<float>(arma::mean(x)-arma::mean(y)))/std::sqrt(((arma::var(x)/nx) + (arma::var(y)/ny)));
  } else {
    double term1 = std::abs(static_cast<float>(arma::mean(x)-arma::mean(y)));
    double term2 = std::sqrt(((nx-1.0)*arma::var(x) + (ny-1.0)*arma::var(y))/(nx+ny-2.0))*std::sqrt((1.0/nx)+(1.0/ny));
    output = term1/term2;
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List samc_cpp_ttest(arma::vec x, arma::vec y, const int niter, int L, int m, int t0,
                          arma::vec frequency_old, arma::vec partition_old,
                          arma::vec thetas_old, arma::vec vecpi_old,
                          bool var_equal){
  // pseudo-initialization
  double count_accept = 0;
  
  arma::vec frequency = frequency_old;
  arma::vec partition = partition_old;
  arma::vec thetas    = thetas_old;
  arma::vec vecpi     = vecpi_old;
  
  double ts_init = cpp_ttest_statistic(x, y, var_equal);
  double ts_old  = ts_init;
  double ts_new  = 0.0;
  
  int id_old = find_mingeq(partition, ts_init);
  int id_new = 0;

  // main iteration
  //   initialization
  Rcpp::List resampler;
  arma::vec x_old = x;
  arma::vec y_old = y;
  arma::vec x_new;
  arma::vec y_new;
  double r, tmpr;
  double runif1;
  double gfactor = 0.0;
  for (int b=0;b<niter;b++){
    // step 1. update proposal
    resampler = cpp_ttest_resample(x, y, L);
    x_new = as<arma::vec>(resampler["xsample"]);
    y_new = as<arma::vec>(resampler["ysample"]);
    
    // step 2. compute test statistic, partition, and acceptance probability
    ts_new = cpp_ttest_statistic(x_new, y_new, var_equal);
    id_new = find_mingeq(partition, ts_new);
    tmpr   = std::exp(static_cast<float>(thetas(id_old) - thetas(id_new)));
    if (tmpr <= 1.0){
      r = tmpr;
    } else {
      r = 1.0;
    }
    
    // step 3. dichotomous branching for accept/reject
    runif1 = R::runif(0.0, 1.0);
    if (runif1 <= r){
      count_accept = count_accept + 1.0;
      x_old = x_new;
      y_old = y_new;
      ts_old = ts_new;
      id_old = id_new;
    }
    
    if (b>t0){
      gfactor = static_cast<double>(t0)/static_cast<double>(b);
    } else {
      gfactor = static_cast<double>(t0)/static_cast<double>(t0);
    }
    arma::vec id_update(m, fill::zeros);
    id_update(id_old) = 1.0;
    thetas = thetas + gfactor*(id_update-vecpi);
    
    // step 4. frequency update
    frequency(id_old) += 1;
  }
  
  // return really many things
  Rcpp::List output;
  output["thetas"] = thetas;
  output["frequency"] = frequency;
  output["count.accept"] = count_accept;
  output["ts.init"] = ts_init;
  return(output);
}