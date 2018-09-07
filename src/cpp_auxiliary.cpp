#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// 1. Sampling from random-walk proposal
arma::vec sampling_rw(arma::vec xold, arma::mat domain, const double stepsize){
  // 1-1. checker
  const int n = xold.n_elem;
  arma::vec xnew(n,fill::zeros);
  
  // 1-2. apply random walk
  xnew = xold + stepsize*arma::vec(n,fill::randn);
  if ((all(xnew<=domain.col(1)))&&(all(domain.col(0)<=xnew))){
    return(xnew);
  } else{
    return(xold);
  }
}
arma::vec sampling_rwvec(arma::vec xold, arma::mat domain, arma::vec stepsize){
  // 1-1. checker
  const int n = xold.n_elem;
  arma::vec xnew(n,fill::zeros);
  arma::vec randomness(n,fill::randn);
  
  // 1-2. apply random walk
  for (int i=0;i<n;i++){
    xnew(i) = xold(i) + stepsize(i)*randomness(i);
  }
  if ((all(xnew<=domain.col(1)))&&(all(domain.col(0)<=xnew))){
    return(xnew);
  } else{
    return(xold);
  }
}



// Auxiliary 2 : find a suitable location for energy vector
int find_location(double x, arma::vec y){
  // 2-1. checker
  const int n = y.n_elem;
  // 2-2. find it
  int location;
  for (int i=0;i<n;i++){
    if ((x<=y(i+1))&&(x>=y(i))){
      location = i;
      break;
    }
  }
  return(location);
}

// Auxiliary 3 : adjust weight update
arma::vec adjust_weights(arma::vec& weight,arma::mat& trange){
  // 3-1. checker
  const int n = weight.n_elem;
  arma::vec output(n,fill::zeros);
  double disc = 0;
  
  // 3-2. adjust discrepancies
  if (any(weight<trange.col(0))){ // left boundary
    for (int i=0;i<n;i++){
      double tmpdisc = trange(i,0)-weight(i);
      if (tmpdisc>disc){
        disc = tmpdisc;
      }
    }
    for (int i=0;i<n;i++){
      output(i) = weight(i)+disc;
    }
  } else if (any(trange.col(1)<weight)){ // right boundary
    for (int i=0;i<n;i++){
      double tmpdisc = weight(i)-trange(i,1);
      if (tmpdisc>disc){
        disc = tmpdisc;
      }
    }
    for (int i=0;i<n;i++){
      output(i) = weight(i) - disc;
    }
  }
  return(output);
}

// Auxiliary 4 : find min, max
arma::vec find_minmax(arma::vec x){
  arma::vec values(2,fill::zeros);
  values(0) = x.min();
  values(1) = x.max();
  return(values);
}
arma::vec find_rowminmax(arma::rowvec x){
  arma::vec values(2,fill::zeros);
  values(0) = x.min();
  values(1) = x.max();
  return(values);
}

// Auxiliary 5 : rescale : vertically and horizontally
// 5-1. simply run by vector
arma::colvec rescale_colvec(arma::colvec x){
  // parameters and ready
  const int n = x.n_elem;
  arma::colvec y((2*n)-1,fill::zeros);
  // original values
  for (int i=0;i<n;i++){
    y(2*i) = x(i);
  }
  // intermediate values
  for (int j=0;j<(n-1);j++){
    y(2*j+1) = (y(2*j)+y(2*(j+1)))/2;
  }
  return(y);
}
arma::rowvec rescale_rowvec(arma::rowvec x){
  // parameters and ready
  const int n = x.n_elem;
  arma::rowvec y((2*n)-1,fill::zeros);
  // original values
  for (int i=0;i<n;i++){
    y(2*i) = x(i);
  }
  // intermediate values
  for (int j=0;j<(n-1);j++){
    y(2*j+1) = (y(2*j)+y(2*(j+1)))/2;
  }
  return(y);
}
// 5-2. vertically
// [[Rcpp::export]]
arma::mat rescale_vert2(arma::mat A){
  // parameters
  const int m = A.n_rows;
  const int n = A.n_cols;
  // get ready
  arma::mat B((2*m)-1,n,fill::zeros);
  for (int i=0;i<n;i++){
    B.col(i) = rescale_colvec(A.col(i));
  }
  return(B);
}
// 5-3. horizontally
// [[Rcpp::export]]
arma::mat rescale_hori2(arma::mat A){
  const int m = A.n_rows;
  const int n = A.n_cols;
  arma::mat B(m,(2*n)-1,fill::zeros);
  arma::vec tmprow;
  arma::vec tmpcol;
  for (int i=0;i<m;i++){
    B.row(i) = rescale_rowvec(A.row(i));
  }
  return(B);
}

// Auxiliary 6 : get grid vector, whow!
arma::vec get_gridvec(const double pstart, const double pend, const int N){
  arma::vec output = linspace<vec>(pstart,pend,N);
  return(output);
}

// Auxiliary 7 : evaludate CSAMC2 g(lambda(x))
double evalcsamc_g2(arma::mat g, arma::vec gvec1, arma::vec gvec2, arma::vec lbdx){
  int idx1, idx2;
  for (idx1=0;idx1<gvec1.n_elem;idx1++){
    if ((lbdx(0)>=gvec1(idx1))&&(lbdx(0)<=gvec1(idx1+1))){
      break;
    }
  }
  for (idx2=0;idx2<gvec2.n_elem;idx2++){
    if ((lbdx(1)>=gvec2(idx2))&&(lbdx(1)<=gvec2(idx2+1))){
      break;
    }
  }
  
  double gij = g(idx1,idx2);
  double gij1= g(idx1,idx2+1);
  double gi1j = g(idx1+1,idx2);
  double gi1j1= g(idx1+1,idx2+1);
  
  double u = (lbdx(0)-gvec1(idx1))/(gvec1(idx1+1)-gvec1(idx1));
  double v = (lbdx(1)-gvec2(idx2))/(gvec2(idx2+1)-gvec2(idx2));
  
  double gval=((1-u)*(1-v)*gij)+((1-u)*v*gij1) + (u*(1-v)*gi1j) + (u*v*gi1j1);
  return(gval);
}

// Auxiliary 8 : kernel density update for 2dimensional case
arma::mat evalcsamc_ks2(arma::vec vec1, arma::vec vec2, arma::mat ysamples, arma::mat H){
  // 8-1. get parameters and ready
  const int L1 = vec1.n_elem;
  const int L2 = vec2.n_elem;
  const int M  = ysamples.n_cols;
  arma::mat zeta(L1,L2,fill::zeros);
  
  double detH = det(H);
  
  arma::colvec tmpz(2,fill::zeros);
  arma::colvec x(2,fill::zeros);
  // 8-2. iterative updates
  for (int i=0;i<L1;i++){
    for (int j=0;j<L2;j++){
      tmpz(0) = vec1(i);
      tmpz(1) = vec2(j);
      
      for (int it=0;it<M;it++){
        x = tmpz-ysamples.col(it);
        if (norm(x)<0.01){
          zeta(i,j) += exp(-dot(x,solve(H,x))/2)/sqrt(detH);
        }
      }
    }
  }
  // 8-3. normalization
  double allsum = accu(zeta);
  zeta = zeta/allsum;
  // 8-4. return output
  return(zeta);
}