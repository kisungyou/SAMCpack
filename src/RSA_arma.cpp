// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <stdlib.h>

using namespace Rcpp;

arma::mat padded_pinv(arma::mat X){
  const int N = X.n_rows;
  arma::mat smallX = X.submat(1,1,N-1,N-1);
  arma::mat smallY = arma::inv(smallX);
  
  // pad column of N-1 first
  arma::colvec colpad(N-1); colpad.fill(0);
  arma::mat fatY = arma::join_horiz(colpad, smallY);
  
  // now pad row of N elemenst
  arma::rowvec rowpad(N); rowpad.fill(0);
  arma::mat Y = arma::join_vert(rowpad, fatY);
  return(Y);
}
arma::ivec arma_subset_sample(int M,int N)
{
  int i, j, k, m;
  double rnormval;
  
  arma::ivec x(M+1); x.fill(0);
  arma::ivec y(N+1); y.fill(0);
  m=N;
  for(i=1; i<=N; i++) y(i)=i;
  for(i=1; i<=M; i++){
    j=0;
    while(j<1 || j>m){
      GetRNGstate();
      rnormval = norm_rand();
      PutRNGstate();
      j=std::floor(rnormval*1.0/RAND_MAX*m)+1;
    } 
    x(i) = y(j);
    for(k=j+1; k<=m; k++) y(k-1)=y(k);
    m--;
  }
  
  return(x);
  
  return 0;
}



// [[Rcpp::export]]
Rcpp::List RSAarma(arma::vec pData,int pDataCol,int pDataNum, int pSampleNum,int pStepscale, int pTotal_Iteration, int pWarm){
  // pre-header
  const double pi = 3.14159265;
  const double rho = 0.01;
  const double eta = 1.0;
  const double Rho = 100;
  const double Eta = 0.55;
  //const double RAND_MAX = 65536;
  
  // [039-044] Parameters
  int sampleNum = pSampleNum;
  int dataCol = pDataCol;
  int dataNum = pDataNum;
  int stepscale=pStepscale;
  int total_iteration=pTotal_Iteration;
  int warm = pWarm;
  
  // Initialize : Others
  int i, j;
  double avephi, avesigmasq, avetausq;
  int iter, truncation, truncount;
  double delta, sum, a, delctrl;
  double phi, kappa, sigmasq, tausq;
  double Dphi, Dsigmasq, Dtausq;//Dkappa,
  
  // [048-051] Initialize : data loading
  arma::mat ppData(pDataNum+1, pDataCol+1); ppData.fill(0);
  for(i = 1;i<=pDataNum;i++){
    for(j = 1;j<=pDataCol;j++)//{
      ppData(i,j) = pData((j-1)*(pDataNum)+(i-1));
  }
  
  // [062-077] initialization : variable setting
  arma::ivec sam(sampleNum+1); sam.fill(0);
  arma::dmat R(sampleNum+1,sampleNum+1); R.fill(0.0);
  arma::dmat V(sampleNum+1,sampleNum+1); V.fill(0.0);
  arma::dmat IV(sampleNum+1,sampleNum+1); IV.fill(0.0);
  arma::dmat Dist(sampleNum+1,sampleNum+1); Dist.fill(0.0);
  arma::dmat DRphi(sampleNum+1,sampleNum+1); DRphi.fill(0.0);
  arma::dmat DRkappa(sampleNum+1,sampleNum+1); DRkappa.fill(0.0);
  arma::dvec mu(sampleNum+1); mu.fill(0.0);
  arma::dmat C(sampleNum+1,dataCol+1); C.fill(0.0);
  arma::dvec zmu(sampleNum+1); zmu.fill(0.0);
  arma::dvec b(sampleNum+1); b.fill(0.0);
  
  arma::dvec Dbeta(dataCol-2+1); Dbeta.fill(0.0);
  arma::dvec beta(dataCol-2+1); beta.fill(0.0);
  arma::dvec avebeta(dataCol-2+1); avebeta.fill(0.0);
  
  // [080-082] initialization : parameters
  arma::vec low(6); low.fill(0.0);
  arma::vec up(6);  up.fill(0.0);
  arma::vec sze(6); sze.fill(0.0);
  arma::vec lower(6); lower.fill(0.0);
  arma::vec upper(6); upper.fill(0.0);
  arma::vec add(pDataCol+2); add.fill(0.0); ///////////////// 1, or 2 to be run ?
  
  low(1)=-2.0; up(1)=2.0; low(2)=-2.0; up(2)=2.0;
  low(3)=0.0;  up(3)=4.0; low(4)=-2.0; up(4)=2.0; low(5)=-2.0; up(5)=2.0;
  sze(1)=1.0; sze(2)=1.0; sze(3)=1.0; sze(4)=1.0; sze(5)=1.0;
  avephi=avesigmasq=avetausq=0.0;
  
  arma::ivec vec1N = arma::regspace<arma::ivec>(1,1,dataNum);
  
  for (j=1;j<=dataCol-2;j++){
    avebeta(j) = 0;
  }
  truncount  = 0;
  truncation = 1;
  
  double rnormval;
  GetRNGstate();
  rnormval = norm_rand();
  PutRNGstate();
  double logphi, logsigmasq, logtausq;//logkappa, 
  
  // Algorithm Iteration  
  for(iter=1; iter<=total_iteration; iter++){
    ABC:
    if(truncation==1){ 
      for(j=1;j<=dataCol-2;j++){
        GetRNGstate();
        rnormval = norm_rand();
        PutRNGstate();
        beta(j)=rnormval*1.0/RAND_MAX*(up(1)-low(1))+low(1);
      }
      GetRNGstate();
      rnormval = norm_rand();
      PutRNGstate();
      logphi=rnormval*1.0/RAND_MAX*(up(3)-low(3))+low(3);

      GetRNGstate();
      rnormval = norm_rand();
      PutRNGstate();
      logsigmasq=rnormval*1.0/RAND_MAX*(up(4)-low(4))+low(4); 
      
      GetRNGstate();
      rnormval = norm_rand();
      PutRNGstate();
      logtausq=rnormval*1.0/RAND_MAX*(up(5)-low(5))+low(5);
      lower(1)=low(1)-sze(1)*truncount; upper(1)=up(1)+sze(1)*truncount;
      lower(2)=low(2)-sze(2)*truncount; upper(2)=up(2)+sze(2)*truncount;
      lower(3)=low(3)-sze(3)*truncount; upper(3)=up(3)+sze(3)*truncount;
      lower(4)=low(4)-sze(4)*truncount; upper(4)=up(4)+sze(4)*truncount;
      lower(5)=low(5)-sze(5)*truncount; upper(5)=up(5)+sze(5)*truncount;
    }
    
    
    if (iter<=stepscale){
      delta=rho; 
      delctrl=Rho;}
    else{
      delta=rho*exp(eta*std::log(1.0*stepscale/iter));
      delctrl=Rho*exp(Eta*std::log(1.0*stepscale/iter));
    } 
    
    // Generate a subset sample of size M from the set 1:N
    // arma::ivec zero1(1); zero1.fill(0);
    // sam = arma::join_vert(zero1,RcppArmadillo::sample(vec1N, sampleNum, TRUE, NumericVector::create()));
    // 
    sam = arma_subset_sample(sampleNum, dataNum);
    phi=std::exp(logphi); kappa=1.0; sigmasq=std::exp(logsigmasq); tausq=std::exp(logtausq);
    
    // V=sigmasq*R+ tausq*I 
    for(i=1; i<=sampleNum; i++){
      R(i,i)=1.0; Dist(i,i)=0.0; V(i,i)=sigmasq+tausq;DRphi(i,i)=0.0; DRkappa(i,i)=0.0;
      for(j=1; j<i; j++){
        Dist(j,i)=std::sqrt(std::pow(ppData(sam(i),1)-ppData(sam(j),1),2.0)+std::pow(ppData(sam(i),2)-ppData(sam(j),2),2.0));
        Dist(i,j)=Dist(j,i);
        
        a=std::pow(Dist(i,j)/phi,kappa);
        R(j,i)=std::exp(-a); R(i,j)=R(j,i);
        V(j,i)=R(i,j)*sigmasq; V(i,j)=V(j,i);
        DRphi(j,i)=a*R(i,j)*kappa/phi; DRphi(i,j)=DRphi(j,i);
        DRkappa(j,i)=-a*R(i,j)*std::log(Dist(i,j)/phi); DRkappa(i,j)=DRkappa(j,i);
      }
    } 
    // this must be changed : now, V /////////////////////////////////////////////////////////////
    IV = padded_pinv(V);
    // IV = padded_pinv(V); // matrix_inverse(V,IV,sampleNum) with Y=inv(X), return d=log(det(X))
    
    for(i=1; i<=sampleNum; i++){
      for(j=1;j<=dataCol-2;j++){
        C(i,j)=ppData(sam(i),j+2);
      } 
      C(i,1)=1;
      mu(i)=0;
      for(j=1;j<=dataCol-2;j++){
        mu(i)+=beta(j)*C(i,j);
      }
      zmu(i)=ppData(sam(i),3)-mu(i);
    }
    
    for(i=1; i<=sampleNum; i++){
      b(i)=0.0;
      for(j=1; j<=sampleNum; j++){
        b(i)+=IV(i,j)*zmu(j);
      }
    }
    
    
    // calculate derivatives    
    for(j=1;j<=dataCol-2;j++){
      Dbeta(j)=0.0;
      for (i=1;i<=sampleNum;i++){
        Dbeta(j)+=b(i)*C(i,j);
      }
    }
    
    Dphi=0.0;
    for(i=1; i<=sampleNum; i++){
      for(j=1; j<=sampleNum; j++){
        Dphi+=b(i)*DRphi(i,j)*b(j);
      }
    }
    for(sum=0.0, i=1; i<=sampleNum; i++){
      for(j=1; j<=sampleNum; j++){
        sum+=IV(i,j)*DRphi(j,i); 
      }
    }
    Dphi-=sum;
    Dphi*=0.5*sigmasq;
    
    Dsigmasq=0.0;
    for(i=1; i<=sampleNum; i++){
      for(j=1; j<=sampleNum; j++){
        Dsigmasq+=b(i)*R(i,j)*b(j); 
      }
    }
    for(sum=0.0, i=1; i<=sampleNum; i++){
      for(j=1; j<=sampleNum; j++){
        sum+=IV(i,j)*R(j,i); 
      }
    }
    Dsigmasq-=sum;
    Dsigmasq*=0.5; 
    
    Dtausq=0.0;
    for(i=1; i<=sampleNum; i++){
      Dtausq+=b(i)*b(i); 
    }
    for(sum=0.0, i=1; i<=sampleNum; i++){
      sum+=IV(i,i);
    }
    Dtausq-=sum;
    Dtausq*=0.5;
    
    
    for(j=1;j<=dataCol-2;j++){
      add(j)=delta*Dbeta(j);  
    }
    add(dataCol-1)=delta*Dphi*phi;
    add(dataCol)=delta*Dsigmasq*sigmasq;
    add(dataCol+1)=delta*Dtausq*tausq;
    
    //beta0+=add[1];
    //beta1+=add[2];
    for(j=1;j<=dataCol-2;j++){
      beta(j)+=add(j);
    }
    logphi+=add(dataCol-1);
    logsigmasq+=add(dataCol);
    logtausq+=add(dataCol+1);
    
    truncation=0;
    //for(sum=0.0, i=1; i<=5; i++) sum+=add[i]*add[i];
    
    for(sum=0.0, j=1; j<=dataCol+1; j++){
      sum+=add(j)*add(j); 
    }
    if(std::sqrt(sum)>delctrl){
      truncation=1; 
    } else {
      for(j=1;j<=dataCol-2;j++){
        if (beta(j)<lower(1) || beta(j)>upper(1)) truncation=1;
        else if(logphi<lower(3) || logphi>upper(3)) truncation=1;
        else if(logsigmasq<lower(4) || logsigmasq>upper(4)) truncation=1;
        else if(logtausq<lower(5) || logtausq>upper(5)) truncation=1; 
      }
    }
    
    if(truncation==1){ truncount++; 
      goto ABC; 
    }
    if(iter>warm){ avephi+=logphi/(total_iteration-warm); avesigmasq+=logsigmasq/(total_iteration-warm); avetausq+=logtausq/(total_iteration-warm); 
    for(j=1;j<=dataCol-2;j++) avebeta(j)+=beta(j)/(total_iteration-warm);
    }
  }
  
  arma::vec pbeta(dataCol-2); pbeta.fill(0);
  for(j=1;j<=dataCol-2;j++) pbeta(j-1)=beta(j);
  //    *pbeta0=beta0; 
  //    *pbeta1=beta1; 
  double pPhi     = std::exp(logphi);
  double pSigmasq = std::exp(logsigmasq);
  double pTausq   = std::exp(logtausq);
  // 
  Rcpp::List output;
  output["pbeta"] = pbeta;
  output["pPhi"]  = pPhi;
  output["pSigmasq"] = pSigmasq;
  output["pTausq"]   = pTausq;
  return(output);
}