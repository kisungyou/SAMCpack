## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache=FALSE, autodep=TRUE, cache.comments=FALSE,
  comment = "#>"
)

## ------------------------------------------------------------------------
myoption = list()
myoption$partition = c(-Inf,seq(from=-8,to=0,length.out=41)) # energy partition
myoption$tau       = 1.0                                     # temperature
myoption$domain    = c(-1.1,1.1)                             # domain for sample space
myoption$vecpi     = as.vector(rep(1/41,41))                 # desired sampling distribution
myoption$niter     = 200000                                  # total number of iterations
myoption$stepsize  = 0.25                                    # s.d for random-walk proposal

## ------------------------------------------------------------------------
func_r = function(x){
x1 = x[1]; x2 = x[2];
val1 = (-(x1*sin(20*x2)+x2*sin(20*x1))^2)*cosh(sin(10*x1)*x1);
val2 = (-(x1*cos(10*x2)-x2*sin(10*x1))^2)*cosh(cos(20*x2)*x2);
return(val1+val2);
}

## ------------------------------------------------------------------------
library(RcppXPtrUtils)
cppscript = "SEXP funcH(arma::vec x){
double x1 = x(0);
double x2 = x(1);
double val1 = (-std::pow((x1*sin(20*x2)+x2*sin(20*x1)),2))*cosh(sin(10*x1)*x1);
double val2 = (-std::pow((x1*cos(10*x2)-x2*sin(10*x1)),2))*cosh(cos(20*x2)*x2);
return Rcpp::wrap(val1+val2);
}"
func_ptr = RcppXPtrUtils::cppXPtr(cppscript,depends="RcppArmadillo") # as a pointer

## ------------------------------------------------------------------------
library(SAMCpack)
res1 = SAMC(2,func_r,options=myoption)       # use R function
res2 = SAMCPLUS(2,func_ptr,options=myoption) # use C++ function

## ----echo=FALSE----------------------------------------------------------
select = seq(from=101,to=myoption$niter,by=100) # 100 burn-in, 1/100 thinning 
par(mfrow=c(1,2))
plot(res1$samples[select,1], res1$samples[select,2],xlab='x',ylab='y',main='R function samples')
plot(res2$samples[select,1], res2$samples[select,2],xlab='x',ylab='y',main='C++ function samples')

## ------------------------------------------------------------------------
library(microbenchmark)
execution = microbenchmark(
  list = alist(
    runR   = SAMC(2,func_r,options=myoption),
    runCPP = SAMCPLUS(2,func_ptr,options=myoption)
  ), times=10
)

## ---- echo=FALSE---------------------------------------------------------
boxplot(execution, xlab="")

## ------------------------------------------------------------------------
data_vec = as.vector(c(10,20)) # vector-valued data=(10,20)
cppscript_vec = "SEXP funcH_vec(arma::vec x, arma::vec data){
double x1 = x(0);
double x2 = x(1);
double par1 = data(0);  // first element  : 10
double par2 = data(1);  // second element : 20
double val1 = (-std::pow((x1*sin(par2*x2)+x2*sin(par2*x1)),2))*cosh(sin(par1*x1)*x1);
double val2 = (-std::pow((x1*cos(par1*x2)-x2*sin(par1*x1)),2))*cosh(cos(par2*x2)*x2);
return Rcpp::wrap(val1+val2);
}"

## ------------------------------------------------------------------------
data_mat = matrix(c(20,10,20,10,10,20),nrow=2) 
cppscript_mat = "SEXP funcH_mat(arma::vec x, arma::mat data){
double x1 = x(0);
double x2 = x(1);
double val1 = (-std::pow((x1*sin(data(0,0)*x2)+x2*sin(data(0,1)*x1)),2))*cosh(sin(data(0,2)*x1)*x1);
double val2 = (-std::pow((x1*cos(data(1,0)*x2)-x2*sin(data(1,1)*x1)),2))*cosh(cos(data(1,2)*x2)*x2);
return Rcpp::wrap(val1+val2);
}"

## ------------------------------------------------------------------------
func_ptrvec = RcppXPtrUtils::cppXPtr(cppscript_vec,depends="RcppArmadillo") # as a pointer
func_ptrmat = RcppXPtrUtils::cppXPtr(cppscript_mat,depends="RcppArmadillo") # as a pointer

## ------------------------------------------------------------------------
res3 = SAMCPLUS(2,func_ptrvec,options=myoption,data=data_vec) # vector data
res4 = SAMCPLUS(2,func_ptrmat,options=myoption,data=data_mat) # matrix data

## ----echo=FALSE----------------------------------------------------------
select = seq(from=101,to=myoption$niter,by=100) # 100 burn-in, 1/100 thinning 
par(mfrow=c(1,2))
plot(res3$samples[select,1], res3$samples[select,2],xlab='x',ylab='y',main='arma::vec data')
plot(res4$samples[select,1], res4$samples[select,2],xlab='x',ylab='y',main='arma::mat data')

## ------------------------------------------------------------------------
data_list = list()
data_list[[1]] = as.vector(c(20,10))
data_list[[2]] = matrix(c(20,10,10,20),nrow=2)

## ------------------------------------------------------------------------
cppscript_list = "SEXP funcH_list(arma::vec x, Rcpp::List datalist){
double x1 = x(0);
double x2 = x(1);

NumericVector lefty = as<NumericVector>(datalist[0]);
NumericMatrix right = as<NumericMatrix>(datalist[1]);

double val1 = (-std::pow((x1*sin(lefty[0]*x2)+x2*sin(right(0,0)*x1)),2))*cosh(sin(right(0,1)*x1)*x1);
double val2 = (-std::pow((x1*cos(lefty[1]*x2)-x2*sin(right(1,0)*x1)),2))*cosh(cos(right(1,1)*x2)*x2);
return Rcpp::wrap(val1+val2);
}"

## ----echo=FALSE----------------------------------------------------------
func_ptrlist = RcppXPtrUtils::cppXPtr(cppscript_list,depends="RcppArmadillo") # as a pointer
res5 = SAMCPLUS(2,func_ptrlist,options=myoption,data=data_list)

## ----echo=FALSE----------------------------------------------------------
 select = seq(from=101,to=myoption$niter,by=100) # 100 burn-in, 1/100 thinning 
 par(mfrow=c(1,2))
 plot(res5$samples[select,1], res5$samples[select,2],xlab='x',ylab='y',main='with List Data')
 barplot(as.vector(res5$frequency/sum(res5$frequency)),
         main="visiting frequency by energy partition",
         names.arg=myoption$partition[-1], xlab="energy")

## ------------------------------------------------------------------------
cppscript_sexp = "SEXP funcH_sexp(arma::vec x, SEXP data){
double x1 = x(0);
double x2 = x(1);

// Wrap the input 'data' as 'datavec' (Rcpp type)
NumericVector datavec = as<NumericVector>(data);

double par1 = datavec[0];  // first element  : 10
double par2 = datavec[1];  // second element : 20
double val1 = (-std::pow((x1*sin(par2*x2)+x2*sin(par2*x1)),2))*cosh(sin(par1*x1)*x1);
double val2 = (-std::pow((x1*cos(par1*x2)-x2*sin(par1*x1)),2))*cosh(cos(par2*x2)*x2);
return Rcpp::wrap(val1+val2);
}"

## ------------------------------------------------------------------------
data_sexp = as.vector(c(10,20)) # vector-valued data=(10,20) to be passed as SEXP
func_ptrsexp = cppXPtr(cppscript_sexp,depends="RcppArmadillo")
res6 = SAMCPLUS(2,func_ptrsexp,options=myoption,data=data_sexp)

## ----echo=FALSE----------------------------------------------------------
 select = seq(from=101,to=myoption$niter,by=100) # 100 burn-in, 1/100 thinning 
 par(mfrow=c(1,2))
 plot(res6$samples[select,1], res6$samples[select,2],xlab='x',ylab='y',main='with SEXP data')
 barplot(as.vector(res6$frequency/sum(res5$frequency)),
         main="visiting frequency by energy partition",
         names.arg=myoption$partition[-1], xlab="energy")

