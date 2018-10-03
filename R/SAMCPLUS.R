#' SAMC Sampler with C++ 
#' 
#' The function \code{SAMCPLUS} is a generic SAMC sampler for distributions on continuous domain. Instead of an \code{R} function, 
#' \code{SAMCPLUS} requires a function pointer to be provided for faster sampling, with all other values and parameters being equal 
#' to its cousin \code{\link{SAMC}}. We limited the flexibility of the function pointer to be passed. See the below for more details or 
#' the vignette.
#' 
#' @section Note on writing your own C++ function:
#' First, the output should be returned as \code{SEXP} rather than \code{double} in evaluating the negative log density. Second, the variable and 
#' extra data should be provided as \code{arma::vec} and \code{arma::mat} type, with an exception for \code{Rcpp::List} for list-valued data. This means, for the \code{data} Even 
#' though we could let data to be passed freely, we believe using  \href{https://cran.r-project.org/web/packages/RcppArmadillo/index.html}{RcppArmadillo}, which is a templated linear algebra library, 
#' enables easier writing of one's own C++ code in a style of R or MATLAB while providing sufficient computational power. Furthermore, limiting extra data to one of 3 types (vector, matrix, and list) 
#' reduces potential type-matching issue in encapsulating of the current environment by removing unexpected errors a user might have incurred.
#' 
#' @param nv number of variables.
#' @param energy a \code{CPP} function pointer for negative log density.
#' @param data extra data to be supplemented. It should be a vector, a matrix, or a list. 
#' @param options a list specifying parameters/options for SAMC algorithm,
#' \tabular{lll}{
#' PARAMETER        \tab SPECIFICATION \tab DESCRIPTION \cr
#' \code{domain}    \tab vector(\eqn{2}) or matrix(\eqn{(nv\times 2)}) \tab domain of sample space \cr
#' \code{partition} \tab vector(\eqn{m}) \tab energy partition \cr
#' \code{vecpi}     \tab vector(\eqn{m-1}) \tab desired sampling distribution \cr
#' \code{tau}       \tab positive number \tab temperature \cr
#' \code{niter}     \tab positive integer \tab number of iterations to be run \cr
#' \code{t0}        \tab \eqn{(0.5,1]}  \tab gain factor sequence value \cr
#' \code{xi}        \tab positive number \tab gain factor sequence exponent \cr
#' \code{stepsize}  \tab positive number \tab stepsize for random-walk sampler \cr
#' \code{trange}    \tab vector(\eqn{2}) or matrix(\eqn{m\times 2}) \tab domain of estimates for \eqn{\log(g_i /\pi_i)}
#' }
#' 
#' @return a named list containing \describe{
#' \item{samples}{an \eqn{(niter\times nv)} samples generated.}
#' \item{frequency}{length-\eqn{m} vector of visiting frequency for energy partition.}
#' \item{theta}{length-\eqn{m} vector of estimates of \eqn{\log(g_i / \pi_i)}}
#' }
#' 
#' @examples 
#' \dontrun{
#' ##### Two-Dimensional Multimodal sampling
#' ## Step 1 : Define negative log-density function as a CPP function
#' cppscript = "SEXP funcH(arma::vec x){
#' double x1 = x(0);
#' double x2 = x(1);
#' double val1 = (-std::pow((x1*sin(20*x2)+x2*sin(20*x1)),2))*cosh(sin(10*x1)*x1);
#' double val2 = (-std::pow((x1*cos(10*x2)-x2*sin(10*x1)),2))*cosh(cos(20*x2)*x2);
#' return Rcpp::wrap(val1+val2);
#' }"
#' func_ptr = RcppXPtrUtils::cppXPtr(cppscript,depends="RcppArmadillo") # as a pointer
#' 
#' ## Step 2 : Prepare a setting option
#' myoption = list()
#' myoption$partition = c(-Inf,seq(from=-8,to=0,length.out=41))
#' myoption$tau       = 1.0
#' myoption$domain    = c(-1.1,1.1)
#' myoption$vecpi     = as.vector(rep(1/41,41))
#' myoption$niter     = 200000
#' myoption$stepsize  = 0.25
#' 
#' ## Step 3 : Run The Code
#' res = SAMCPLUS(2,func_ptr,options=myoption)
#' 
#' ## Step 4 : Visualize
#' select = seq(from=101,to=myoption$niter,by=100) # 100 burn-in, 1/100 thinning 
#' par(mfrow=c(1,2))
#' plot(res$samples[select,1], res$samples[select,2],xlab='x',ylab='y',main='samples')
#' barplot(as.vector(res$frequency/sum(res$frequency)),
#'         main="visiting frequency by energy partition",
#'         names.arg=myoption$partition[-1], xlab="energy")
#'
#' ##### (2) Use Extra Data
#' ## Define negative log-density function as CPP function
#' cppscript2 = "SEXP funcH2(arma::vec x, arma::vec data){
#' double x1 = x(0);
#' double x2 = x(1);
#' double par1 = data(0);
#' double par2 = data(1);
#' 
#' double val1 = (-std::pow((x1*sin(par2*x2)+x2*sin(par2*x1)),2))*cosh(sin(par1*x1)*x1);
#' double val2 = (-std::pow((x1*cos(par1*x2)-x2*sin(par1*x1)),2))*cosh(cos(par2*x2)*x2);
#' return Rcpp::wrap(val1+val2);
#' }"
#' func_ptr2 = RcppXPtrUtils::cppXPtr(cppscript2,depends="RcppArmadillo") # as a pointer
#' 
#' ## Run The Code
#' vecdata = as.vector(c(10,20)) 
#' res2 = SAMCPLUS(2,func_ptr2,data=vecdata, options=myoption)
#' select = seq(from=101,to=ex_niter,by=100) # 100 burn-in, 1/100 thinning 
#' par(mfrow=c(1,2))
#' plot(res2$samples[select,1], res2$samples[select,2],xlab='x',ylab='y',main='samples')
#' barplot(as.vector(res2$frequency/sum(res2$frequency)),
#'         main="visiting frequency by energy partition",
#'         names.arg=ex_part[2:(m+1)], xlab="energy")
#' }
#' 
#' @references 
#' \insertRef{SAMC}{SAMCpack}
#' 
#' @author Kisung You
#' @export
SAMCPLUS <- function(nv,energy,data=NA,options=list()){
  ##-------------------------------------------------------------------
  # PREPROCESSING
  #   1. (int)     nv        : number of variables
  #   2. (func)    energy    : NOT CHECKING FOR CPP 
  #   3. (vec/mat) domain    : for samples to reside (nv-by-2) or length-2 vector
  #   3. (vec)     partition : energy-level partition for dividing the sample space into partitions (m+1)
  #   4. (vec)     vecpi     : desired sampling distribution (m)
  #   5. (double)  tau       : temperature
  #   6. (int)     niter     : number of iterations for sampling
  #   7. (double)  t0       : gain factor sequence
  #   8. (double)  xi       : gain factor exponent
  #   9. (num/vec)  stepsize : normal proposal
  #  10. (vec/mat) trange   : theta space of (m-by-2)
  if (!is.list(options)){
    stop("* SAMCPLUS : 'options' must be a list.")
  }
  
  pruned    = optionlist(options,nv,"SAMCPLUS")
  domain    = pruned$domain
  partition = pruned$partition
  vecpi     = pruned$vecpi
  tau       = pruned$tau
  niter     = pruned$niter
  t0        = pruned$t0
  xi        = pruned$xi
  stepsize  = pruned$stepsize
  trange    = pruned$trange
  m         = pruned$m  # added for fun
  # parnames = names(options)
  # if (!("domain" %in% parnames)){    domain = c(-Inf,Inf)  } else {domain=options$domain}
  # if (!("partition" %in% parnames)){ partition = seq(-1e+2,1e+2,length.out=9)} else {partition=options$partition}
  # if (!("vecpi" %in% parnames)){     vecpi=rep(1/10,10)} else {vecpi=options$vecpi}
  # if (!("tau" %in% parnames)){       tau=1.0} else {tau=options$tau}
  # if (!("niter" %in% parnames)){     niter=20000} else {niter=options$niter}
  # if (!("t0" %in% parnames)){        t0=200} else {t0=options$t0}
  # if (!("xi" %in% parnames)){        xi=2/3} else {xi=options$xi}
  # if (!("stepsize" %in% parnames)){  stepsize=1.0} else {stepsize=options$stepsize}
  # if (!("trange" %in% parnames)){    trange=c(-Inf,Inf)} else {trange=options$trange}
  # 
  # 
  # if ((length(nv)>1)||(!is.numeric(nv))||(nv<1)||(is.infinite(nv))||(is.na(nv))){
  #   stop("* SAMCplus : 'nv' should be a positive integer.")
  # }
  # nv = as.integer(nv)
  # if (is.function(energy)){
  #   stop("* SAMCplus : for any R-compiled function, use SAMC instead.")
  # }
  # if (is.vector(domain)){
  #   if ((length(domain)!=2)||(any(!is.numeric(domain)))||(any(is.na(domain)))){
  #     stop("* SAMCplus : 'domain' should be a vector of length 2.")
  #   }
  #   domain = matrix(sort(domain),nrow=nv,ncol=length(domain),byrow=TRUE)
  # } else if (is.matrix(domain)){
  #   if ((dim(domain)[1]!=nv)||(dim(domain)[2]!=2)||(any(is.na(domain)))||(any(!is.numeric(domain)))){
  #     stop("* SAMCplus : 'domain' should be a matrix of size (nv-by-2).")
  #   }
  #   for (i in 1:nv){
  #     domain[i,] = sort(domain[i,])
  #   }
  # } else {
  #   stop("* SAMCplus : 'domain' should be either a vector or a matrix.")
  # }
  # if (any(is.infinite(domain))){
  #   message("* SAMCplus : 'domain' would be better without infinite values.")
  #   message("*      : automatically replacing Inf's with suitably large numbers.")
  #   domain = adjust_inf(domain)
  # }
  # if ((!is.vector(partition))||(length(partition)<2)||(any(!is.numeric(partition)))){
  #   stop("* SAMCplus : 'energy' should be a vector whose length is greater than 1.")
  # }
  # partition = as.vector(sort(partition))
  # if ((any(is.null(vecpi)))||(length(vecpi)!=(length(partition)-1))||(any(vecpi<=0))||(abs(sum(vecpi)-1)>1e-10)){
  #   stop("* SAMCplus : desired sampling distribution 'vecpi' is invalid.")
  # }
  # vecpi = as.vector(sort(vecpi))
  # tau = as.double(tau)
  # if ((!is.numeric(niter))||(niter<=1)||(is.infinite(niter))||(is.na(niter))||(length(niter)>1)){
  #   stop("* SAMCplus : 'niter' should be a positive integer as an iteration number.")
  # }
  # t0 = as.double(t0)
  # if ((!is.numeric(xi))||(xi<=0.5)||(xi>1)||(length(xi)>1)||(is.null(xi))){
  #   stop("* SAMCplus : 'xi' should be in (0.5,1].")
  # }
  # xi = as.double(xi)
  # if ((!is.numeric(stepsize))||(stepsize<=0)||(is.infinite(stepsize))||(length(stepsize)>1)){
  #   stop("* SAMCplus : 'stepsize' is a standard deviation term for normal proposal density")
  # }
  # stepsize = as.double(stepsize)
  # m = length(vecpi)
  # if (is.vector(trange)){
  #   if ((length(trange)!=2)||(any(!is.numeric(trange)))||(any(is.na(trange)))){
  #     stop("* SAMCplus : 'trange' should be a vector of length 2.")
  #   }
  #   trange = matrix(sort(trange),nrow=m,ncol=length(trange),byrow=TRUE)
  # } else if (is.matrix(trange)){
  #   if ((dim(trange)[1]!=m)||(dim(trange)[2]!=2)||(any(is.na(trange)))||(any(!is.numeric(trange)))){
  #     stop("* SAMCplus : 'trange' should be a matrix of size (m-by-2).")
  #   }
  #   for (i in 1:m){
  #     trange[i,] = sort(trange[i,])
  #   }
  # } else {
  #   stop("* SAMCplus : 'trange' should be either a vector or a matrix.")
  # }
  
  ##-------------------------------------------------------------------
  # INITIALIZATION
  init = array(0,nv)
  for (i in 1:nv){
    subdomain = domain[i,]
    if (is.infinite(subdomain[1])){
      if (is.infinite(subdomain[2])){
        init[i] = 0
      }else{
        init[i] = runif(1,subdomain[2]-10,subdomain[2])
      }
    } else{
      if (is.infinite(subdomain[2])){
        init[i] = runif(1,subdomain[1],subdomain[1]+10)
      }else{
        init[i] = (subdomain[1]+subdomain[2])/2
      }
    }
  }
  
  ##-------------------------------------------------------------------
  # MAIN COMPUTATION : here I'm case-branching the data
  sexpflag = TRUE
  sexpflag = tryCatch(checkXPtr(energy,"SEXP",c("arma::vec","SEXP")),
                     error = function(e){FALSE})
  if (is.null(sexpflag)){ # case we have SEXP data input
    output = exec_samcfast_sexpdata(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data)
  } else if (sexpflag==FALSE){ # specified cases
    if ((is.list(data))&&(is.vector(data))){
      checkXPtr(energy,"SEXP",c("arma::vec","Rcpp::List"))
      output = exec_samcfast_type3(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data)  
    } else if ((is.na(data))||(any(is.na(data)))){
      checkXPtr(energy,"SEXP",c("arma::vec"))
      output = exec_samcfast_type0(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init)
    } else if ((is.vector(data))&&(!is.list(data))){
      checkXPtr(energy,"SEXP",c("arma::vec","arma::vec"))
      output = exec_samcfast_type1(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data)
    } else if (is.matrix(data)){
      checkXPtr(energy,"SEXP",c("arma::vec","arma::mat"))
      output = exec_samcfast_type2(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data)
    }   
  }
  
  
  ##-------------------------------------------------------------------
  # END
  return(output)
}