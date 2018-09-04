#' SAMC Sampler with C++ 
#' 
#' @examples 
#' \dontrun{
#' ##### (1) Two-Dimensional Multimodal sampling
#' ## Define negative log-density function as CPP function
#' cppscript = "SEXP funcH(arma::vec x){
#' double x1 = x(0);
#' double x2 = x(1);
#' double val1 = (-std::pow((x1*sin(20*x2)+x2*sin(20*x1)),2))*cosh(sin(10*x1)*x1);
#' double val2 = (-std::pow((x1*cos(10*x2)-x2*sin(10*x1)),2))*cosh(cos(20*x2)*x2);
#' return Rcpp::wrap(val1+val2);
#' }"
#' func_ptr = RcppXPtrUtils::cppXPtr(cppscript,depends="RcppArmadillo") # as a pointer
#' 
#' ## Other Experimental Setting
#' ex_part   = c(-Inf,seq(from=-8,to=0,by=0.2))
#' m         = length(ex_part)-1
#' ex_temp   = 1.0
#' ex_step   = (0.25)
#' ex_vecpi  = as.vector(array(1/m,c(1,m)))
#' ex_domain = c(-1.1,1.1)
#' ex_niter  = 200000
#' 
#' ## Run The Code
#' res = SAMCplus(2,func_ptr,partition=ex_part,tau=ex_temp,stepsize=ex_step,vecpi=ex_vecpi,domain=ex_domain,niter=ex_niter)
#' 
#' ## Visualize
#' select = seq(from=101,to=ex_niter,by=100) # 100 burn-in, 1/100 thinning 
#' par(mfrow=c(1,2))
#' plot(res$samples[select,1], res$samples[select,2],xlab='x',ylab='y',main='samples')
#' barplot(as.vector(res$frequency/sum(res$frequency)),
#'         main="visiting frequency by energy partition",
#'         names.arg=ex_part[2:(m+1)], xlab="energy")
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
#' res2 = SAMCplus(2,func_ptr2,partition=ex_part,tau=ex_temp,stepsize=ex_step,vecpi=ex_vecpi,domain=ex_domain,niter=ex_niter,data=vecdata)
#' select = seq(from=101,to=ex_niter,by=100) # 100 burn-in, 1/100 thinning 
#' par(mfrow=c(1,2))
#' plot(res2$samples[select,1], res2$samples[select,2],xlab='x',ylab='y',main='samples')
#' barplot(as.vector(res2$frequency/sum(res2$frequency)),
#'         main="visiting frequency by energy partition",
#'         names.arg=ex_part[2:(m+1)], xlab="energy")
#' }
#' 
#' @author Kisung You
#' @rdname SAMCplus_generic
#' @export
SAMCplus <- function(nv,energy,domain=c(-Inf,Inf),partition=seq(-1e+2,1e+2,length.out=9),vecpi=rep(1/10,10),
                     tau=1.0,niter=20000,t0=200,xi=2/3,stepsize=1.0,trange=c(-Inf,Inf),data=NA){
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
  #   9. (double)  stepsize : normal proposal
  #  10. (vec/mat) trange   : theta space of (m-by-2)
  if ((length(nv)>1)||(!is.numeric(nv))||(nv<1)||(is.infinite(nv))||(is.na(nv))){
    stop("* SAMCplus : 'nv' should be a positive integer.")
  }
  nv = as.integer(nv)
  if (is.vector(domain)){
    if ((length(domain)!=2)||(any(!is.numeric(domain)))||(any(is.na(domain)))){
      stop("* SAMCplus : 'domain' should be a vector of length 2.")
    }
    domain = matrix(sort(domain),nrow=nv,ncol=length(domain),byrow=TRUE)
  } else if (is.matrix(domain)){
    if ((dim(domain)[1]!=nv)||(dim(domain)[2]!=2)||(any(is.na(domain)))||(any(!is.numeric(domain)))){
      stop("* SAMCplus : 'domain' should be a matrix of size (nv-by-2).")
    }
    for (i in 1:nv){
      domain[i,] = sort(domain[i,])
    }
  } else {
    stop("* SAMCplus : 'domain' should be either a vector or a matrix.")
  }
  if (any(is.infinite(domain))){
    message("* SAMCplus : 'domain' would be better without infinite values.")
    message("*      : automatically replacing Inf's with suitably large numbers.")
    domain = adjust_inf(domain)
  }
  if ((!is.vector(partition))||(length(partition)<2)||(any(!is.numeric(partition)))){
    stop("* SAMCplus : 'energy' should be a vector whose length is greater than 1.")
  }
  partition = as.vector(sort(partition))
  if ((any(is.null(vecpi)))||(length(vecpi)!=(length(partition)-1))||(any(vecpi<=0))||(abs(sum(vecpi)-1)>1e-10)){
    stop("* SAMCplus : desired sampling distribution 'vecpi' is invalid.")
  }
  vecpi = as.vector(sort(vecpi))
  tau = as.double(tau)
  if ((!is.numeric(niter))||(niter<=1)||(is.infinite(niter))||(is.na(niter))||(length(niter)>1)){
    stop("* SAMCplus : 'niter' should be a positive integer as an iteration number.")
  }
  t0 = as.double(t0)
  if ((!is.numeric(xi))||(xi<=0.5)||(xi>1)||(length(xi)>1)||(is.null(xi))){
    stop("* SAMCplus : 'xi' should be in (0.5,1].")
  }
  xi = as.double(xi)
  if ((!is.numeric(stepsize))||(stepsize<=0)||(is.infinite(stepsize))||(length(stepsize)>1)){
    stop("* SAMCplus : 'stepsize' is a standard deviation term for normal proposal density")
  }
  stepsize = as.double(stepsize)
  m = length(vecpi)
  if (is.vector(trange)){
    if ((length(trange)!=2)||(any(!is.numeric(trange)))||(any(is.na(trange)))){
      stop("* SAMCplus : 'trange' should be a vector of length 2.")
    }
    trange = matrix(sort(trange),nrow=m,ncol=length(trange),byrow=TRUE)
  } else if (is.matrix(trange)){
    if ((dim(trange)[1]!=m)||(dim(trange)[2]!=2)||(any(is.na(trange)))||(any(!is.numeric(trange)))){
      stop("* SAMCplus : 'trange' should be a matrix of size (m-by-2).")
    }
    for (i in 1:m){
      trange[i,] = sort(trange[i,])
    }
  } else {
    stop("* SAMCplus : 'trange' should be either a vector or a matrix.")
  }
  
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
  if ((length(data)==1)&&(is.na(data))){
    checkXPtr(energy,"SEXP",c("arma::vec"))
    output = exec_samcfast_type0(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init)
  } else if (is.vector(data)){
    checkXPtr(energy,"SEXP",c("arma::vec","arma::vec"))
    output = exec_samcfast_type1(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data)
  } else if (is.matrix(data)){
    checkXPtr(energy,"SEXP",c("arma::vec","arma::mat"))
    output = exec_samcfast_type2(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data)
  } else if (is.list(data)){
    checkXPtr(energy,"SEXP",c("arma::vec","Rcpp::List"))
    output = exec_samcfast_type3(energy,nv,partition,domain,tau,niter,vecpi,t0,xi,stepsize,trange,init,data)  
  } else {
    stop("* SAMCplus : 'data' should be one of NA, vector, matrix, or list.")
  }
      
      
  
  
  ##-------------------------------------------------------------------
  # END
  return(output)
}