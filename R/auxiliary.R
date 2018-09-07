# COMMON FUNCTIONS --------------------------------------------------------
#   1. adjust_inf
#   2. check_inclusion
#   3. optionlist
#   4. dist_dinvgamma
#   5. vec_init
#   6. compute_pinv

#' @keywords internal
#' @noRd
adjust_inf <- function(pdomain){
  IdxPInf = which(pdomain==Inf)
  IdxNInf = which(pdomain==-Inf)
  IdxConst= which(!is.infinite(pdomain))
  if (length(IdxConst)==0){
    output = matrix(rep(c(-10000,10000),nrow(pdomain)),nrow=nrow(pdomain),byrow=TRUE)
  } else {
    rangep  = range(pdomain[IdxConst])
    rgapad  = abs(rangep[2]-rangep[1])/10
    
    maxvalue = max(pdomain[IdxConst])
    minvalue = min(pdomain[IdxConst])
    
    output = pdomain
    output[IdxPInf] = maxvalue+rgapad
    output[IdxNInf] = minvalue-rgapad
  }
  return(output)
}

#' @keywords internal
#' @noRd
check_inclusion <- function(dvec,pvec){
  dvec = as.vector(dvec)
  pvec = as.vector(pvec)
  
  if ((dvec[1]<=pvec[1])&&(pvec[2]<=dvec[2])){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# OPTIONLIST covers 9 parameters
#   - 'domain','partition','vecpi','tau','niter','t0','xi','stepsize','trange'
#' @keywords internal
#' @noRd
optionlist <- function(options, nv, fname){
  parnames = names(options)
  if (!("domain" %in% parnames)){    domain = c(-Inf,Inf)  } else {domain=options$domain}
  if (!("partition" %in% parnames)){ partition = seq(-1e+2,1e+2,length.out=10)} else {partition=options$partition}
  if (!("vecpi" %in% parnames)){     vecpi=rep(1/9,9)} else {vecpi=options$vecpi}
  if (!("tau" %in% parnames)){       tau=1.0} else {tau=options$tau}
  if (!("niter" %in% parnames)){     niter=20000} else {niter=options$niter}
  if (!("t0" %in% parnames)){        t0=200} else {t0=options$t0}
  if (!("xi" %in% parnames)){        xi=2/3} else {xi=options$xi}
  if (!("stepsize" %in% parnames)){  stepsize=1.0} else {stepsize=options$stepsize}
  if (!("trange" %in% parnames)){    trange=c(-Inf,Inf)} else {trange=options$trange}
  
  if ((length(nv)>1)||(!is.numeric(nv))||(nv<1)||(is.infinite(nv))||(is.na(nv))){
    stop(cat("* ",fname," : 'nv' should be a positive integer.",sep=""))
  }
  nv = as.integer(nv)
  if (is.vector(domain)){
    if ((length(domain)!=2)||(any(!is.numeric(domain)))||(any(is.na(domain)))){
      stop(cat("* ",fname," : 'domain' should be a vector of length 2.",sep=""))
    }
    domain = matrix(sort(domain),nrow=nv,ncol=length(domain),byrow=TRUE)
  } else if (is.matrix(domain)){
    if ((dim(domain)[1]!=nv)||(dim(domain)[2]!=2)||(any(is.na(domain)))||(any(!is.numeric(domain)))){
      stop(cat("* ",fname," : 'domain' should be a matrix of size (nv-by-2).",sep=""))
    }
    for (i in 1:nv){
      domain[i,] = sort(domain[i,])
    }
  } else {
    stop(cat("* ",fname," : 'domain' should be either a vector or a matrix.",sep=""))
  }
  if (any(is.infinite(domain))){
    #message(cat("* ",fname," : 'domain' would be better without infinite values.",sep=""))
    #message("*      : automatically replacing Inf's with suitably large numbers.")
    domain = adjust_inf(domain)
  }
  if ((!is.vector(partition))||(length(partition)<2)||(any(!is.numeric(partition)))){
    stop(cat("* ",fname," : 'partition' should be a vector whose length is greater than 1.",sep=""))
  }
  partition = as.vector(sort(partition))
  if ((any(is.null(vecpi)))||(length(vecpi)!=(length(partition)-1))||(any(vecpi<=0))||(abs(sum(vecpi)-1)>1e-10)){
    stop(cat("* ",fname," : desired sampling distribution 'vecpi' is invalid.",sep=""))
  }
  vecpi = as.vector(sort(vecpi))
  tau = as.double(tau)
  if ((!is.numeric(niter))||(niter<=1)||(is.infinite(niter))||(is.na(niter))||(length(niter)>1)){
    stop(cat("* ",fname," : 'niter' should be a positive integer as an iteration number.",sep=""))
  }
  t0 = as.double(t0)
  if ((!is.numeric(xi))||(xi<=0.5)||(xi>1)||(length(xi)>1)||(is.null(xi))){
    stop(cat("* ",fname," : 'xi' should be in (0.5,1].",sep=""))
  }
  xi = as.double(xi)
  if ((any(!is.numeric(stepsize)))||(any(stepsize<=0))||(any(is.infinite(stepsize)))){
    stop(cat("* ",fname," : 'stepsize' is a standard deviation term for normal proposal density.",sep=""))
  }
  stepsize = as.double(stepsize)
  if ((length(stepsize)==1)&&(length(stepsize)<nv)){
    stepsize = rep(stepsize, nv)
  } 
  if (length(stepsize)!=nv){
    stop(cat("* ",fname," : 'stepsize' should have a length of 'nv'."))
  }
  m = length(vecpi)
  if (is.vector(trange)){
    if ((length(trange)!=2)||(any(!is.numeric(trange)))||(any(is.na(trange)))){
      stop(cat("* ",fname," : 'trange' should be a vector of length 2.",sep=""))
    }
    trange = matrix(sort(trange),nrow=m,ncol=length(trange),byrow=TRUE)
  } else if (is.matrix(trange)){
    if ((dim(trange)[1]!=m)||(dim(trange)[2]!=2)||(any(is.na(trange)))||(any(!is.numeric(trange)))){
      stop(cat("* ",fname," : 'trange' should be a matrix of size (m-by-2).",sep=""))
    }
    for (i in 1:m){
      trange[i,] = sort(trange[i,])
    }
  } else {
    stop(cat("* ",fname," : 'trange' should be either a vector or a matrix.",sep=""))
  }
  
  ## let's try to return the results
  #  'domain','partition','vecpi','tau','niter','t0','xi','stepsize','trange'
  pruned = list()
  pruned$domain = domain
  pruned$partition = partition
  pruned$vecpi = vecpi
  pruned$tau = tau
  pruned$niter = niter
  pruned$t0 = t0
  pruned$xi = xi
  pruned$stepsize = stepsize
  pruned$trange = trange
  pruned$m = m
  return(pruned)
}



#  4. dist_dinvgamma ------------------------------------------------------
#' @keywords internal
#' @noRd
dist_dinvgamma <- function(x, alpha, beta){
  return(((beta^alpha)/gamma(alpha))*(x^(-alpha-1))*exp(-beta/x))
}


# 5. vec_init -------------------------------------------------------------
#' @keywords internal
#' @noRd
vec_init <- function(nv,domain){
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
  return(init)
}

#  6. compute_pinv --------------------------------------------------------
#' @keywords internal
#' @noRd
compute_pinv <- function(mat){
  return(solve(mat,diag(nrow(mat))))
}
