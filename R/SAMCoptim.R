#' Stochastic Approximation Annealing
#' 
#' 
#' @export
SAMCoptim <- function(fn, lower, upper, options=list(), ...){
  ## -------------------------------------------------------------------
  ## Checking the arguments
  if (!is.function(fn)){  stop("* SAMCoptim : 'fn' needs a function.")}
  if (!is.list(options)){ stop("* SAMCoptim : 'options' must be a list.")}
  if (!is.vector(lower)){ stop("* SAMCoptim : 'lower' must be a vector.")}
  if (!is.vector(upper)){ stop("* SAMCoptim : 'upper' must be a vector.")}
  if (length(lower)!=length(upper)){
    stop("* SAMCoptim : 'lower' and 'upper' must be of same length.")
  }
  ## Now, we need a bit different approach on parameter checking
  options$nv     = length(lower)
  options$domain = cbind(lower,upper) # generate conventional domain matrix on my own.
  
  ## separate out the parameters
  pruned    = optionlist.optim(options,options$nv,"SAMCoptim")
  domain    = pruned$domain
  partition = pruned$partition
  vecpi     = pruned$vecpi
  niter     = pruned$niter
  t0        = pruned$t0
  xi        = pruned$xi
  stepsize  = pruned$stepsize
  trange    = pruned$trange
  m         = pruned$m  # added for fun
  temp.tau.h    = pruned$temp.tau.h # five extra arguments
  temp.tau.star = pruned$temp.tau.star
  temp.t0       = pruned$temp.t0
  trunc.init    = pruned$trunc.init
  trunc.r       = pruned$trunc.r
  
  ##------------------------------------------------)-------------------
  # INITIALIZATION
  #   Domain
  nv = options$nv
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
  #   Wrapped Function
  Ux <- function(x){
    return(fn(x, ...))
  }
  return(Ux)
}