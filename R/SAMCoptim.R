#' Stochastic Approximation Annealing
#' 
#' 
#' 
#' @examples
#' \donttest{
#' ##### Multidimensional Objective Function 
#' ## Step 1 : Define negative log-density function as an R function
#' func_r = function(x){
#' x1 = x[1]; x2 = x[2];
#' val1 = (-(x1*sin(20*x2)+x2*sin(20*x1))^2)*cosh(sin(10*x1)*x1);
#' val2 = (-(x1*cos(10*x2)-x2*sin(10*x1))^2)*cosh(cos(20*x2)*x2);
#' return(val1+val2);
#' }
#' 
#' ## Step 2 : Prepare a setting option
#' myoption = list()
#' myoption$partition = c(-Inf,seq(from=-10,to=0,length.out=41))
#' myoption$tau       = 1.0
#' myoption$domain    = c(-1.1,1.1)
#' myoption$vecpi     = as.vector(rep(1/41,41))
#' myoption$niter     = 20000
#' myoption$stepsize  = c(0.25, 10)
#' 
#' ## Step 3 : Run The Code
#' res = SAMCoptim(func_r,lower=c(-1.1,-1.1), upper=c(1.1,1.1), options=myoption)
#' 
#' ## Step 4 : Visualize
#' par(mfrow=c(1,2))
#' #  4-1. plot all samples
#' select = seq(from=101,to=myoption$niter,by=100) # 100 burn-in, 1/100 thinning 
#' sam.x  = res$samples[select,1] # 1st coordinate
#' sam.y  = res$samples[select,2] # 2nd coordinate
#' plot(sam.x,sam.y,xlab='x',ylab='y',main='samples',xlim=c(-1.2,1.2),ylim=c(-1.2,1.2))
#' 
#' #  4-2. plot samples with minimal values
#' nval   = 3
#' idmin  = which(res$fnval<=sort(res$fnval))
#' opt.x  = res$samples[idmin,1] # 1st coordinate
#' opt.y  = res$samples[idmin,2] # 2nd coordinate
#' plot(opt.x,opt.y,xlab='x',ylab='y',main='optimal samples',xlim=c(-1.2,1.2),ylim=c(-1.2,1.2))
#' 
#' #  4-3. evoluation
#' plot(1:length(select),res$fnval[select],"b",xlab="iteration",ylab="function value")
#' }
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
  
  ##------------------------------------------------)-------------------
  # RUN
  output = exec_SAMCoptim(Ux,nv,partition,domain,niter,vecpi,t0,xi,stepsize,init,
                          temp.tau.h, temp.tau.star, temp.t0, trunc.init, trunc.r)
  # EXTRA STEP TO EVALUATE
  fnval        = apply(output$samples, 1, Ux)
  output$fnval = fnval
  output$bestmem = output$samples[which.min(fnval),]
  output$bestval = min(fnval)
  ##-------------------------------------------------------------------
  # END
  return(output)   
}