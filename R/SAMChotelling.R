#' Two-sample Hotelling's T-squared test with Resampling-based SAMC
#' 
#' SAMC is used to estimate (possibly) extremal \eqn{p}-value (e.g. smaller than \eqn{1e-6}) 
#' for two-sample Hotelling's T-squared test with multivariate samples as an alternative to naive permutation test. 
#' 
#' @param x matrix of \emph{1st} sample.
#' @param y matrix of \emph{2nd} sample.
#' @param niter number of SAMC iterations.
#' @param m number of partition refinement.
#' @param t0 gain factor control at which sequence starts to decrease.
#' @param sample.ratio ratio in \eqn{(0,1)} at which resampling is performed.
#' 
#' @return a named list containing \describe{
#' \item{p.val}{estimated \eqn{p}-value.}
#' \item{statistic}{observed statistic.}
#' \item{acceptance}{ratio of accepted runs from \code{niter} iterations.}
#' \item{frequency}{visiting frequency of energy partition.}
#' \item{max.error}{maximal discrepancy. Rule of thumb to determine convergence when the value is \eqn{< 0.2}.}
#' }
#' 
#' @examples
#' \dontrun{
#' ## generate samples from normal distribution for small p-value
#' x = matrix(rnorm(100*3),nrow=100)          # standard normal
#' y = matrix(rnorm(100*3,mean=0.3),nrow=100) # mean are (0.5,0.5,0.5)
#' 
#' ## run SAMCttest with different partition scales
#' #    default iteration of 1e+5 subsamples are used.
#' res1 = SAMChotelling(x,y,m=100)
#' res2 = SAMChotelling(x,y,m=300)
#' res3 = SAMChotelling(x,y,m=500)

#' ## we also use 1e+4 permutations 
#' require(Hotelling) 
#' hot.asym = Hotelling::hotelling.test(x,y)
#' hot.perm = Hotelling::hotelling.test(x,y,perm=TRUE,B=1e+5)
#' 
#' ## compare with true p-value from 't.test' and permutation test
#' printer <- paste(
#'   "=========================================\n",
#'   "      Compare Different p-values\n",
#'   "=========================================\n",
#'   "1. Hotelling     : ",hot.asym$pval,
#'   "\n2. SAMCttest",
#'   "\n      2-1. m=100 : ",res1$p.val,
#'   "\n      2-2. m=300 : ",res2$p.val,
#'   "\n      2-3. m=500 : ",res3$p.val,
#'   "\n3. Permutation   : ",hot.perm$pval, sep="")
#' writeLines(printer)
#' 
#' ## visualize visiting frequency
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' barplot(res1$frequency/1e+5, main=paste("max.error=",res1$max.error,sep=""))
#' barplot(res2$frequency/1e+5, main=paste("max.error=",res2$max.error,sep=""))
#' barplot(res3$frequency/1e+5, main=paste("max.error=",res3$max.error,sep=""))
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{yu_efficient_2011}{SAMCpack}
#' 
#' @export
SAMChotelling = function(x, y, niter=1e+5, m=300, t0=1000, sample.ratio=0.05){
  # niter = B  : number of iterations 
  # m     = m  : length of partition
  # t0    = t0 : iteration at which the updates begin to decay
  ##-------------------------------------------------------------------
  # PREPROCESSING
  x.old = x
  y.old = y
  if (any(is.infinite(x.old))||any(is.na(x.old))){
    stop("* SAMChotelling : no Inf or NaN values are allowed for 1st sample 'x'.")
  }
  if (any(is.infinite(y.old))||any(is.na(y.old))){
    stop("* SAMChotelling : no Inf or NaN values are allowed for 2nd sample 'y'.")
  }
  if (ncol(x.old)!=ncol(y.old)){
    stop("* SAMChotelling : dimension of two samples must match.")
  }
  m = as.integer(m)
  t0 = as.integer(t0)
  if (length(sample.ratio)>1||(sample.ratio<=0)||(sample.ratio>=1)){
    stop("* SAMChotelling : 'sample ratio' should be a scalar in (0,1).")
  }

  ##-------------------------------------------------------------------
  # INITIALIZATION
  ts.init      = hotelling.statistic(x.old, y.old)  # initial test statistic
  count.accept = 0         # number of acceptance
  frequency    = rep(0, m) # vector to record visiting frequency
  partition    = seq(0, ts.init, length.out=m) # energy partition
  
  thetas = rep(0,m)    # estimate of a theta vector 
  vecpi  = rep(1,m)/m  # desired sampling distribution
  
  L = ceiling(min(nrow(x.old), nrow(y.old))*sample.ratio) # number of resampling batch
  
  ts.old = ts.init                     # initial value for test statistic
  id.old = sum(ts.old >= partition)    # initial index
  
  ##-------------------------------------------------------------------
  # MAIN ITERATION
  for (b in 1:niter){
    # Step 1. update proposal
    resampler = hotelling.resample(x.old, y.old, L)
    x.new = resampler$xsample
    y.new = resampler$ysample
    
    # Step 2. compute test statistic, partition, and acceptance probability
    ts.new = hotelling.statistic(x.new, y.new)
    id.new = sum(ts.new >= partition)
    r = min(exp(thetas[id.old] - thetas[id.new]), 1.0)
    
    # Step 3. dichotomous branching for accept/reject
    if (runif(1)<=r){
      count.accept = count.accept + 1
      x.old = x.new
      y.old = y.new
      ts.old = ts.new
      id.old = id.new
    }
    id.update = rep(0, m)
    id.update[id.old] = 1
    thetas = thetas + (t0 / max(b, t0)) * (id.update-vecpi)
    
    # Step 4. update frequency
    frequency[id.old] = frequency[id.old] + 1
  }
  
  ##-------------------------------------------------------------------
  # POSTPROCESSING
  #   1. compute p-value
  p.val = exp(thetas[m]) * vecpi[m] / sum(exp(thetas) * vecpi)
  
  #   2. maximum error : max < 0.2 : well converged
  m0 = sum(frequency==0)
  frequency.norm = frequency/niter
  max.error = max(abs(((frequency.norm-1/(m-m0))/(1/(m-m0))*(frequency>0))))
  
  
  # return values
  output = list()
  output$p.val = p.val
  output$statistic = ts.init
  output$acceptance = count.accept/niter
  output$frequency  = frequency
  output$max.error  = max.error
  return(output)
}


# extra functions for SAMChotelling ---------------------------------------
#' @keywords internal
#' @noRd
hotelling.resample <- function(x, y, L){
  idx = sample(1:nrow(x), L)
  idy = sample(1:nrow(y), L)
  
  xsample = rbind(x[-idx,], y[idy,])
  ysample = rbind(x[idx,], y[-idy,])
  
  output = list()
  output$xsample = xsample
  output$ysample = ysample
  return(output)
}
#' @keywords internal
#' @noRd
hotelling.statistic <- function(x, y){
  nx = nrow(x)
  ny = nrow(y)
  
  xbar = colMeans(x)
  ybar = colMeans(y)
  dbar = as.vector(xbar-ybar)
  
  Sigma = ((nx-1)*stats::cov(x) + (ny-1)*stats::cov(y))/(nx+ny-2)

  output =  (sum(as.vector(solve(Sigma, dbar))*dbar)*((nx*ny)/(nx+ny)))
  return(output)
}
# x = matrix(rnorm(100*3),nrow=100)
# y = matrix(rnorm(100*3,mean=0.5),nrow=100)
# myiter=5000;
# c(SAMChotelling(x,y,niter=myiter, sample.ratio = 0.075)$p.val,
#   Hotelling::hotelling.test(x,y,perm=TRUE,B=myiter)$pval,
#   Hotelling::hotelling.test(x,y)$pval)