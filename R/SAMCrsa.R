#' A Resampling-based Stochastic Approximation Method for Analysis of Large Geostatitical data
#' 
#' Performs parameter estimation using a resampling-based Stochastic Approximation (RSA) method. 
#' It is a stochatic approximation method. At every iteration, only a subset of the data is drawn and used to update the estimation of the parameters. 
#' The data are assumed to have a powered exponential correlation structure.
#' 
#' @param coords an \eqn{(n\times 2)} matrix. 2D location coordinates.
#' @param y      a length-\eqn{n} vector of response value.
#' @param X      an \eqn{(n\times k)} matrix of extra covariates.
#' @param nsubset the size of the subset drawn from the data. It is recommended to be set to 300 or higher.
#' @param stepscale gain factor control. It specifies the number of iterations when the gain factor begins to shrink. For example, one can be set it equal to 2 times the burn-in steps.
#' @param niter the total number of iterations for stochastic approximation. In practice, it is recommended to be set to 2500 or higher.
#' @param warm the number of burn-in iterations
#' 
#' @return a named list containing \describe{
#' \item{beta}{the coefficient estimates of the mean effect. It is a vector of length equal to the number of coefficients plus 1.}
#' \item{phi}{the shape estimate in the powered exponential correlation matrix.}
#' \item{sigmasq}{the estimate of error variance.}
#' \item{tausq}{the estimate of nugget variance.}
#' }
#' 
#' @examples
#' ##### load example data pre-loaded
#' data(gdata)
#' 
#' ##### run RSA
#' output = SAMCrsa(gdata$coords, gdata$y, gdata$X, nsubset=50, stepscale=40, niter=100, warm=20)
#' 
#' @references
#' \insertRef{SAMCrsa}{SAMCpack}
#' 
#' @author Yichen Cheng, Faming Liang, Kisung You
#' @export
SAMCrsa <- function(coords,y,X=NULL,nsubset=max(ceiling(length(y)/5),10),
                    stepscale=200,niter=2500,warm=100){
  if ((!is.matrix(coords))||(ncol(coords)!=2)){
    stop("* SAMCrsa : 'coords' should be a matrix of 2 columns.")
  }
  if (is.vector(y)){
    y = matrix(y)
  }
  N = nrow(coords)
  if (N!=nrow(y)){
    stop("* SAMCrsa : nrow(coords) should be equal to length(y).")
  }
  if ((is.null(X))&&(length(X)==0)){
    data = cbind(coords,y)
  } else {
    if (is.vector(X)){
      X = matrix(X)
    }
    if ((!is.matrix(X))||(nrow(X)!=N)){
      stop("* SAMCrsa : if given covariates 'X', it must have 'nrow(coords)' number of rows.")
    }
    data = cbind(coords,y,X)
  }
  
  if ((nsubset < 2)||(nsubset >= nrow(data))){
    stop("* SAMCrsa : 'nsubset' is invalid.")
  }
  if (stepscale <= 0){
    stop("* SAMCrsa : 'stepscale' should be a positive number.")
  }
  if (niter < 5){
    stop("* SAMCrsa : 'niter' should be a decently large positive number.")
  }
  if (warm < 5){
    stop("* SAMCrsa : 'warm' controls the burn-in's.")
  }
  
  dataCol = ncol(data);
  dataNum = nrow(data);
  beta = rep(0,dataCol-2);
  phi = 0;
  sigmasq = 0;
  tausq = 0;
  # RetVec2 = .C("RSAc",
  # as.numeric(data),
  # as.integer(dataCol),
  # as.integer(dataNum),
  # as.integer(nsubset),
  # as.integer(stepscale),
  # as.integer(niter),
  # as.integer(warm),
  # as.numeric(beta),
  # as.numeric(phi),
  # as.numeric(sigmasq),
  # as.numeric(tausq)
  # )
  # 
  # beta = RetVec2[[8]];
  # phi = RetVec2[[9]];
  # sigmasq = RetVec2[[10]];
  # tausq = RetVec2[[11]];
  
  RetVec2 = RSAarma(as.numeric(data),as.integer(dataCol),as.integer(dataNum),as.integer(nsubset),
                    as.integer(stepscale),as.integer(niter),as.integer(warm))
  # beta    = RetVec2$pbeta
  # phi     = RetVec2$pPhi
  # sigmasq = RetVec2$pSigmasq
  # tausq   = RetVec2$pTausq
  beta    = RetVec2[[8]]
  phi     = RetVec2[[9]]
  sigmasq = RetVec2[[10]]
  tausq   = RetVec2[[11]]
  
  Z = NULL
  z = list(beta = beta,phi=phi,sigmasq=sigmasq,tausq=tausq)
  #return(c(beta0,beta1,phi,sigmasq,tausq));
  return(z)
}