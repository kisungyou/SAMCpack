#' A Resampling-based Stochastic Approximation Method for Analysis of Large Geostatitical data
#' 
#' Performs parameter estimation using a resampling-based Stochastic Approximation (RSA) method. 
#' It is a stochatic approximation method. At every iteration, only a subset of the data is drawn and used to update the estimation of the parameters. 
#' The data are assumed to have a powered exponential correlation structure.
#' 
#' @param data an \eqn{(n\times p)} matrix. The first two columns of the data are 
#' 2D location coordinates (\eqn{x_1, x_2}) for each observation. The third column gives the 
#' response value (\eqn{y}). Extra covariates can be included from the 4th column and beyond.
#' @param nsubset the size of the subset drawn from the data. It is recommended to be set to 300 or higher.
#' @param stepscale gain factor control. It specifies the number of iterations when the gain factor begins to shrink. For example, one can be set it equal to 2 times the burn-in steps.
#' @param niter the total number of iterations for stochastic approximation. In practice, it is recommended to be set to 2500 or higher.
#' @param warm the number of burn-in iterations
#' 
#' @return a named list containing \describe{
#' \item{beta}{the coefficient estimates of the mean effect. It is a vector of lenght equal to the number of coefficients plus 1.}
#' \item{phi}{the shape estimate in the powered exponential correlation matrix.}
#' \item{sigmasq}{the estimate of error variance.}
#' \item{tausq}{the estimate of nugget variance.}
#' }
#' 
#' @examples 
#' ## Load Data
#' data(gdata)
#' 
#' ## Run RSA
#' SAMCrsa(gdata, nsubset=50, stepscale=40, niter=100, warm=20)
#' 
#' @references
#' \insertRef{SAMCrsa}{SAMCpack}
#' 
#' @author Yichen Cheng, Faming Liang, Kisung You
#' @export
SAMCrsa <- function(data,nsubset=max(ceiling(nrow(data)/5),10),
                    stepscale=200,niter=2500,warm=100){
  if ((!is.matrix(data))||(ncol(data)<3)){
    stop("* SAMCrsa : 'data' must be a matrix of at least 3 columns.")
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
  RetVec2 = .C("RSAc",
               as.numeric(data),
               as.integer(dataCol),
               as.integer(dataNum),
               as.integer(nsubset),
               as.integer(stepscale),
               as.integer(niter),
               as.integer(warm),             
               as.numeric(beta),
               as.numeric(phi),
               as.numeric(sigmasq),
               as.numeric(tausq)
  )

  beta = RetVec2[[8]];
  phi = RetVec2[[9]];
  sigmasq = RetVec2[[10]];
  tausq = RetVec2[[11]];
  Z = NULL
  z = list(beta = beta,phi=phi,sigmasq=sigmasq,tausq=tausq)
  #return(c(beta0,beta1,phi,sigmasq,tausq));
  return(z)
}