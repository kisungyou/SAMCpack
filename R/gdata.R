#' Dataset for SAMCrsa
#' 
#' Simulated data list with response variable distributed with 'powered exponential' covariance matrix, 
#' with \eqn{\phi=25, \tau=1, \kappa=1, \sigma=1}.
#'
#' @usage 
#' data(gdata)
#' 
#' @format 
#' A sample dataset with 1000 observations. It's save as an R's native \code{list}, and see the table for description of variables in the \code{gdata} list,
#' \tabular{lll}{
#' VARIABLE \tab SIZE \tab DESCRIPTION \cr
#' \code{coords} \tab \eqn{(1000\times 2)} matrix \tab 2D location for each observation \cr
#' \code{y}      \tab length-\eqn{1000} vector \tab response value \cr
#' \code{X}      \tab length-\eqn{1000} vcetor \tab covariate
#' }
#' This is the default example data for \code{\link{SAMCrsa}}.
#' 
#' @details 
#' Below is the code used to generate \emph{gdata}:
#' \preformatted{
#' require("geoR")
#' require("RandomFields")
#' DataNum=1000
#' gData=grf(DataNum,grid="irreg",DataNum,DataNum,xlims=c(0,100),ylims=c(0,100),nsim=1,mean=0,
#'           cov.mode="powered.exponential",cov.par=c(1,25),nugget=1,kappa=1)
#' x=rnorm(DataNum)
#' gdata = list(y=gData$data+.5+x, X=x, coords=gData$coords)
#' }
"gdata"

