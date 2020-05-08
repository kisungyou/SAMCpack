#' Stochastic Approximation Monte Carlo (SAMC) Sampler and Methods
#'
#' Stochastic Approximation Monte Carlo (SAMC) is one of the celebrated Markov chain Monte Carlo (MCMC) algorithms. It is known to be capable of sampling from 
#' multimodal or doubly intractable distributions. We provide generic SAMC samplers for continuous distributions. User-specified densities in R and C++ are both supported.
#' We also provide functions for specific problems that exploit SAMC computation. See Liang et al (2010) <doi:10.1002/9780470669723> for complete introduction to the method.
#' 
#' @author Kisung You
#' @docType package
#' @name SAMCpack-package
#' @import RcppXPtrUtils
#' @import Rdpack
#' @importFrom stats runif var cov rnorm t.test
#' @importFrom utils packageVersion tail
#' @importFrom Rcpp evalCpp
#' @useDynLib SAMCpack, .registration=TRUE
NULL


