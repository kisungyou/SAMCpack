SAMCpack: Stochastic Approximation Monte Carlo (SAMC) Sampler and Methods
=========================================================================

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/SAMCpack?color=green)](https://CRAN.R-project.org/package=SAMCpack) [![Travis build status](https://travis-ci.org/kisungyou/SAMCpack.svg?branch=master)](https://travis-ci.org/kisungyou/SAMCpack) [![](https://cranlogs.r-pkg.org/badges/SAMCpack)](https://cran.r-project.org/package=SAMCpack)

Stochastic Approximation Monte Carlo (SAMC) is one of the celebrated Markov chain Monte Carlo (MCMC) algorithms. It is known to be capable of sampling from multimodal or doubly intractable distributions. We provide generic SAMC samplers for continuous distributions. User-specified densities in R and C++ are both supported. We also provide functions for specific problems that exploit SAMC computation. See [Liang et al (2010)](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470669723) for complete introduction to the method.

Installation
------------

You can install the released version of SAMCpack from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SAMCpack")
```

or the development version from github:

``` r
## install.packages("devtools")
devtools::install_github("kisungyou/SAMCpack")
```
