# COMMON FUNCTIONS --------------------------------------------------------
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