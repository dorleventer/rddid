rdrobust_vce = function(d, s, RX, res, C) {
  k = ncol(as.matrix(RX))
  M = matrix(0,k,k)
  n  = length(C)
  M  = crossprod(c(res)*RX)
  return(M)
}
