rdrobust_res = function(y, m, d) {
  n = length(y)
  res = matrix(NA,n,1)
  w = sqrt(n/(n-d))
  res[,1] = w*(y-m[,1])
  return(res)
}
