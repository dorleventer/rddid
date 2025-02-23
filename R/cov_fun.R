cov_fun = function(ids1, ids2, res1, res2) {
  sigma = matrix(NA, length(ids1), length(ids2))
  for(i in 1:length(ids1)) {
    for(j in 1:length(ids2)) {
      if(ids1[i,1] == ids2[j,1]) {
        sigma[i,j] <- res1[i,1] * res2[j,1]
      } else {
        sigma[i,j] <- 0
      }
    }
  }
  return(sigma)
}
