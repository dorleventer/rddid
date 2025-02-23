qrXXinv = function(x) {
  chol2inv(chol(crossprod(x)))
}
