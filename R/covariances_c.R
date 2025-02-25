#' The four covariances between 2 time periods for conventional estimation
#'
#' @param x1 running variable in period 1
#' @param x2 running variable in period 2
#' @param y1 outcome variable in period 1
#' @param y2 outcome variable in period 2
#' @param h bandwidth
#' @param b bias correction bandwidth
#' @param ids1 id vector in period 1
#' @param ids2 id vector in period 2
#'
#' @return four covariances
covariances_c = function(x1, x2, y1, y2, h, b, ids1, ids2) {
  mat1 = mat_fun(x1,y1,h,b, ids = ids1)
  mat2 = mat_fun(x2,y2,h,b, ids = ids2)
  A1_l = mat1$invG_p_l %*% t(mat1$R_p_l*mat1$W_h_l)
  A2_l = mat2$invG_p_l %*% t(mat2$R_p_l*mat2$W_h_l)
  A1_r = mat1$invG_p_r %*% t(mat1$R_p_r*mat1$W_h_r)
  A2_r = mat2$invG_p_r %*% t(mat2$R_p_r*mat2$W_h_r)

  sigma_1_2_l_l = cov_fun(mat1$ids_l, mat2$ids_l, mat1$res_h_l, mat2$res_h_l)
  sigma_1_2_r_r = cov_fun(mat1$ids_r, mat2$ids_r, mat1$res_h_r, mat2$res_h_r)
  sigma_1_2_r_l = cov_fun(mat1$ids_r, mat2$ids_l, mat1$res_h_r, mat2$res_h_l)
  sigma_1_2_l_r = cov_fun(mat1$ids_l, mat2$ids_r, mat1$res_h_l, mat2$res_h_r)

  c_1_2_l_l = A1_l %*% sigma_1_2_l_l %*% t(A2_l)
  c_1_2_r_r = A1_r %*% sigma_1_2_r_r %*% t(A2_r)
  c_1_2_l_r = A1_l %*% sigma_1_2_l_r %*% t(A2_r)
  c_1_2_r_l = A1_r %*% sigma_1_2_r_l %*% t(A2_l)

  return(
    list(
      "c_1_2_l_l" = c_1_2_l_l,
      "c_1_2_r_r" = c_1_2_r_r,
      "c_1_2_l_r" = c_1_2_l_r,
      "c_1_2_r_l" = c_1_2_r_l
    )
  )
}
