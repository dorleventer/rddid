covariances_bc = function(x1, x2, y1, y2, h, b, ids1, ids2) {

  mat1 = mat_fun(x1,y1,h,b, ids = ids1)
  mat2 = mat_fun(x2,y2,h,b, ids = ids2)
  A1_l = mat1$invG_p_l %*% t(mat1$Q_q_l)
  A2_l = mat2$invG_p_l %*% t(mat2$Q_q_l)
  A1_r = mat1$invG_p_r %*% t(mat1$Q_q_r)
  A2_r = mat2$invG_p_r %*% t(mat2$Q_q_r)

  sigma_1_2_l_l = cov_fun(mat1$ids_l, mat2$ids_l, mat1$res_b_l, mat2$res_b_l)
  sigma_1_2_r_r = cov_fun(mat1$ids_r, mat2$ids_r, mat1$res_b_r, mat2$res_b_r)
  sigma_1_2_r_l = cov_fun(mat1$ids_r, mat2$ids_l, mat1$res_b_r, mat2$res_b_l)
  sigma_1_2_l_r = cov_fun(mat1$ids_l, mat2$ids_r, mat1$res_b_l, mat2$res_b_r)

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
