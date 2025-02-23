#' RD estimates for single time period
#'
#' @param x running variable
#' @param y outcome variable
#' @param h bandwidth
#' @param b bias correction bandwidth
#' @param v derivative
#' @param p polynomial
#' @param q bias correction polynomial
#' @param c cutoff
#' @param kernel type of kernel
#'
#' @return Outcome discontinuity estimates and standard errors
single_period_est = function(x,y,h,b,v=0,p=1,q=2,c=0,kernel="") {

  h_l = h_r = h
  b_l = b_r = b

  x = as.matrix(x)
  y = as.matrix(y)

  order_x = order(x)
  x = x[order_x,,drop=FALSE]
  y = y[order_x,,drop=FALSE]

  ind_l = x<c
  ind_r = x>=c
  X_l = x[ind_l,,drop=FALSE]
  X_r = x[ind_r,,drop=FALSE]
  Y_l = y[ind_l,,drop=FALSE]
  Y_r = y[ind_r,,drop=FALSE]

  w_h_l <- rdrobust_kweight(X_l,c,h_l,kernel)
  w_h_r <- rdrobust_kweight(X_r,c,h_r,kernel)
  w_b_l <- rdrobust_kweight(X_l,c,b_l,kernel)
  w_b_r <- rdrobust_kweight(X_r,c,b_r,kernel)

  ind_h_l <- w_h_l> 0
  ind_h_r <- w_h_r> 0
  ind_b_l <- w_b_l> 0
  ind_b_r <- w_b_r> 0

  ind_l = ind_b_l
  ind_r = ind_b_r

  eY_l  = Y_l[ind_l,,drop=FALSE]
  eY_r  = Y_r[ind_r,,drop=FALSE]
  eX_l  = X_l[ind_l,,drop=FALSE]
  eX_r  = X_r[ind_r,,drop=FALSE]
  W_h_l = w_h_l[ind_l]
  W_h_r = w_h_r[ind_r]
  W_b_l = w_b_l[ind_l]
  W_b_r = w_b_r[ind_r]

  u_l <- (eX_l-c)/h_l
  u_r <-(eX_r-c)/h_r
  R_q_l = matrix(NA,length(eY_l),(q+1))
  R_q_r = matrix(NA,length(eY_r),(q+1))
  for (j in 1:(q+1))  {
    R_q_l[,j] = (eX_l-c)^(j-1)
    R_q_r[,j] = (eX_r-c)^(j-1)
  }
  R_p_l = R_q_l[,1:(p+1)]
  R_p_r = R_q_r[,1:(p+1)]

  theta_l = crossprod(R_p_l*W_h_l,u_l^(p+1))
  theta_r = crossprod(R_p_r*W_h_r,u_r^(p+1))
  invG_q_l  = qrXXinv((sqrt(W_b_l)*R_q_l))
  invG_q_r  = qrXXinv((sqrt(W_b_r)*R_q_r))
  invG_p_l  = qrXXinv((sqrt(W_h_l)*R_p_l))
  invG_p_r  = qrXXinv((sqrt(W_h_r)*R_p_r))
  e_p1 = matrix(0,(q+1),1); e_p1[p+2]=1
  e_v  = matrix(0,(p+1),1); e_v[v+1]=1

  Q_q_l = t(t(R_p_l*W_h_l) - h_l^(p+1)*(theta_l%*%t(e_p1))%*%t(t(invG_q_l%*%t(R_q_l))*W_b_l))
  Q_q_r = t(t(R_p_r*W_h_r) - h_r^(p+1)*(theta_r%*%t(e_p1))%*%t(t(invG_q_r%*%t(R_q_r))*W_b_r))
  D_l = eY_l
  D_r = eY_r

  beta_p_l  = invG_p_l%*%crossprod(R_p_l*W_h_l,D_l)
  beta_p_r  = invG_p_r%*%crossprod(R_p_r*W_h_r,D_r)
  beta_q_l  = invG_q_l%*%crossprod(R_q_l*W_b_l,D_l)
  beta_q_r  = invG_q_r%*%crossprod(R_q_r*W_b_r,D_r)
  beta_bc_l = invG_p_l%*%crossprod(Q_q_l,D_l)
  beta_bc_r = invG_p_r%*%crossprod(Q_q_r,D_r)
  beta_p    = beta_p_r  - beta_p_l
  beta_bc   = beta_bc_r - beta_bc_l

  tau_cl = beta_p[(v+1),1]
  tau_bc = beta_bc[(v+1),1]

  predicts_p_l=R_p_l%*%beta_p_l
  predicts_p_r=R_p_r%*%beta_p_r
  predicts_q_l=R_q_l%*%beta_q_l
  predicts_q_r=R_q_r%*%beta_q_r

  res_h_l = rdrobust_res(eY_l, predicts_p_l, p+1)
  res_h_r = rdrobust_res(eY_r, predicts_p_r, p+1)

  res_b_l = rdrobust_res(eY_l, predicts_q_l, q+1)
  res_b_r = rdrobust_res(eY_r, predicts_q_r, q+1)

  V_Y_cl_l = invG_p_l%*%crossprod(c(res_h_l)*R_p_l*W_h_l)%*%invG_p_l
  V_Y_cl_r = invG_p_r%*%crossprod(c(res_h_r)*R_p_r*W_h_r)%*%invG_p_r
  V_Y_rb_l = invG_p_l%*%crossprod(c(res_b_l)*Q_q_l)%*%invG_p_l
  V_Y_rb_r = invG_p_r%*%crossprod(c(res_b_r)*Q_q_r)%*%invG_p_r
  V_tau_cl = (V_Y_cl_l+V_Y_cl_r)[v+1,v+1]
  V_tau_rb = (V_Y_rb_l+V_Y_rb_r)[v+1,v+1]
  se_tau_cl = sqrt(V_tau_cl)
  se_tau_rb = sqrt(V_tau_rb)

  tau = c(tau_cl, tau_bc)
  se  = c(se_tau_cl,se_tau_rb)

  tibble(
    type = c("C","BC"),
    est = tau,
    se = se
  )
}
