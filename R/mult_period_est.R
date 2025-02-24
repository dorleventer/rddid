#' RDDID Estimation
#'
#' This function estimates the ATU or ATT using a RDDID framework, assuming potential outcome discontinuities are constant across time.
#'
#' @param df dataframe
#' @param t_star time period with RD treatment
#' @param t_vec time period for bias estimation
#' @param w_vec weights for t_vec time periods
#' @param h bandwidth
#' @param b bias correction bandwidth
#' @param tname the name of the time period variable in df
#' @param xname the name of the running variable in df
#' @param yname the name of the outcome variable in df
#' @param idname the name of the id variable in df
#'
#' @return Conventional (C) and Bias Corrected (BC) treatment effect estimates and standard errors (SE). SEs are provided for sampling under repeated cross section (CS), panel with constant running variable (PC) and varying running variable (PV).
#' @export
mult_period_est = function(df, t_star, t_vec, w_vec, h, b, tname = "time", xname = "R", yname = "Y", idname = "id") {
  # prep the data
  time = df[,tname]
  x = df[,xname]
  y = df[,yname]
  ids = df[,idname]
  df = data.frame(time, R = x, Y = y, id = ids)
  # calculate point estimate and CS variance
  df_t = df[df[,tname] == t_star,]
  mod_t_star = single_period_est(x = df_t$R, y = df_t$Y, h, b)
  tau_c = mod_t_star$est[1]
  tau_bc = mod_t_star$est[2]
  V_c_cs = mod_t_star$se[1]^2
  V_bc_cs = mod_t_star$se[2]^2
  for(i in 1:length(t_vec)) {
    df_t = df[df[,tname] == t_vec[i],]
    mod_t= single_period_est(x = df_t$R, y = df_t$Y, h, b)
    tau_c = tau_c - mod_t$est[1] * w_vec[i]
    tau_bc = tau_bc - mod_t$est[2] * w_vec[i]
    V_c_cs = V_c_cs + (mod_t$se[1]^2) * (w_vec[i]^2)
    V_bc_cs = V_bc_cs + (mod_t$se[2]^2) * (w_vec[i]^2)
  }

  # calculate PC variance
  w_vec = c(w_vec, -1)
  t_vec = c(t_vec, t_star)
  V_c_pv = V_c_pc = V_c_cs
  V_bc_pv = V_bc_pc = V_bc_cs
  for(i in 1:length(t_vec)) {
    for(j in 1:length(t_vec)) {
      if(t_vec[i] == t_vec[j]) next
      if(t_vec[i] != t_vec[j]) {
        df1 = df[df[,tname] == t_vec[i],]
        df2 = df[df[,tname] == t_vec[j],]

        x1 = df1$R
        x2 = df2$R
        y1 = df1$Y
        y2 = df2$Y
        ids1 = df1$id
        ids2 = df2$id

        cov12 = covariances_c(x1, x2, y1, y2, h, b, ids1, ids2)
        cov12_bc = covariances_bc(x1, x2, y1, y2, h, b, ids1, ids2)

        C_c_pc = cov12$c_1_2_l_l[1,1] + cov12$c_1_2_r_r[1,1]
        C_bc_pc = cov12_bc$c_1_2_l_l[1,1] + cov12_bc$c_1_2_r_r[1,1]
        C_c_pv = cov12$c_1_2_r_l[1,1] + cov12$c_1_2_l_r[1,1]
        C_bc_pv = cov12_bc$c_1_2_r_l[1,1] + cov12_bc$c_1_2_l_r[1,1]

        V_c_pc = V_c_pc + C_c_pc * w_vec[i] * w_vec[j]
        V_bc_pc = V_bc_pc + C_bc_pc * w_vec[i] * w_vec[j]

        V_c_pv = V_c_pv + (C_c_pc - C_c_pv) * w_vec[i] * w_vec[j]
        V_bc_pv = V_bc_pv + (C_bc_pc - C_bc_pv) * w_vec[i] * w_vec[j]
      }
    }
  }

  data.frame(
    method = c("C", "BC"),
    est = c(tau_c, tau_bc),
    se_cs = c(sqrt(V_c_cs), sqrt(V_bc_cs)),
    se_pc = c(sqrt(V_c_pc), sqrt(V_bc_pc)),
    se_pv = c(sqrt(V_c_pv), sqrt(V_bc_pv))
  )
}
