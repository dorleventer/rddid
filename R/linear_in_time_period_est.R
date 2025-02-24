#' RDDID Estimation
#'
#' This function estimates the ATU or ATT using a RDDID framework, assuming potential outcome discontinuities evolve linearly across time.
#'
#' @param df dataframe
#' @param t_star time period with RD treatment
#' @param h bandwidth
#' @param b bias correction bandwidth
#' @param t_vec vector of time periods to estimate linear trend over
#' @param tname the name of the time period variable in df
#' @param xname the name of the running variable in df
#' @param yname the name of the outcome variable in df
#' @param idname the name of the id variable in df
#'
#' @return Conventional (C) and Bias Corrected (BC) treatment effect estimates and standard errors (SE). SEs are provided for sampling under repeated cross section (CS), panel with constant running variable (PC) and varying running variable (PV).
#' @export
linear_in_time_period_est = function(df, t_star, h, b, t_vec, tname = "time", xname = "R", yname = "Y", idname = "id") {
  # prep the data
  time = df[,tname]
  x = df[,xname]
  y = df[,yname]
  ids = df[,idname]
  df = data.frame(time, R = x, Y = y, id = ids)
  # calculate single time period models
  df_1 = df[df$time == t_vec[1],]
  mod_1 = single_period_est(x = df_1$R, y = df_1$Y, h, b)
  df_2 = df[df$time == t_vec[2],]
  mod_2 = single_period_est(x = df_2$R, y = df_2$Y, h, b)
  df_t_star = df[df$time == t_star,]
  mod_t_star = single_period_est(df_t_star$R, df_t_star$Y, h, b)

  # calculate slope and point estimate
  m_c = mod_2$est[1] - mod_1$est[1]
  m_bc = mod_2$est[2] - mod_1$est[2]

  alpha_t_star_c = (t_star - t_vec[2])*m_c + mod_2$est[1]
  alpha_t_star_bc = (t_star - t_vec[2])*m_bc + mod_2$est[2]

  tau_c = mod_t_star$est[1] - alpha_t_star_c
  tau_bc = mod_t_star$est[2] - alpha_t_star_bc

  # CS variance
  V_c_cs = mod_t_star$se[1]^2 +
    ((t_star - t_vec[1])^2) * mod_2$se[1]^2 +
    ((t_star - t_vec[2])^2) * mod_1$se[1]^2
  V_bc_cs = mod_t_star$se[2]^2 +
    ((t_star - t_vec[1])^2) * mod_2$se[2]^2 +
    ((t_star - t_vec[2])^2) * mod_1$se[2]^2

  # PC and PV variance

  cov_t_star_1 = covariances_c(df_t_star$R, df_1$R, df_t_star$Y, df_1$Y, h, b, df_t_star$id, df_1$id)
  cov_t_star_1_bc = covariances_bc(df_t_star$R, df_1$R, df_t_star$Y, df_1$Y, h, b, df_t_star$id, df_1$id)

  cov_t_star_2 = covariances_c(df_t_star$R, df_2$R, df_t_star$Y, df_2$Y, h, b, df_t_star$id, df_2$id)
  cov_t_star_2_bc = covariances_bc(df_t_star$R, df_2$R, df_t_star$Y, df_2$Y, h, b, df_t_star$id, df_2$id)

  cov_1_2 = covariances_c(df_1$R, df_2$R, df_1$Y, df_2$Y, h, b, df_1$id, df_2$id)
  cov_1_2 = covariances_bc(df_1$R, df_2$R, df_1$Y, df_2$Y, h, b, df_1$id, df_2$id)

  C_c_pc_t_star_1 = cov_t_star_1$c_1_2_l_l[1,1] + cov_t_star_1$c_1_2_r_r[1,1]
  C_c_pc_t_star_2 = cov_t_star_2$c_1_2_l_l[1,1] + cov_t_star_2$c_1_2_r_r[1,1]
  C_c_pc_1_2 = cov_1_2$c_1_2_l_l[1,1] + cov_1_2$c_1_2_r_r[1,1]

  C_bc_pc_t_star_1 = cov_t_star_1_bc$c_1_2_l_l[1,1] + cov_t_star_1_bc$c_1_2_r_r[1,1]
  C_bc_pc_t_star_2 = cov_t_star_2_bc$c_1_2_l_l[1,1] + cov_t_star_2_bc$c_1_2_r_r[1,1]
  C_bc_pc_1_2 = cov_1_2_bc$c_1_2_l_l[1,1] + cov_1_2_bc$c_1_2_r_r[1,1]

  C_c_pv_t_star_1 = cov_t_star_1$c_1_2_r_l[1,1] + cov_t_star_1$c_1_2_l_r[1,1]
  C_c_pv_t_star_2 = cov_t_star_2$c_1_2_r_l[1,1] + cov_t_star_2$c_1_2_l_r[1,1]
  C_c_pv_1_2 = cov_1_2$c_1_2_r_l[1,1] + cov_1_2$c_1_2_l_r[1,1]

  C_bc_pv_t_star_1 = cov_t_star_1_bc$c_1_2_r_l[1,1] + cov_t_star_1_bc$c_1_2_l_r[1,1]
  C_bc_pv_t_star_2 = cov_t_star_2_bc$c_1_2_r_l[1,1] + cov_t_star_2_bc$c_1_2_l_r[1,1]
  C_bc_pv_1_2 = cov_1_2_bc$c_1_2_r_l[1,1] + cov_1_2_bc$c_1_2_l_r[1,1]

  V_c_pc = V_c_cs +
    2*C_c_pc_t_star_1*(t_star - t_vec[2]) +
    (-2)*C_c_pc_t_star_2*(t_star - t_vec[1]) +
    (-2)*C_c_pc_1_2*(t_star - t_vec[1])*(t_star - t_vec[2])
  V_bc_pc = V_bc_cs +
    2*C_bc_pc_t_star_1*(t_star - t_vec[2]) +
    (-2)*C_bc_pc_t_star_2*(t_star - t_vec[1]) +
    (-2)*C_bc_pc_1_2*(t_star - t_vec[1])*(t_star - t_vec[2])
  V_c_pv = V_c_pc +
    (-2)*C_c_pv_t_star_1*(t_star - t_vec[2]) +
    (2)*C_c_pv_t_star_2*(t_star - t_vec[1]) +
    (2)*C_c_pv_1_2*(t_star - t_vec[1])*(t_star - t_vec[2])
  V_bc_pv = V_bc_pc +
    (-2)*C_bc_pv_t_star_1*(t_star - t_vec[2]) +
    (2)*C_bc_pv_t_star_2*(t_star - t_vec[1]) +
    (2)*C_bc_pv_1_2*(t_star - t_vec[1])*(t_star - t_vec[2])

  data.frame(
    method = c("C", "BC"),
    est = c(tau_c, tau_bc),
    se_cs = c(sqrt(V_c_cs), sqrt(V_bc_cs)),
    se_pc = c(sqrt(V_c_pc), sqrt(V_bc_pc)),
    se_pv = c(sqrt(V_c_pv), sqrt(V_bc_pv))
  )
}
