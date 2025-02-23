#' RDDID Estimation
#'
#' This function estimates the ATU or ATT using a RDDID framework, assuming potential outcome discontinuities evolve linearly across time.
#'
#' @param df dataframe
#' @param t_star time period with RD treatment
#' @param h bandwidth
#' @param b bias correction bandwidth
#'
#' @return Conventional (C) and Bias Corrected (BC) treatment effect estimates and standard errors (SE). SEs are provided for sampling under repeated cross section (CS), panel with constant running variable (PC) and varying running variable (PV).
#' @export
linear_in_time_period_est = function(df, t_star, h, b) {

  # calculate single time period models
  df_1999 = df |> filter(year == 1999)
  mod_1999 = single_period_est(x = df_1999$R, y = df_1999$Y, h, b)
  df_2000 = df |> filter(year == 2000)
  mod_2000 = single_period_est(x = df_2000$R, y = df_2000$Y, h, b)
  df_t_star = df |> filter(year == t_star)
  mod_t_star = single_period_est(df_t_star$R, df_t_star$Y, h, b)

  # calculate slope and point estimate
  m_c = mod_2000$est[1] - mod_1999$est[1]
  m_bc = mod_2000$est[2] - mod_1999$est[2]

  alpha_t_star_c = (t_star - 2000)*m_c + mod_2000$est[1]
  alpha_t_star_bc = (t_star - 2000)*m_bc + mod_2000$est[2]

  tau_c = mod_t_star$est[1] - alpha_t_star_c
  tau_bc = mod_t_star$est[2] - alpha_t_star_bc

  # CS variance
  V_c_cs = mod_t_star$se[1]^2 +
    ((t_star - 1999)^2) * mod_2000$se[1]^2 +
    ((t_star - 2000)^2) * mod_1999$se[1]^2
  V_bc_cs = mod_t_star$se[2]^2 +
    ((t_star - 1999)^2) * mod_2000$se[2]^2 +
    ((t_star - 2000)^2) * mod_1999$se[2]^2

  # PC and PV variance

  cov_t_star_1999 = covariances_c(df_t_star$R, df_1999$R, df_t_star$Y, df_1999$Y, h, b, df_t_star$id, df_1999$id)
  cov_t_star_1999_bc = covariances_bc(df_t_star$R, df_1999$R, df_t_star$Y, df_1999$Y, h, b, df_t_star$id, df_1999$id)

  cov_t_star_2000 = covariances_c(df_t_star$R, df_2000$R, df_t_star$Y, df_2000$Y, h, b, df_t_star$id, df_2000$id)
  cov_t_star_2000_bc = covariances_bc(df_t_star$R, df_2000$R, df_t_star$Y, df_2000$Y, h, b, df_t_star$id, df_2000$id)

  cov_1999_2000 = covariances_c(df_1999$R, df_2000$R, df_1999$Y, df_2000$Y, h, b, df_1999$id, df_2000$id)
  cov_1999_2000_bc = covariances_bc(df_1999$R, df_2000$R, df_1999$Y, df_2000$Y, h, b, df_1999$id, df_2000$id)

  C_c_pc_t_star_1999 = cov_t_star_1999$c_1_2_l_l[1,1] + cov_t_star_1999$c_1_2_r_r[1,1]
  C_c_pc_t_star_2000 = cov_t_star_2000$c_1_2_l_l[1,1] + cov_t_star_2000$c_1_2_r_r[1,1]
  C_c_pc_1999_2000 = cov_1999_2000$c_1_2_l_l[1,1] + cov_1999_2000$c_1_2_r_r[1,1]

  C_bc_pc_t_star_1999 = cov_t_star_1999_bc$c_1_2_l_l[1,1] + cov_t_star_1999_bc$c_1_2_r_r[1,1]
  C_bc_pc_t_star_2000 = cov_t_star_2000_bc$c_1_2_l_l[1,1] + cov_t_star_2000_bc$c_1_2_r_r[1,1]
  C_bc_pc_1999_2000 = cov_1999_2000_bc$c_1_2_l_l[1,1] + cov_1999_2000_bc$c_1_2_r_r[1,1]

  C_c_pv_t_star_1999 = cov_t_star_1999$c_1_2_r_l[1,1] + cov_t_star_1999$c_1_2_l_r[1,1]
  C_c_pv_t_star_2000 = cov_t_star_2000$c_1_2_r_l[1,1] + cov_t_star_2000$c_1_2_l_r[1,1]
  C_c_pv_1999_2000 = cov_1999_2000$c_1_2_r_l[1,1] + cov_1999_2000$c_1_2_l_r[1,1]

  C_bc_pv_t_star_1999 = cov_t_star_1999_bc$c_1_2_r_l[1,1] + cov_t_star_1999_bc$c_1_2_l_r[1,1]
  C_bc_pv_t_star_2000 = cov_t_star_2000_bc$c_1_2_r_l[1,1] + cov_t_star_2000_bc$c_1_2_l_r[1,1]
  C_bc_pv_1999_2000 = cov_1999_2000_bc$c_1_2_r_l[1,1] + cov_1999_2000_bc$c_1_2_l_r[1,1]

  V_c_pc = V_c_cs +
    2*C_c_pc_t_star_1999*(t_star - 2000) +
    (-2)*C_c_pc_t_star_2000*(t_star - 1999) +
    (-2)*C_c_pc_1999_2000*(t_star - 1999)*(t_star - 2000)
  V_bc_pc = V_bc_cs +
    2*C_bc_pc_t_star_1999*(t_star - 2000) +
    (-2)*C_bc_pc_t_star_2000*(t_star - 1999) +
    (-2)*C_bc_pc_1999_2000*(t_star - 1999)*(t_star - 2000)
  V_c_pv = V_c_pc +
    (-2)*C_c_pv_t_star_1999*(t_star - 2000) +
    (2)*C_c_pv_t_star_2000*(t_star - 1999) +
    (2)*C_c_pv_1999_2000*(t_star - 1999)*(t_star - 2000)
  V_bc_pv = V_bc_pc +
    (-2)*C_bc_pv_t_star_1999*(t_star - 2000) +
    (2)*C_bc_pv_t_star_2000*(t_star - 1999) +
    (2)*C_bc_pv_1999_2000*(t_star - 1999)*(t_star - 2000)

  tibble(
    method = c("C", "BC"),
    est = c(tau_c, tau_bc),
    se_cs = c(sqrt(V_c_cs), sqrt(V_bc_cs)),
    se_pc = c(sqrt(V_c_pc), sqrt(V_bc_pc)),
    se_pv = c(sqrt(V_c_pv), sqrt(V_bc_pv))
  )
}
