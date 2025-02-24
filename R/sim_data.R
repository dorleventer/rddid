#' The Data Generating Process for RDDID
#'
#' @param dgp_type CS, PC, PV
#' @param seed random seed
#' @param params dgp parameters
#'
#' @return A data-frame of simulated RDDID data with 2 periods based on the Grembi et al (2016) application.
simulate_dgp <- function(dgp_type, seed = 1, params = list()) {
  set.seed(seed)

  if (dgp_type == "cs") {
    df <- data.frame(
      id = rep(1:params$n, times = 2),
      time = rep(c(1, 2), each = params$n),
      R = (stats::rbeta(n = 2 * params$n, shape1 = 2, shape2 = 4) - 0.375) * 5000
    ) |>
      dplyr::mutate(unit_fixed_effect = stats::rnorm(2 * params$n, params$mean_FE, params$sd_FE))
  }
  if (dgp_type == "pc") {
    df <- data.frame(
      id = rep(1:params$n, times = 2),
      time = rep(c(1, 2), each = params$n),
      R = rep((stats::rbeta(n = params$n, shape1 = 2, shape2 = 4) - 0.375) * 5000, times = 2)
    ) %>%
      dplyr::mutate(unit_fixed_effect = stats::rnorm(2 * params$n, params$mean_FE, params$sd_FE)) %>%
      dplyr::arrange(id, time) %>%
      dplyr::mutate(unit_fixed_effect = ifelse(time == 2, dplyr::lag(unit_fixed_effect), unit_fixed_effect))
  }
  if (dgp_type == "pv") {
    R_t1 <- (stats::rbeta(n = params$n, shape1 = 2, shape2 = 4) - 0.375) * 5000
    R_t2 <- 0.97*R_t1 + stats::rnorm(params$n, 153, 410)

    df <- data.frame(
      id = rep(1:params$n, times = 2),
      time = rep(c(1, 2), each = params$n),
      R = c(R_t1, R_t2)
    ) %>%
      dplyr::mutate(unit_fixed_effect = stats::rnorm(2 * params$n, params$mean_FE, params$sd_FE)) %>%
      dplyr::arrange(id, time) %>%
      dplyr::mutate(unit_fixed_effect = ifelse(time == 2, dplyr::lag(unit_fixed_effect), unit_fixed_effect))
  }

  df = df |>
    dplyr::mutate(
      u = stats::rnorm(2 * params$n, 0, params$se_error),
      mu = dplyr::case_when(
        time == 1 & R >= 0 ~ R * params$mu1_t1_right + R^2 * params$mu2_t1_right + R^3 * params$mu3_t1_right + R^4 * params$mu4_t1_right + R^5 * params$mu5_t1_right,
        time == 1 & R < 0 ~ R * params$mu1_t1_left + R^2 * params$mu2_t1_left + R^3 * params$mu3_t1_left + R^4 * params$mu4_t1_left + R^5 * params$mu5_t1_left,
        time == 2 & R >= 0 ~ R * params$mu1_t2_right + R^2 * params$mu2_t2_right + R^3 * params$mu3_t2_right + R^4 * params$mu4_t2_right + R^5 * params$mu5_t2_right,
        time == 2 & R < 0 ~ R * params$mu1_t2_left + R^2 * params$mu2_t2_left + R^3 * params$mu3_t2_left + R^4 * params$mu4_t2_left + R^5 * params$mu5_t2_left
      ),
      Y = ifelse(
        time == 1,
        mu + params$jump_t1*(R >= 0) + params$time_fixed_effect_t1 + unit_fixed_effect + u,
        mu + params$jump_t2*(R >= 0) + params$time_fixed_effect_t2 + unit_fixed_effect + u
      )
    )

  return(df |> dplyr::select(id, time, R, Y))
}
#' A simulated data-frame
#'
#' @param n number of observations
#' @param seed random seed
#' @param dgp_type CS, PC, PV
#' @param jump_t1,jump_t2 Outcome discontinuity in time periods 1 and 2
#' @export
#' @return A data-frame of simulated RDDID data with 2 periods based on the Grembi et al (2016) application.
sim_data = function(n, seed, dgp_type, jump_t1 = 63, jump_t2 = -63) {

  se_error = 40
  time_fixed_effect_t1 <- -46
  time_fixed_effect_t2 <- 0
  mean_FE <- 116
  sd_FE <- 155

  data("grembi_parameters")
  grembi_parameters =  grembi_parameters|>
    dplyr::filter(deriv != 0) |>
    dplyr::arrange(deriv, side, time)
  for(t in 1:2) {
    for(s in c("left", "right")) {
      for(v in 1:5) {
        temp = grembi_parameters |>
          dplyr::filter(time == t, side == s, deriv == v) |>
          dplyr::pull(estimate)
        name = glue::glue("mu{v}_t{t}_{s}")
        assign(name, temp)
      }
    }
  }

  dgp_parameters <- list(
    n = n,
    se_error = se_error,
    jump_t1 = jump_t1,
    jump_t2 = jump_t2,
    mean_FE = mean_FE,
    sd_FE= sd_FE,
    time_fixed_effect_t2 = time_fixed_effect_t2,
    time_fixed_effect_t1 = time_fixed_effect_t1,
    mu1_t1_left = mu1_t1_left,
    mu1_t1_right = mu1_t1_right,
    mu2_t1_left = mu2_t1_left,
    mu2_t1_right = mu2_t1_right,
    mu3_t1_left = mu3_t1_left,
    mu3_t1_right = mu3_t1_right,
    mu4_t1_left = mu4_t1_left,
    mu4_t1_right = mu4_t1_right,
    mu5_t1_left = mu5_t1_left,
    mu5_t1_right = mu5_t1_right,
    mu1_t2_left = mu1_t2_left,
    mu1_t2_right = mu1_t2_right,
    mu2_t2_left = mu2_t2_left,
    mu2_t2_right = mu2_t2_right,
    mu3_t2_left = mu3_t2_left,
    mu3_t2_right = mu3_t2_right,
    mu4_t2_left = mu4_t2_left,
    mu4_t2_right = mu4_t2_right,
    mu5_t2_left = mu5_t2_left,
    mu5_t2_right = mu5_t2_right
  )

  data = simulate_dgp(dgp_type = dgp_type, seed = seed, params = dgp_parameters)
  return(data)
}
