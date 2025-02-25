% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_data.R
\name{simulate_dgp}
\alias{simulate_dgp}
\title{The Data Generating Process for RDDID}
\usage{
simulate_dgp(dgp_type, seed = 1, params = list())
}
\arguments{
\item{dgp_type}{CS, PC, PV}

\item{seed}{random seed}

\item{params}{dgp parameters}
}
\value{
A data-frame of simulated RDDID data with 2 periods based on the Grembi et al (2016) application.
}
\description{
The Data Generating Process for RDDID
}
