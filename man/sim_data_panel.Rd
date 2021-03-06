% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{sim_data_panel}
\alias{sim_data_panel}
\title{generate simulated panel data for testing purposes}
\usage{
sim_data_panel(k = 20, t = 10, x_mean = 10, x_sd = 1, z_mean = 0.5,
  z_sd = 1, frontier_coef = c(10, 6, 3), z_intercept = 1,
  ineff_covar_coef = c(0.1, -0.4), ineff_sigma = 2, frontier_sigma = 3,
  ineff_sign = -1, aslist = FALSE)
}
\arguments{
\item{k}{number of individuals - cross section}

\item{t}{number of observations per individual}

\item{frontier_coef}{coefficients for explanatory variables}

\item{ineff_covar_coef}{coefficients for exogeneous varaibles}

\item{ineff_sigma}{variance of the inefficiency term}

\item{frontier_sigma}{variance of the random noise}

\item{ineff_sign}{production (-1) or cost inefficiency form (1),}

\item{aslist}{TRUE/FALSE to return list or array}
}
\description{
generate simulated panel data for testing purposes
}
