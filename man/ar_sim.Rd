% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ar_sim.R
\name{ar_sim_gamma}
\alias{ar_sim_gamma}
\alias{ar_sim_lnorm}
\alias{ar_sim_norm}
\alias{ar_sim}
\title{FUNCTION_TITLE}
\usage{
ar_sim_gamma(l = 10, y0 = rgamma(1, 1, 1), ar_coef = 0.8,
  ar_const = 0.2, innov = numeric(l), shape = 3)

ar_sim_lnorm(l = 10, y0 = rlnorm(1, 0, 1), ar_coef = 0.8,
  ar_const = 0.2, innov = numeric(l), sdlog = 0.2, ...)

ar_sim_norm(l = 10, y0 = rlnorm(1, 0, 1), ar_coef = 0.8, ar_const = 0.2,
  innov = numeric(l), ...)

ar_sim(dist = c("norm", "gamma", "lnorm"), ...)
}
\arguments{
\item{l}{PARAM_DESCRIPTION, Default: 10}

\item{y0}{PARAM_DESCRIPTION, Default: rlnorm(1, 0, 1)}

\item{...}{PARAM_DESCRIPTION}

\item{dist}{PARAM_DESCRIPTION, Default: c("norm", "gamma", "lnorm")}

\item{innov_data}{PARAM_DESCRIPTION, Default: NULL}

\item{innov_coef}{PARAM_DESCRIPTION, Default: NULL}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
DETAILS
}
