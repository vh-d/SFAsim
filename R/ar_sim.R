#' @export
#' @rdname ar_sim
ar_sim_gamma <- function(l = 10, 
                         y0 = rgamma(1, 1, 1),
                         ar_coef = 0.8, 
                         ar_const = 0.2,
                         innov = numeric(l),
                         shape = 3) {
  
  y <- numeric(l)
  y[1] <- y0
  
  for (i in 2:l) {
    y[i] <- rgamma(1, shape = shape, rate = max(0.001, shape/(ar_const + ar_coef * y[i-1] + innov[i])))
  }
  
  return(y)
}


#' @export
#' @rdname ar_sim
ar_sim_lnorm <- function(l = 10, 
                         y0 = rlnorm(1, 0, 1),
                         ar_coef = 0.8, 
                         ar_const = 0.2,
                         innov = numeric(l),
                         sdlog = 0.2, 
                         ...) {
  
  y <- numeric(l)
  y[1] <- y0
  
  for (i in 2:l) {
    y[i] <- rlnorm(1, meanlog = log(ar_const + ar_coef * y[i-1] + innov[i]) - 0.5*sdlog^2, sdlog = sdlog)
  }
  
  return(y)
}


#' @export 
#' @rdname ar_sim
ar_sim_norm <- function(l = 10, 
                        y0 = rlnorm(1, 0, 1),
                        ar_coef = 0.8, 
                        ar_const = 0.2,
                        innov = numeric(l),
                        ...) {
  
  y <- numeric(l)
  y[1] <- y0
  
  for (i in 2:l) {
    y[i] <- rnorm(1, ar_const + ar_coef * y[i-1] + innov[i], ...)
  }
  
  return(y)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dist PARAM_DESCRIPTION, Default: c("norm", "gamma", "lnorm")
#' @param l PARAM_DESCRIPTION, Default: 10
#' @param y0 PARAM_DESCRIPTION, Default: rlnorm(1, 0, 1)
#' @param innov_data PARAM_DESCRIPTION, Default: NULL
#' @param innov_coef PARAM_DESCRIPTION, Default: NULL
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @export 
#' @rdname ar_sim
ar_sim <- function(dist = c("norm", "gamma", "lnorm"),
                   ...){
  dist <- match.arg(dist)
  do.call(paste0("ar_sim_", dist), args = list(...))
}

