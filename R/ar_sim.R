#' @export
ar_sim_gamma <- function(l = 10, 
                         ar_coef = 0.8, 
                         ar_const = 0.2,
                         scale = 0.05,
                         innov = numeric(l),
                         ...) {
  
  y0 <- rgamma(1, 1, 1)
  y <- numeric(l)
  y[1] <- y0
  
  for (i in 2:l) {
    y[i] <- rgamma(1, shape = max(0, (ar_const + ar_coef * y[i-1] + innov[i])) / scale, scale = scale)
  }
  
  return(y)
}


#' @export
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

#' @export
ar_sim <- function(dist = c("norm", "gamma", "lnorm"),
                   l = 10,
                   y0 = rlnorm(1, 0, 1),
                   innov_data = NULL,
                   innov_coeffs = NULL,
                   ...){
  dist <- match.arg(dist)

  innov <- if (!is.null(innov_data)) as.vector(as.matrix(innov_data) %*% innov_coeffs) else numeric(l)
  
  switch (dist,
    norm = ar_sim_norm(l = l,   innov = innov, ...),
    lnorm = ar_sim_lnorm(l = l, innov = innov, ...),
    gamma = ar_sim_gamma(l = l, innov = innov, ...)
  )
}

