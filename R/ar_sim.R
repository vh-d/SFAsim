#' @export
sim_ar_gamma <- function(l = 10, 
                         ar_coef = 0.8, 
                         ar_const = 0.2,
                         ar_scale = 0.05, ...) {
  
  y0 <- rgamma(1, 1, 1)
  y <- numeric(l)
  y[1] <- y0
  for (i in 2:l) {
    y[i] <- rgamma(1, shape = max(0, (ar_const + ar_coef * y[i-1])) / ar_scale, scale = ar_scale)
  }
  
  return(y)
}


#' @export
sim_ar_lnorm <- function(l = 10, 
                         ar_coef = 0.8, 
                         ar_const = 0.2,
                         sdlog = 0.2, ...) {
  
  y0 <- rlnorm(1, 0, 1)
  y <- numeric(l)
  y[1] <- y0
  for (i in 2:l) {
    y[i] <- rlnorm(1, meanlog = log(ar_const + ar_coef * y[i-1]) - 0.5*sdlog^2, sdlog = sdlog)
  }
  
  return(y)
}


#' @export
sim_ar_norm <- function(l = 10, 
                        ar_coef = 0.8, 
                        ar_const = 0.2,
                        sd = 0.2, ...) {
  
  y0 <- rlnorm(1, 0, 1)
  y <- numeric(l)
  y[1] <- y0
  for (i in 2:l) {
    y[i] <- rnorm(1, ar_const + ar_coef * y[i-1], sd = sd)
  }
  
  return(y)
}



# test
# plot(sim_ar_gamma(l= 30), ylim = c(0, 3), type = "l")
# plot(sim_ar_lnorm(l= 30), ylim = c(0, 3), type = "l")
# plot(sim_ar_norm(l= 30), ylim = c(0, 3), type = "l")

