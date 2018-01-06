#' generate simulated cross-section data using translog production function
#' @param N number of observations
#' @param frontier_coef coefficients for explanatory variables
#' @param ineff_covar_coef coefficients for exogeneous varaibles
#' @param ineff_sigma variance of the inefficiency term
#' @param frontier_sigma variance of the random noise
#' @param ineff_sign production (-1) or cost inefficiency form (1),
#' @param aslist TRUE/FALSE to return list or array
#' @export
sim_data_cs <- function(N = 500,
                        n_of_inputs = 2,
                        frontier_coef = c(a_0 = 3, 
                                          a_1 = 0.5, 
                                          a_2 = 0.3, 
                                          b_1_1 = -0.1, 
                                          b_2_2 = 0.2, 
                                          b_1_2 = -0.04),
                        ineff_covar_coef = c(0.5, 0.1, -0.1),
                        ineff_sigma = 2,
                        frontier_sigma = 3,
                        ineff_sign  = -1,
                        aslist  = FALSE) {
  
  
  # handle arguments ---------------
  
  # frontier
  front_coef_len  <- length(frontier_coef) - 1
  stopifnot(front_coef_len == n_of_inputs + n_of_inputs*(n_of_inputs + 1)/2)
  
  # inefficiency covariates
  ineff_covar_len <- length(ineff_covar_coef) - 1
  stopifnot(ineff_covar_len == 0)
  
  
  # inputs data ------------------------ 
  X <- matrix(rtnorm(n     = N*n_of_inputs,
                     mean  = 5,
                     sd    = 3,
                     lower = 0),
              nrow = N,
              ncol = n_of_inputs)
  colnames(X) <- paste0("x", 1:n_of_inputs)
  
  # frontier (optimal output)
  frontier <- micEcon::translogCalc(xNames = colnames(X), 
                                    data = as.data.frame(X),
                                    coef = frontier_coef,
                                    dataLogged = FALSE)
  
  est <- suppressWarnings(micEcon::translogEst(yName = "y", 
                                               xNames = colnames(X),
                                               data = cbind(y = frontier, as.data.frame(X))))
  
  # inefficiency data -------------------
  if (ineff_covar_len > 0 ) {
    INEFF_COVAR <- matrix(rnorm(n    = N*ineff_covar_len,
                                mean = 3,
                                sd   = 2),
                          nrow = N,
                          ncol = ineff_covar_len)
    
    colnames(INEFF_COVAR) <- paste0("ineff_covar_", 1:ineff_covar_len)
    
    ineff_term <- sapply(as.vector(cbind(1, INEFF_COVAR) %*% ineff_covar_coef),
                         function(x) {
                           rtnorm(n     = 1,
                                  mean  = x,
                                  lower = 0,
                                  sd    = ineff_sigma)
                         }
    )
  } else {
    ineff_term <- rgamma(N,
                         shape = ineff_covar_coef/ineff_sigma,
                         scale = ineff_sigma)
    
    # ineff_term <- rtnorm(N,
    #                      mean = ineff_covar_coef,
    #                      lower = 0,
    #                      sd = ineff_sigma)
    Z <- NULL
  }
  
  
  # error terms ----------------------
  
  # symmetric error term
  frontier_error_term <- rnorm(N,
                               mean = 0,
                               sd = frontier_sigma)
  
  # composite error term
  total_error_term <- ineff_sign*ineff_term + frontier_error_term
  
  # real output
  y <- log(frontier) + total_error_term
  
  
  if (aslist) {
    return(list(y = y, 
                X = est$model.matrix, 
                Z = Z, 
                frontier = log(frontier), 
                ineff_term = ineff_term, 
                frontier_error_term = frontier_error_term, 
                total_error_term = total_error_term))
  } else
    return(cbind(y = y, 
                 X = est$model.matrix, 
                 Z, 
                 ineff_term, 
                 frontier_error_term, 
                 total_error_term))
}




#' generate simulated panel data for testing purposes
#' @param k number of individuals - cross section
#' @param t number of observations per individual
#' @param frontier_coef coefficients for explanatory variables
#' @param ineff_covar_coef coefficients for exogeneous varaibles
#' @param ineff_sigma variance of the inefficiency term
#' @param frontier_sigma variance of the random noise
#' @param ineff_sign production (-1) or cost inefficiency form (1),
#' @param aslist TRUE/FALSE to return list or array
#' @export
sim_data_panel <- function(k = 20, # number of individuals - cross section
                           t = 10, # number of observations per individual
                           x_mean = 10,
                           x_sd = 1,
                           z_mean = 0.5,
                           z_sd = 1,
                           frontier_coef = c(10, 6, 3),
                           z_intercept = 1,
                           ineff_covar_coef = c(0.1, -0.4),
                           ineff_sigma = 2, # variance of the inefficiency term
                           frontier_sigma = 3, # variance of the random noise
                           ineff_sign = -1,
                           aslist = FALSE) {
  
  N <- t * k
  
  x_ncols <- length(frontier_coef) - 1
  ineff_covar_len <- length(ineff_covar_coef)
  
  # random explanatory data
  X <- matrix(rtnorm(n     = N*x_ncols,
                     mean  = x_mean,
                     sd    = z_sd,
                     lower = 0),
              nrow = N,
              ncol = x_ncols)
  
  colnames(X) <- paste0("x", 1:x_ncols)
  
  # random exogeneous data
  
  # intercept for each individual (cross section)
  mu_unique <- rtnorm(n = k, mean = z_intercept, sd = 1, lower = 0)
  mu <- rep(mu_unique, each = t)
  K <- cbind(k = rep(1:k, each = t), t = rep(1:t, times = k))
  
  if (ineff_covar_len > 0 ) {
    Z <- matrix(rnorm(n    = N*ineff_covar_len,
                      mean = z_mean,
                      sd   = z_sd),
                nrow = N,
                ncol = ineff_covar_len)
    
    colnames(Z) <- paste0("z", 1:ineff_covar_len)
    
    ineff_term <- sapply(mu + as.vector(as.vector(Z %*% ineff_covar_coef)),
                         function(x) {
                           rtnorm(n     = 1,
                                  mean  = x,
                                  lower = 0,
                                  sd    = ineff_sigma)
                         }
    )
    
  } else {
    
    ineff_term <- rtnorm(n     = N,
                         mean  = mu,
                         lower = 0,
                         sd    = ineff_sigma)
    Z <- NULL
  }
  
  # error terms
  frontier_error_term <- rnorm(n    = N,
                               mean = 0,
                               sd   = frontier_sigma)
  
  total_error_term <- ineff_sign*ineff_term + frontier_error_term
  
  y <- as.vector( cbind(1, X) %*% frontier_coef) + total_error_term
  
  if (aslist) {
    return(list(y = y, X = X, Z = Z, K = K,
                ineff_term = ineff_term, frontier_error_term = frontier_error_term, total_error_term = total_error_term,
                N = N, t = t, k = k, mu = mu_unique))
  } else
    return(cbind(K, y, X, Z, ineff_term, frontier_error_term, total_error_term))
}

# inputs data ---------------------------
inputs_random <- function(t, k, m, 
                          mean = 0, 
                          sd = 1, 
                          intercept = TRUE, 
                          prefix = "x", 
                          ...) {
  
  N <- t*k
  # random explanatory data
  X <- cbind(
    if (intercept) 1 else NULL,
    matrix(rnorm(n     = N*(m-as.integer(intercept)),
                 mean  = mean,
                 sd    = sd),
           nrow = N,
           ncol = m-as.integer(intercept)))
  
  colnames(X) <- c( if (intercept) paste0(prefix, "_intercept") else NULL, paste0(prefix, 1:(m-as.integer(intercept))))
  
  return(X)
}

# ineff determinants data funcionts ------------------
deter_random <- function(...) {
  inputs_random(..., prefix = "d")
}

# ineff term data functions -------------------
ineff_ar <- function(k, t, 
                     innov = NULL,
                     start_dist  = "gamma",
                     start_dist_opts = list(shape = 1, rate = 2),
                     ar_dist      = c("gamma", "lnorm", "norm"),
                     ar_dist_args = list(shape = 5),
                     ar_const    = 0.3,
                     ar_coef     = 0.7, 
                     ...) {
  
  ar_dist <- match.arg(ar_dist)
  
  # starting value for each individual (cross section)
  start_dist_opts_alt <- start_dist_opts
  start_dist_opts_alt$n = k
  start_ineff <- do.call(what = paste0("r", start_dist),
                         args = start_dist_opts_alt)
  
  ar_sim_args <-
    c(
      list(
        dist = ar_dist, 
        l = t, 
        ar_coef  = ar_coef,
        ar_const = ar_const
      ),
      ar_dist_args
    )
  
  ineff_term <- 
    as.vector(
      sapply(1:k,
             FUN =  
               function(k_i) {
                 ar_sim_args$y0 = start_ineff[k_i]
                 if (!is.null(innov)) ar_sim_args$innov = innov[((k_i-1)*t):((k_i)*t)]
                 do.call(what = "ar_sim", args = ar_sim_args)
               }))
  
  return(ineff_term)
}


deter_ar <- function(m, intercept = TRUE, ...) {
  ret <- sapply(1:(m-as.integer(intercept)), function(i) ineff_ar(...))
  if (intercept) ret <- cbind(1, ret)
  
  colnames(ret) <- c(if (intercept) "d_intercept" else NULL, paste0("d",1:(m-as.integer(intercept)) ))
  
  ret
}

inputs_ar <- function(m, intercept = TRUE, ...) {
  ret <- sapply(1:(m-as.integer(intercept)), function(i) ineff_ar(...))
  if (intercept) ret <- cbind(1, ret)
  
  colnames(ret) <- c(if (intercept) "x_intercept" else NULL, paste0("x",1:(m-as.integer(intercept)) ))
  
  ret
}

#' generate simulated (cross-section or balanced panel) data for testing purposes
#' @param k number of individuals (cross section dimension)
#' @param t number of observations per individual (time dimension)
#' @param ineff_sign production (-1) or cost inefficiency form (1),
#' @param frontier_coef coefficients for explanatory variables
#' @param inputs model for inputs data
#' @param inputs_params parameters for model for inputs data
#' @param inefficiency model for inefficiency
#' @param ineff_params parameters for model for inefficiency
#' @param determinants model for inefficiency exogeneous determinants
#' @param deter_params parameters for model for inefficiency exogeneous determinants
#' @param deter_coef coefficients on exogeneous determinants
#' @param ineff_covar_coef coefficients for exogeneous varaibles
#' @param ineff_sigma variance of the inefficiency term
#' @param frontier_sigma sqrt of variance of the symmetric error
#' @param aslist TRUE/FALSE to return list or array
#' @export
sim_data <- function(k              = 40, # number of individuals - cross section
                     t              = 10, # number of observations per individual
                     ineff_sign     = -1,
                     inputs         = c("random", "fe", "ar"),
                     inputs_params  = list(),
                     frontier_coef  = c(10, 0.5, -0.1),
                     inefficiency   = c("random", "fe", "ar"),
                     ineff_params   = list(),
                     frontier_sigma = 1.5, # variance of the random noise
                     determinants   = c("none", "random", "fe", "ar"),
                     deter_params   = list(),
                     deter_coef     = c(0.1, 0.4),
                     aslist         = FALSE) {
  
  inputs       <- match.arg(inputs)
  inefficiency <- match.arg(inefficiency)
  determinants <- match.arg(determinants)
  
  inputs_params$k <- k
  inputs_params$t <- t
  
  ineff_params$k <- k
  ineff_params$t <- t
  
  deter_params$k <- k
  deter_params$t <- t
  
  # total number of observations = n of periods * n of entities 
  N <- t * k
  
  # prepare inputs data
  inputs_params[["m"]] <- length(frontier_coef)
  X <- do.call(what = paste0("inputs_", inputs), 
               args = inputs_params)
  
  # prepare determinants of inefficiency (if requested)
  deter_params[["m"]] <- length(deter_coef)
  Z <- if (determinants == "none") matrix(0, N, 1) else 
    do.call(what = paste0("deter_", determinants), 
            args = deter_params)
  
  # prepare inefficiency term
  ineff_params[["innov"]] <- Z %*% deter_coef
  ineff_term <- do.call(paste0("ineff_", inefficiency), args = ineff_params)
  
  # matrix of panel data indices
  K <- cbind(k = rep(1:k, each = t), t = rep(1:t, times = k))
  
  # error terms
  frontier_error_term <- rnorm(n    = N,
                               mean = 0,
                               sd   = frontier_sigma)
  
  total_error_term <- ineff_sign*ineff_term + frontier_error_term
  
  yhat <- as.vector( X %*% frontier_coef) 
  y <- yhat + total_error_term
  
  # return list or matrix as requested
  if (aslist) {
    return(list(y = y, yhat = yhat, X = X, Z = Z, K = K,
                ineff_term = ineff_term,
                frontier_error_term = frontier_error_term,
                total_error_term = total_error_term,
                N = N, t = t, k = k))
  } else
    return(cbind(K, y, yhat, X, Z, ineff_term, frontier_error_term, total_error_term))
}


sim_data_cs_homo <- function(N, ineff_sign = -1) {
  
  x1 <- rnorm(N, 8)
  x2 <- rnorm(N, 15)
  
  ineff_term <- rtnorm(N, mean = 2, lower = 1, sd = 3)
  frontier_error_term <- rnorm(N, 0, 1)
  total_error_term <- ineff_sign*ineff_term + frontier_error_term
  
  y <- 10 + 6*x1 + 3*x2 + total_error_term
  
  return(cbind(y, x1, x2, z1, z2, ineff_term, frontier_error_term, total_error_term))
}

#' draw random values from truncated normal distribution
#' @export
rtnorm <- function(n = 1, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  rtnorm_cpp(n, mean, sd, lower, upper)
}
