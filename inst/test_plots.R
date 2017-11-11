# test
require(VHtools)
# plot(ar_sim_gamma(l= 30), ylim = c(0, 3), type = "l")
# plot(ar_sim_lnorm(l= 30), ylim = c(0, 3), type = "l")
# plot(ar_sim_norm(l= 30), ylim = c(0, 3), type = "l")
len <- 30

covariate <- 
  cbind(
    ar_sim(l = len + 1) + 0.2*(1:(len + 1)),
    ar_sim(l = len + 1) +     (1:(len + 1)))
    

ar <- ar_sim(dist = "gamma", 
             l = len,
             ar_coeff = 0.9,
             innov_data = diff(covariate),
             innov_coeffs = c(0.1, 0.3),
             scale = 0.01)

plot(as.ts(ar), ylim = c(0, 3), type = "l")
lines(as.ts(1 + rescale(covariate[-1], -1, 1)), col = "red")

cor(ar, covariate[-1,])
lm(ar~covariate[-1, ])

arima(ar, order = c(1, 0, 0), xreg = covariate[-1, ])
