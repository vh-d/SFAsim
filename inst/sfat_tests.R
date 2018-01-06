library(ggplot2)
library(data.table)
library(SFAt)
library(SFAsim)
library(lmtest)
library(plm)


# CROSS-SECTION DATA ------------------------------------------------------

test_data_cs <- sim_data_cs(N = 1000, frontier_sigma = 1, z_coef = c(1, -0.3, 0.1))
colnames(test_data_cs)


ggplot(data = as.data.frame(test_data_cs),
       aes(x = a_2,
           y = log(y)))+
  geom_point() +
  # facet_wrap(~group, ncol = , scales = "free") +
  theme(legend.position = "bottom") +
  ggtitle("title")

SFA(y ~ a_1 + a_2 + b_1_1 + b_1_2 + b_2_2, 
    cm = ~ z1 + z2,
    # dist = "exp",
    data = as.data.frame(test_data_cs))
lrtest(.Last.value)

# PANEL DATA --------------------------------------------------------------

test_data <- as.data.table(sim_data_panel(k = 30, 
                                          x_coef = c(8, 2, -1),
                                          z_sd = 0.5,
                                          z_intercept = 0.5,
                                          z_coef = c(1, -0.5), 
                                          sigma_v = 1,
                                          sigma_u = 1))

ggplot(data = melt(test_data, id.vars = c("k", "t")),
       aes(x = t,
           y = value,
           color = as.factor(k))) +
  geom_line() +
  facet_wrap(~variable, ncol = , scales = "free") +
  theme(legend.position = "bottom") +
  ggtitle("title")

hist(test_data$eps)

# test_data$
fit <- SFA(y ~ x1 + x2, 
           # dist = "tnorm",
           data = test_data)
print(fit)
lrtest(fit)
ineff_fit <- inefficiencyTerm(fit, "JLMS")
plot(test_data$u ~ ineff_fit)
cor(test_data$u, ineff_fit)
summary(plm(ineff_fit ~ z1 + z2, index = c("k", "t"), data = test_data))



fit_cm <- SFA(y ~ x1 + x2, 
              data = test_data,
              dist = "tnorm",
              grad = "analytic",
              cm = ~ as.factor(k) + z1 + z2,
              optim_control = list(trace = 1, maxit = 1000))
print(fit_cm)
lrtest(fit_cm)
ineff_fit_cm <- inefficiencyTerm(fit_cm, "BC")
plot(test_data$u ~ ineff_fit_cm)
cor(test_data$u, ineff_fit_cm)
cor(ineff_fit, ineff_fit_cm, use = "pairwise.complete.obs")

# library(frontier)
# fit_front <- frontier::sfa(y ~ x1 + x2 | z1 + z2, 
#               data = test_data, maxit = 20)
# summary(fit_front)
# cor(test_data$u, log(1/efficiencies(fit_front)))


# PANEL DATA WITH AR ------------------------------------------------------


test_data <- as.data.table(sim_data_panel_ar(k = 20, 
                                             x_coef = c(8, 2, -1),
                                             z_sd = 0.5,
                                             # z_intercept = 0.5,
                                             z_coef = c(1, -0.5), 
                                             sigma_v = 1, 
                                             sigma_u = 1,
                                             dist = "gamma"))

ggplot(data = melt(test_data, id.vars = c("k", "t")),
       aes(x = t,
           y = value,
           color = as.factor(k))) +
  geom_line() +
  facet_wrap(~variable, ncol = , scales = "free") +
  theme(legend.position = "bottom") +
  ggtitle("title")

hist(test_data$eps)

# test_data$
fit <- SFA(y ~ x1 + x2, 
           data = test_data)
print(fit)
lrtest(fit)
ineff_fit <- inefficiencyTerm(fit, "JLMS")
plot(test_data$u ~ ineff_fit)
cor(test_data$u, ineff_fit)

summary(plm(ineff_fit ~ cumsum(z1) + cumsum(z2), data = test_data, index = c("k", "t")))

fit_cm <- SFA(y ~ x1 + x2, 
              data = test_data,
              dist = "tnorm",
              cm = ~ as.factor(k) + cumsum(z1) + cumsum(z2), 
              optim_control = list(trace = 1))
print(fit_cm)
lrtest(fit_cm)
ineff_fit_cm <- inefficiencyTerm(fit_cm, "JLMS")
plot(test_data$u ~ ineff_fit_cm)
cor(test_data$u, ineff_fit_cm)

