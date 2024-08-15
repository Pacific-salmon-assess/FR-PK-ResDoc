library(tidyverse)
library(rstan)
library(here)
library(bayesplot)
set.seed(123)
# load data ------------------------------------------------------------------------------
data <- read.csv(here("analysis/data/raw/fr_pk_spw_har.csv"))

stan.data <- list("T" = length(data$year),
                  "S_obs" = data$spawn/1000000,
                  "H_obs" = data$harvest/1000000,
                  "S_cv" = data$spawn_cv,
                  "H_cv" = data$harvest_cv,
                  "pSmax_mean" = (max(data$spawn)/1000000)*.75, #Smax prior - 75% of max spawners
                  "pSmax_sig" = (max(data$spawn)/1000000)*.75)

# ----------------------------------------------------------------------------------------
# fit model AR1 model used to calc benchmarks, estimate states, etc. ---------------------
# ----------------------------------------------------------------------------------------
AR1.stan.fit <- rstan::stan(file = here("analysis/Stan/ss-sr-ar1.stan"),
                           model_name = "SS-SR-AR1",
                           data = stan.data,
                           chains = 4,
                           iter = 2000,
                           seed = 1,
                           thin = 1,
                           control = list(adapt_delta = 0.99, max_treedepth = 20))

#shinystan::launch_shinystan(AR1.stan.fit) #there if ya want it...
#saveRDS(AR1.stan.fit, file=here("analysis/data/generated/AR1-SS-SR.stan.fit.rds"))

# basic diagnostics ----------------------------------------------------------------------
AR1.model.summary <- as.data.frame(rstan::summary(AR1.stan.fit)$summary)
model.pars.AR1 <- rstan::extract(AR1.stan.fit)

#plot PPC
R <- (data$harvest+data$spawn)/1000000 #recruit data
R_rep <- model.pars.AR1$H_rep[1:500,] + model.pars.AR1$S_rep[1:500,]

ppc_dens_overlay(R, R_rep) +
  xlim(0, 70) +
  theme(legend.position = "none") +
  labs(y = "density", x = "y_est")


# Ideally n_eff for individual parameters are > 0.1 (i.e., 10%) of the iter
# values at zero can be ignored as these are unsampled parameters.
min(AR1.model.summary$n_eff, na.rm = TRUE)/2000
ggplot(AR1.model.summary, aes(n_eff/2000)) +
  geom_histogram() +
  labs(y = "frequency",
       x = "ESS/iter") +
  theme_minimal()

# If chains have not mixed well (i.e., the between- and within-chain estimates don't agree),
# R-hat is > 1. Only use the sample if R-hat is less than 1.05.
ggplot(AR1.model.summary, aes(Rhat)) +
  geom_histogram() +
  labs(y = "frequency",
       x =  expression(hat("R"))) +
  theme_minimal()

max(AR1.model.summary$Rhat, na.rm = T)

# check the chains directly
mcmc_combo(AR1.stan.fit, pars = c("beta", "ln_alpha", "sigma_R_corr", "phi"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())

# how do correlations in lnalpha and beta posteriors look?
pairs(AR1.stan.fit, pars = c("beta", "ln_alpha", "sigma_R_corr", "phi"))


# ----------------------------------------------------------------------------------------
# fit time varying alpha model used to inform biological submodel in fwd sims ------------
# ----------------------------------------------------------------------------------------
TV.stan.fit <- rstan::stan(file = here("analysis/Stan/ss-sr-tv-alpha.stan"),
                 model_name = "SS-SR-TV-alpha",
                 data = stan.data,
                 chains = 4,
                 iter = 2000,
                 seed = 1,
                 thin = 1,
                 control = list(adapt_delta = 0.99, max_treedepth = 20))

saveRDS(TV.stan.fit, file=here("analysis/data/generated/TV-alpha-SS-SR.stan.fit.rds"))

# basic diagnostics ----------------------------------------------------------------------
TV.model.summary <- as.data.frame(rstan::summary(TV.stan.fit)$summary)

# Ideally n_eff for individual parameters are > 0.1 (i.e., 10%) of the iter
  # values at zero can be ignored as these are unsampled parameters.
min(TV.model.summary$n_eff, na.rm = TRUE)/2000
ggplot(TV.model.summary, aes(n_eff/2000)) +
  geom_histogram() +
  labs(y = "frequency",
       x = "ESS/iter") +
  theme_minimal()

# If chains have not mixed well (i.e., the between- and within-chain estimates don't agree),
  # R-hat is > 1. Only use the sample if R-hat is less than 1.05.
ggplot(TV.model.summary, aes(Rhat)) +
  geom_histogram() +
  labs(y = "frequency",
       x =  expression(hat("R"))) +
  theme_minimal()

max(TV.model.summary$Rhat, na.rm = T)

# check the chains directly
mcmc_combo(TV.stan.fit, pars = c("beta", "ln_alpha0", "sigma_R", "sigma_alpha"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())
#want to see alpha but theres 33 of them...

# how do correlations in lnalpha and beta posteriors look?
pairs(TV.stan.fit, pars = c("beta", "ln_alpha0",  "sigma_R", "sigma_alpha"))
