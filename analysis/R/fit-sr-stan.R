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

# fit model ------------------------------------------------------------------------------
stan.fit <- rstan::stan(file = here("analysis/Stan/ss-sr-tv-alpha.stan"),
                 model_name = "SS-SR-TV-alpha",
                 data = stan.data,
                 chains = 4,
                 iter = 2000,
                 seed = 1,
                 thin = 1,
                 control = list(adapt_delta = 0.99, max_treedepth = 20))

#shinystan::launch_shinystan(stan.fit)
saveRDS(stan.fit, file=here("analysis/data/generated/SS-SR-TV-alpha.stan.fit.rds"))

# basic diagnostics ----------------------------------------------------------------------
model.summary <- as.data.frame(rstan::summary(stan.fit)$summary)

# Ideally n_eff for individual parameters are > 0.1 (i.e., 10%) of the iter
  # values at zero can be ignored as these are unsampled parameters.
min(model.summary$n_eff, na.rm = TRUE)/2000
ggplot(model.summary, aes(n_eff/2000)) +
  geom_histogram() +
  labs(y = "frequency",
       x = "ESS/iter") +
  theme_minimal()

# If chains have not mixed well (i.e., the between- and within-chain estimates don't agree),
  # R-hat is > 1. Only use the sample if R-hat is less than 1.05.
ggplot(model.summary, aes(Rhat)) +
  geom_histogram() +
  labs(y = "frequency",
       x =  expression(hat("R"))) +
  theme_minimal()

max(model.summary$Rhat, na.rm = T)

# check the chains directly
mcmc_combo(stan.fit, pars = c("beta", "ln_alpha0", "sigma_R", "sigma_alpha"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())
#want to see alpha but theres 33 of them...

# how do correlations in lnalpha and beta posteriors look?
pairs(stan.fit, pars = c("beta", "ln_alpha0",  "sigma_R", "sigma_alpha"))
