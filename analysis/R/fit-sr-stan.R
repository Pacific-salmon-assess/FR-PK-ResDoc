library(tidyverse)
library(rstan)
library(here)
library(bayesplot)
set.seed(123)
# load data ------------------------------------------------------------------------------
data <- read.csv(here("analysis/data/raw/fr_pk_spw_har.csv"))


stan.data <- list("nyrs" = length(data$year),
                  "nRyrs" = length(data$year) ,
                  "S_obs" = data$spawn/1000000,
                  "H_obs" = data$harvest/1000000,
                  "S_cv" = data$spawn_cv,
                  "H_cv" = data$harvest_cv,
                  "pSmax_mean" = (max(data$spawn)/1000000)*.75, #Smax prior - 75% of max spawners
                  "pSmax_sig" = (max(data$spawn)/1000000)*.75)

# fit model ------------------------------------------------------------------------------
stan.fit <- rstan::stan(file = here("analysis/Stan/ss-sr-ar1.stan"),
                 model_name = "SS-SR_AR1",
                 data = stan.data,
                 chains = 4,
                 iter = 2000,
                 seed = 1,
                 thin = 1,
                 control = list(adapt_delta = 0.99, max_treedepth = 20))

#shinystan::launch_shinystan(stan.fit)
#saveRDS(stan.fit, file=here("analysis/data/generated/SS-SR_AR1.stan.fit.rds")) #uncomment when good

# basic diagnostics ----------------------------------------------------------------------
model.summary <- as.data.frame(rstan::summary(stan.fit)$summary)

# Ideally n_eff for individual parameters are > 0.1 (i.e., 10%) of the iter
  # values at zero can be ignored as these are unsampled parameters.
min(model.summary$n_eff, na.rm = TRUE)
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
mcmc_combo(stan.fit, pars = c("ln_alpha", "ln_beta", "beta", "sigma_R", "phi", "lnresid_0"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())

mcmc_combo(stan.fit, pars = c("ln_alpha", "ln_alpha_c"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())

# how do correlations in lnalpha and beta posteriors look?
pairs(stan.fit, pars = c("ln_alpha", "ln_beta", "sigma_R", "phi"))

model.pars <- rstan::extract(stan.fit) #extract pars for later

#-----------------------------------------------------------------------------------------
# Now also fit it to shortened recent data------------------------------------------------
data.93 <- filter(data, year >= 1993)

stan.data.93 <- list("nyrs" = length(data.93$year),
                  "nRyrs" = length(data.93$year) ,
                  "S_obs" = data.93$spawn/1000000,
                  "H_obs" = data.93$harvest/1000000,
                  "S_cv" = data.93$spawn_cv,
                  "H_cv" = data.93$harvest_cv,
                  "pSmax_mean" = (max(data.93$spawn)/1000000)*.75,
                  "pSmax_sig" = ((max(data.93$spawn)/1000000)*.75))

# fit model ------------------------------------------------------------------------------
stan.fit.93 <- rstan::stan(file = here("analysis/Stan/ss-sr-ar1.stan"),
                           model_name = "SS-SR_AR1",
                           data = stan.data.93,
                           chains = 4,
                           iter = 2000,
                           seed = 1,
                           thin = 1,
                           control = list(adapt_delta = 0.99, max_treedepth = 20))
#shinystan::launch_shinystan(stan.fit)
#saveRDS(stan.fit.93, file=here("analysis/data/generated/SS-SR_AR1.stan.fit.93.rds")) #uncomment when good

# basic diagnostics ----------------------------------------------------------------------
model.summary.93 <- as.data.frame(rstan::summary(stan.fit.93)$summary)

# Ideally n_eff for individual parameters are > 0.1 (i.e., 10%) of the iter
# values at zero can be ignored as these are unsampled parameters.
min(model.summary.93$n_eff, na.rm = TRUE)
ggplot(model.summary.93, aes(n_eff/2000)) +
  geom_histogram() +
  labs(y = "frequency",
       x = "ESS/iter") +
  theme_minimal()

# If chains have not mixed well (i.e., the between- and within-chain estimates don't agree),
# R-hat is > 1. Only use the sample if R-hat is less than 1.05.
ggplot(model.summary.93, aes(Rhat)) +
  geom_histogram() +
  labs(y = "frequency",
       x =  expression(hat("R"))) +
  theme_minimal()

max(model.summary.93$Rhat, na.rm = T)

# check the chains directly
mcmc_combo(stan.fit.93, pars = c("ln_alpha", "ln_beta", "sigma_R", "phi", "lnresid_0"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())
mcmc_combo(stan.fit.93, pars = c("ln_alpha", "ln_alpha_c"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())
# how do correlations in lnalpha and beta posteriors look?
pairs(stan.fit.93, pars = c("ln_alpha", "ln_beta", "sigma_R", "phi"))

model.pars.93 <- rstan::extract(stan.fit.93)

#checking for extreme ln_alpha_c values that will break fwd sim.
head(sort(model.pars.93$ln_alpha_c, decreasing = TRUE))
head(sort(model.pars$ln_alpha_c, decreasing = TRUE))
