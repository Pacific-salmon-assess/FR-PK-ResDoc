#checking only using data after '93 (a cutoff that rapid status evaluation has used)
library(tidyverse)
library(here)
library(gsl)
library(rstan)
source(here("analysis/R/functions.R"))
source(here("analysis/R/fwd-sim.R")) #mostly to get benchmarks

data <- read.csv(here("analysis/data/raw/fr_pk_spw_har.csv")) |>
  filter(year >= 1993)

stan.data <- list("nyrs" = length(data$year),
                  "nRyrs" = length(data$year) ,
                  "S_obs" = data$spawn/1000000,
                  "H_obs" = data$harvest/1000000,
                  "S_cv" = data$spawn_cv,
                  "H_cv" = data$harvest_cv)

# fit model ------------------------------------------------------------------------------
stan.fit.93 <- rstan::stan(file = here("analysis/Stan/ss-sr-ar1.stan"),
                        model_name = "SS-SR_AR1",
                        data = stan.data,
                        chains = 4,
                        iter = 2000,
                        seed = 1,
                        thin = 1,
                        control = list(adapt_delta = 0.99, max_treedepth = 20))

model.pars.93 <- rstan::extract(stan.fit.93)

#get benchmarks---------------------------------------------------------------------------
bench.93 <- matrix(NA,1000,3,
                dimnames = list(seq(1:1000), c("Sgen","Smsy","Umsy")))

for(i in 1:1000){
  r <- sample(seq(1,1000),1,replace=TRUE)
  a <- model.pars.93$lnalpha[r]
  b <- model.pars.93$beta[r]
  bench.93[i,2] <- get_Smsy(a,b) #Smsy
  bench.93[i,1] <- get_Sgen(exp(a),b,-1,1/b*2,bench.93[i,2]) #Sgen
  bench.93[i,3] <- (1 - lambert_W0(exp(1 - a))) #Umsy
}

bench.93[,2] <- bench.93[,2]*0.8 #correct to 80% Smsy
bench.quant.93 <- apply(bench.93, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)

percentiles.93 <- quantile(data$spawn, probs=c(0.25, 0.5))/1000000

benchmarks.93 <- matrix(NA,5,3)
benchmarks.93[1,] <- c(bench.quant.93[2,2],bench.quant.93[1,2],bench.quant.93[3,2])
benchmarks.93[2,] <- c(bench.quant.93[2,1],bench.quant.93[1,1],bench.quant.93[3,1])
benchmarks.93[3,] <- c(bench.quant.93[2,3],bench.quant.93[1,3],bench.quant.93[3,3])
benchmarks.93[4,1] <- percentiles.93[1]
benchmarks.93[5,1] <- percentiles.93[2]

rownames(benchmarks.93) <- c("80% Smsy","Sgen","Umsy","25th percentile (spawners)",
                          "50th percentile (spawners)")
colnames(benchmarks.93) <- c("median","lower 95% CI","upper 95% CI")

benchmarks <- benchmarks |>
  as.data.frame() |>
  mutate(par = rownames(benchmarks),
         data = "full")

both.benchmarks <- benchmarks.93 |>
  as.data.frame() |>
  mutate(par = rownames(benchmarks),
         data = "post '93") |>
  bind_rows(benchmarks) |>
  arrange(par, data) |>
  mutate(par = case_when(par == '25th percentile (spawners)' ~ '25% Sp.',
                         par == '50th percentile (spawners)' ~ '50% Sp.',
                         TRUE ~ par))

rownames(both.benchmarks) <- NULL

ggplot(both.benchmarks, aes(x = par, y = median)) +
  geom_pointrange(aes(ymin = `lower 95% CI`, ymax = `upper 95% CI`, color = data),
                  position = "jitter") +
  labs(title = "Comparing estimates from the full data vs. post '93",
       x = "benchmark", y = "value")
