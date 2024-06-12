#fwd sim and some wrangling
library(tidyverse)
library(here)
library(gsl)
source(here("analysis/R/functions.R"))
set.seed(123)

# read in data ---------------------------------------------------------------------------
data <- read.csv(here("analysis/data/raw/fr_pk_spw_har.csv")) |>
  mutate(harvest = round(harvest/1000000, 2),
         spawn = round(spawn/1000000, 2),
         ER = harvest/(harvest+spawn))

AR1.stan.fit <- readRDS(here("analysis/data/generated/AR1-SS-SR.stan.fit.rds"))
AR1.model.pars <- rstan::extract(AR1.stan.fit)

TV.stan.fit <- readRDS(here("analysis/data/generated/TV-alpha-SS-SR.stan.fit.rds"))
TV.model.pars <- rstan::extract(TV.stan.fit)

# calc forecast error ---
for.error <- read.csv(here("analysis/data/raw/PinkSalmonPrediction&ObservationDataFile(from Merran).csv")) |>
  mutate(error = (abs(PreSeasonForecast-FinalRunSize)/FinalRunSize)) #by year

if(FALSE){ #check trends and distribution of error
  ggplot(for.error, aes(YearX, error)) +
    geom_point() +
    stat_smooth(method = "lm") +
    labs(y="forecast error (CV)", x = "year")
  ggplot(lm(error~YearX, data = for.error), aes(x = .fitted, y = .resid)) +
    geom_point() +
    stat_smooth()
}

for.error <- for.error |>
  pull(error) |>
  mean() #mean of annual means (MAPE)

# get benchmarks & pars ------------------------------------------------------------------
bench <- matrix(NA,1000,4,
                dimnames = list(seq(1:1000), c("Sgen","Smsy","Umsy", "Seq")))

for(i in 1:1000){
  r <- sample(seq(1,1000),1,replace=TRUE)
  ln_a <- AR1.model.pars$ln_alpha[r]
  b <- AR1.model.pars$beta[r]

  bench[i,2] <- get_Smsy(ln_a, b) #S_MSY
  bench[i,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[i,2]) #S_gen
  bench[i,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
  bench[i,4] <- ln_a/b #S_eq
}

bench[,2] <- bench[,2]*0.8 #make it 80% Smsy

#get some "posteriors" for plotting later
  #quotes because it's 1000 random not 4000; will generate the same plot...
Smsy.8.post <- bench[,2]
Sgen.post <- bench[,1]

bench.quant <- apply(bench, 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
  t()

mean <- apply(bench,2,mean, na.rm=T) #get means of each

benchmarks <- cbind(bench.quant, mean) |>
  as.data.frame() |>
  relocate('50%', 1)

if(FALSE){ #some benchmarks to pass Sue for salmon scanner
  bench.quant <- t(apply(bench, 2, quantile, probs=c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975), na.rm=T))
  rownames(bench.quant) <- c("80% Smsy","Sgen","Umsy")
  mean <- apply(bench, 2, mean, na.rm=T)
  sd <- apply(bench, 2, sd, na.rm=T)
  bench.sue <- cbind(bench.quant,mean,sd)
  write.csv(bench.sue, "C:/Users/GLASERD/Desktop/Sues.bench.csv")
}

#pull some values to use later and in text
Sgen <- benchmarks[1,1]
Smsy.8 <- benchmarks[2,1]
Umsy <- benchmarks[3,1]
Seq <- benchmarks[4,1]
R.Smsy.8 <- Smsy.8/(1-Umsy) #recruitment @ Smsy to be used as BB USR

percentiles <- quantile(data$spawn, probs=c(0.25, 0.5))
lower.25th.Sp <- percentiles[1]
lower.50th.Sp <- percentiles[2]

alpha.CI <- exp(quantile(AR1.model.pars$ln_alpha, probs = c(.1, .5, .9)))
beta.CI <- quantile(AR1.model.pars$beta, probs = c(.1, .5, .9))
sigma.CI <- quantile(AR1.model.pars$sigma_R_corr, probs = c(.1, .5, .9))
phi.CI <- quantile(AR1.model.pars$phi, probs = c(.1, .5, .9))

par.quants <- rbind(alpha.CI, sigma.CI, beta.CI, phi.CI)

#make big table of bench and pars
par.summary <- as.data.frame(rstan::summary(AR1.stan.fit)$summary) |>
  select(mean, n_eff, Rhat)

#summarise not alphas...
par.summary <- filter(par.summary, row.names(par.summary) %in% c('ln_alpha', 'beta',
                                                                  'phi', 'sigma_R')) |>
  slice(1,4,3,2)

pars <- cbind(par.quants, par.summary)

bench.par.table <- bind_rows(benchmarks, pars) |>
  mutate(n_eff = round(n_eff, 0),
         Rhat = round(Rhat, 4))

bench.par.table[,1:4] <- round(bench.par.table[,1:4], 2)
bench.par.table[is.na(bench.par.table)] <- ""

rownames(bench.par.table) <- c("$S_{gen}$", "80% $S_{MSY}$", "$U_{MSY}$", "$S_{eq}$",
                               "$\\alpha$", "$\\beta$", "$\\phi$", "$\\sigma_{R}$")
colnames(bench.par.table) <- c("Median", "10th percentile", "90th percentile", "Mean",
                               "$N_{eff}$", "$\\hat{R}$")

# initialize the sim ---------------------------------------------------------------------
last.yr <- max(data$year) #final yr of model fit
sim.gens <- 1+5 #final state in model + nyrs (i.e. gens for pinks) to fwd sim
n.sims <- 1000
states <- c("R", "S", "C", "U", "under_Sgen", "over_Smsy") #don't THINK I need to track residual?
last.yr.ind <- ncol(TV.model.pars$lnR) #index the last year of data
scenarios <- c("base", "low_prod")
HCRs <- c("current", "PA_alt", "no_fishing")
OU.CV <- 0.1 #assumed outcome uncertainty, did forecast error earlier

#subset posterior for low productivity
recent.a <- apply(TV.model.pars$ln_alpha[,(last.yr.ind-2):last.yr.ind], 1, median) #median of last 3 gens
recent.a.CIs <- exp(quantile(recent.a, probs = c(.1, .5, .9))) #CIs for reporting
recent.a.50s <- exp(quantile(recent.a, probs = c(.25, .75))) #range to trigger reassessment
recent.a.10th <- quantile(recent.a, .1)

low.a.rows <- which(recent.a<recent.a.10th) #rows to subset later
#overwrite low.a pop dynamics to subset later
low.a.ln.alpha <- recent.a[low.a.rows]
low.a.beta <- TV.model.pars$beta[low.a.rows]
low.a.sigma_R_corr <- TV.model.pars$sigma_R_corr[low.a.rows]

fwd.states <- array(NA, dim = c(length(scenarios), length(HCRs), n.sims, sim.gens, length(states)))

for(i in 1:length(scenarios)){
  for(j in 1:length(HCRs)){
    scenario <- scenarios[i]
    HCR <- HCRs[j]
    for(k in 1:n.sims){
      #take a random correlated slice
      r <- sample(length(TV.model.pars$beta), 1, replace = TRUE) #random draw
      #draw FULLY CORRELATED pars, states, and benchmarks for the sim
      ln_alpha <- recent.a[r]
      beta <- TV.model.pars$beta[r]
      sigma_R_corr <- TV.model.pars$sigma_R_corr[r]
      #estimate draw-specific benchmarks for relative performance measures later
      sub.Smsy <- get_Smsy(ln_alpha, beta)
      sub.Smsy.8 <- sub.Smsy*0.8
      sub.Sgen <- get_Sgen(exp(ln_alpha), beta, -1, 1/beta*2, sub.Smsy)
      #draw final states from model to start fwd sim from
      R <- TV.model.pars$R[r, last.yr.ind]
      S <- TV.model.pars$S[r, last.yr.ind]
      C <- TV.model.pars$C[r, last.yr.ind]
      U <- TV.model.pars$U[r, last.yr.ind]

      fwd.states[i,j,k,1, ] <- c(R,S,C,U,NA,NA) #states and holders for above/below bench

      if(scenario == "low_prod"){#subset the low.a dynamics only
        r.2 <- sample(length(low.a.rows), 1, replace = TRUE)
        ln_alpha <- low.a.ln.alpha[r.2]
        beta <- low.a.beta[r.2]
        sigma_R_corr <- low.a.sigma_R_corr[r.2]
      }
      for(l in 2:sim.gens){
        #forward estimate recruits
        last.S <- fwd.states[i,j,k,l-1,2]
        R <- exp(ln_alpha)*last.S*exp(-beta*last.S+sigma_R_corr) #predict R
        R <- R*rlnorm(1, 0, for.error) #add forecast error
        #then apply the HCR
        if(HCR == "current"){post_HCR <- current_HCR(R, OU=1+rnorm(1, 0, OU.CV))}
        if(HCR == "PA_alt"){post_HCR <- PA_HCR(R, OU=1+rnorm(1, 0, OU.CV),
                                                Sgen=Sgen, R.Smsy=R.Smsy.8, Umsy=Umsy)}
        if(HCR == "no_fishing"){post_HCR <- c(R, 0, 0)}
        under.Sgen <- post_HCR[1] < sub.Sgen
        over.Smsy.8 <- post_HCR[1] > sub.Smsy.8
        fwd.states[i,j,k,l, ] <- c(R, post_HCR, under.Sgen, over.Smsy.8) #write state to array
      }
    }
  }
}

# pull some random individual sims for plotting...
ind.sims <- NULL

for(i in 1:length(scenarios)){
  for(j in 1:length(HCRs)){
    for(k in 1:2){ #can toggle number of sims here
      sub <- fwd.states[i,j,,,]
      r <- sample(n.sims, 1, replace = TRUE)
      single.sim <- sub[r,,] |>
        as.data.frame() |>
        select(2,3) |>
        mutate(year = seq(max(data$year), max(data$year)+((sim.gens-1)*2), by = 2),
               scenario = scenarios[i],
               HCR = HCRs[j],
               sim = as.factor(k))
      ind.sims <- rbind(ind.sims, single.sim)
    }
  }
}
colnames(ind.sims)[1:2] <- c("spawners", "catch")

#wrangle fwd.sim into summary array of [yr, states(+CIs), HCR] for plotting...
fwd.sim <- NULL
yrs <- as.character(seq(from = last.yr,
                        to = last.yr+((sim.gens-1)*2), by=2)) #final yr of sim + fwd sims

for(i in 1:length(scenarios)){
  for(j in 1:length(HCRs)){
    scenario <- scenarios[i]
    HCR <- HCRs[j]
    sub <- fwd.states[i,j,,,]
     for(k in 1:sim.gens){
      sub_sub <- as.data.frame(sub[,k,]) |>
        reframe(across(1:4, quantile_df, .unpack = TRUE)) |> #could use a better fun.?
        mutate(quant = c("lwr", "med", "upr")) |>
        pivot_wider(values_from = 1:4, names_from = quant) |>
        mutate(scenario = scenario,
               HCR = HCR,
               year = as.numeric(yrs[k])) |>
        relocate(c(year, HCR), 1)
      fwd.sim <- rbind(fwd.sim, sub_sub)
    }
  }
}

fwd.sim <- fwd.sim |>
  relocate(scenario, .after=2) |>
  as.data.frame()

colnames(fwd.sim) <- c("year", "HCR", "scenario", "R_lwr", "R", "R_upr", "S_lwr",  "S",
                       "S_upr", "C_lwr", "C", "C_upr", "U_lwr", "U", "U_upr")

# summarise performance metrics------------------------------------------------------
perf.metrics <- NULL

#get relative catch index
rel.catch.index <- filter(data, year >= 2000) |>
  arrange(desc(harvest)) |>
  slice(1:3) |>
  summarise(mean(harvest)) |>
  pull()

for(i in 1:length(scenarios)){
  for(j in 1:length(HCRs)){
  sub.data <- fwd.states[i,j,,2:(sim.gens),] #slice HCR to NOT include final year of observed data.
  Cs <- sub.data[,,3] #catch sims

  #getting total count of points that dipped above/below lines
  below.Sgen <- (length(which(sub.data[,,5]==1))/length(sub.data[,,5]))*100
  above.Smsy.8 <- (length(which(sub.data[,,6]==1))/length(sub.data[,,6]))*100

  #then treating catch more like a distribution and describing the distribution of medians
    #by draw, kind of like describing the distribution of the intercepts...
  catch <- 1:nrow(Cs) |>
    map_dbl(\(x) median(Cs[x,])) |>
    quantile(probs = c(.1, .5, .9)) |>
    round(2)
  catch <- as.character(paste0(catch[2], " (", catch[1], "-", catch[3], ")"))

  catch.stability <- 1:nrow(Cs) |>
    map_dbl(\(x) median(Cs[x,])/sd(Cs[x,])) |>
    quantile(probs = c(.1, .5, .9), na.rm=T) |> #no NAs, but breaks without na.rm=T, maybe because many 0's?
    round(2)

  catch.stability <- as.character(paste0(catch.stability[2], " (", catch.stability[1], "-",
                                          catch.stability[3], ")"))

  catch.index <- length(which(Cs > rel.catch.index))/length(Cs)*100 |> #total points that go above index
    round(2)

  perf.metrics <- rbind(perf.metrics, data.frame(scenario = rep(scenarios[i],5),
                                                 HCR = rep(HCRs[j],5),
                                                 value = c(below.Sgen, above.Smsy.8,
                                                           catch, catch.stability, catch.index),
                                     metric = c("below.Sgen", "above.Smsy.8", "median annual catch",
                                                "catch.stability", "catch index")))
  }
}

perf.metrics <- perf.metrics |>
  pivot_wider(names_from = metric, values_from = value) |>
  as.data.frame() |>
  mutate(catch.stability = ifelse(HCR == "no_fishing", "NA", catch.stability),
         scenario = gsub("_", " ", scenario),
         HCR = gsub("_", " ", HCR))

write.csv(perf.metrics, here("analysis/data/generated/perf-metrics.csv"), row.names = FALSE)

#take the trash out --
rm(beta,ln_a, ln_alpha, C, Cs, catch, catch.stability, fwd.states, bench, bench.quant,
   HCR, HCRs, i,j,k,last.S, last.yr, sub.data, n.sims, post_HCR, r, R, S, sigma_R_corr,
   sim.gens, states, below.Sgen, U, last.yr.ind, above.Smsy.8, over.Smsy.8, catch.index,
   sub.Sgen, under.Sgen, par.quants, par.summary, pars, sub, scenario, scenarios,
   single.sim, sub_sub, b, l, low.a.beta, low.a.ln.alpha, low.a.rows, low.a.sigma_R_corr,
   mean, percentiles, r.2, sub.Smsy.8, TV.stan.fit, AR1.stan.fit, recent.a, sub.Smsy, yrs)
