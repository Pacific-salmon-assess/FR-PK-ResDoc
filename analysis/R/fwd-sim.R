#fwd sim and some wrangling
library(tidyverse)
library(here)
library(gsl)
source(here("analysis/R/functions.R"))
set.seed(123)

# read in data ---------------------------------------------------------------------------
stan.fit <- readRDS(here("analysis/data/generated/SS-SR_AR1.stan.fit.rds"))
model.pars <- rstan::extract(stan.fit)

stan.fit.93 <- readRDS(here("analysis/data/generated/SS-SR_AR1.stan.fit.93.rds"))
model.pars.93 <- rstan::extract(stan.fit.93)

data <- read.csv(here("analysis/data/raw/fr_pk_spw_har.csv")) |>
  mutate(harvest = round(harvest/1000000, 2),
         spawn = round(spawn/1000000, 2),
         ER = harvest/(harvest+spawn))

# get benchmarks & pars ------------------------------------------------------------------
bench <- matrix(NA,1000,4,
                dimnames = list(seq(1:1000), c("Sgen","Smsy","Umsy", "Seq")))

for(i in 1:1000){
  r <- sample(seq(1,1000),1,replace=TRUE)
  ln_a <- model.pars$ln_alpha[r]
  b <- model.pars$beta[r]
  bench[i,2] <- get_Smsy(ln_a,b) #S_MSY
  bench[i,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2,bench[i,2]) #S_gen
  bench[i,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
  bench[i,4] <- ln_a/b #S_eq
}

bench[,2] <- bench[,2]*0.8 #correct to 80% Smsy
bench.quant <- apply(bench, 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
  t()

mean <- apply(bench,2,mean, na.rm=T)

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
low.ln.a <- quantile(model.pars$ln_alpha[which(model.pars$ln_alpha <= quantile(model.pars$ln_alpha, probs = 0.1))],.5)[]
alpha.CI <- quantile(exp(model.pars$ln_alpha), probs = c(.1, .5, .9))
beta.CI <- quantile(model.pars$beta, probs = c(.1, .5, .9))
phi.CI <- quantile(model.pars$phi, probs = c(.1, .5, .9))
sigma.CI <- quantile(model.pars$sigma_R_corr, probs = c(.1, .5, .9))

par.quants <- rbind(alpha.CI, beta.CI, phi.CI, sigma.CI)

#make big table of bench and pars
par.summary <- as.data.frame(rstan::summary(stan.fit)$summary) |>
  select(mean, n_eff, Rhat)

par.summary <- filter(par.summary, row.names(par.summary) %in% c('ln_alpha', 'beta', 'phi', 'sigma_R'))
par.summary[1,1] <- exp(par.summary[1,1])

pars <- cbind(par.quants, par.summary)

bench.par.table <- bind_rows(benchmarks, pars) |>
  mutate(n_eff = round(n_eff, 0),
         Rhat = round(Rhat, 4))

bench.par.table[,1:4] <- round(bench.par.table[,1:4], 2)
bench.par.table[is.na(bench.par.table)] <- ""


rownames(bench.par.table) <- c("$S_{gen}$", "80% $S_{MSY}$", "$U_{MSY}$", "$S_{eq}$",
                               "$\\alpha$", "$\\beta$", "$\\phi$", "$\\sigma$")
colnames(bench.par.table) <- c("Median", "10th percentile", "90th percentile", "Mean",
                               "$N_{eff}$", "$\\hat{R}$")

# initialize the sim ---------------------------------------------------------------------
last.yr <- max(data$year) #final yr of model fit
sim.gens <- 1+5 #final state in model + nyrs (i.e. gens for pinks) to fwd sim
n.sims <- 1000
states <- c("R", "S", "C", "U", "lnresid")
last.yr.ind <- ncol(model.pars$lnR) #index the last year of data
HCRs <- c("current", "PA_alt",
          "low_a_current", "low_a_PA_alt",
          "recent_current", "recent_PA_alt")
OU.CV <- 0.1 #assumed outcome uncertainty

# calc forecast error
  #using mean absolute percent error approach (MAPE)
for.error <- read.csv(here("analysis/data/raw/PinkSalmonPrediction&ObservationDataFile(from Merran).csv")) |>
  mutate(error = (abs(PreSeasonForecast-FinalRunSize)/FinalRunSize)) #divided by 2 to account for in-season adjustment

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
  mean()

# do fwd sim -----------------------------------------------------------------------------
fwd.states <- array(NA, dim = c(n.sims, length(states)+2, sim.gens, length(HCRs)))
ref.pts <- array(NA, dim = c(n.sims, 2, sim.gens-1, length(HCRs)))

for(i in 1:length(HCRs)){
  HCR <- HCRs[i]
  sub.pars <- model.pars #overwrite each time so it doesn't break when subsetting in loop
  if(grepl("low_a", HCR)){ #overwrite posterior & final state to low productivity if "low_a"
    low_a_rows <- which(sub.pars$ln_alpha <= quantile(sub.pars$ln_alpha, probs = 0.1))
    #subset posterior to only draw parms and states from low productivity draws
    sub.pars$ln_alpha <- sub.pars$ln_alpha[low_a_rows]
    sub.pars$beta <- sub.pars$beta[low_a_rows]
    sub.pars$sigma_R_corr <- sub.pars$sigma_R_corr[low_a_rows]
    sub.pars$phi <- sub.pars$phi[low_a_rows]
    #draw sub.pars states to start from
    sub.pars$R <- sub.pars$R[low_a_rows,]
    sub.pars$S <- sub.pars$S[low_a_rows,]
    sub.pars$C <- sub.pars$C[low_a_rows,]
    sub.pars$U <- sub.pars$U[low_a_rows,]
  }
  if(grepl("recent", HCR)){sub.pars <- model.pars.93
  last.yr.ind <- ncol(model.pars.93$lnR)
  }
  for(j in 1:n.sims){
    r <- sample(length(sub.pars$ln_alpha), 1, replace = TRUE) #random draw
    #draw parms for the sim
    ln_alpha <- sub.pars$ln_alpha[r]
    beta <- sub.pars$beta[r]
    sigma_R_corr <- sub.pars$sigma_R_corr[r]
    phi <- sub.pars$phi[r]
    #estimate draw-specific benchmarks for relative performance measures later
    sub.Smsy.8 <- get_Smsy(ln_alpha, beta)*.8
    sub.Sgen <- get_Sgen(exp(ln_alpha), beta, -1, 1/beta*2, sub.Smsy.8)
    #draw final states from model to start fwd sim from
    R <- sub.pars$R[r, last.yr.ind]
    S <- sub.pars$S[r, last.yr.ind]
    C <- sub.pars$C[r, last.yr.ind]
    U <- sub.pars$U[r, last.yr.ind]
    last.lnresid <- sub.pars$lnresid[r, last.yr.ind]
    fwd.states[j,,1,i] <- c(R,S,C,U,last.lnresid, NA, NA) #write final state from model to array
    for(k in 2:sim.gens){
      #get previous state
      last.S <- fwd.states[j,2,k-1,i]
      last.lnresid <- fwd.states[j,5,k-1,i]
      #calc R and pred R
      R <- exp(ln_alpha)*last.S*exp(-beta*last.S+phi*last.lnresid+log(sigma_R_corr))
      R <- R*rlnorm(1, 0, for.error) #add forecast error
      pred.R <- exp(ln_alpha)*last.S*exp(-beta*last.S+phi*last.lnresid)
      last.lnresid <- log(R)-log(pred.R)
      #apply HCR
      if(grepl("current", HCR)){post_HCR <- current_HCR(R, OU=1+rnorm(1, 0, OU.CV))}
      if(grepl("PA", HCR)){post_HCR <- PA_HCR(R, OU=1+rnorm(1, 0, OU.CV),
                                              Sgen=Sgen, R.Smsy=R.Smsy.8, Umsy=Umsy)}
      #make binary obs for some summaries later
      under.Sgen <- post_HCR[1] < sub.Sgen
      over.Smsy.8 <- post_HCR[1] > sub.Smsy.8

      fwd.states[j,,k,i] <- c(R,post_HCR,last.lnresid, under.Sgen, over.Smsy.8) #write it
    }
  }
}
#t(fwd.states[j,,,i]) #see sim year(rows) where it potentially breaks

# pull some individual sims for plotting...
ind.sims <- NULL

for(i in 1:2){
  sub <- fwd.states[,,,i] #subset 1st or 2nd HCR
  for(j in 1:2){
    k <- sample(n.sims, 1, replace = TRUE)
    single.sim <- t(sub[k,,]) |>
      as.data.frame() |>
      select(2,3) |>
      mutate(year = seq(max(data$year), max(data$year)+((sim.gens-1)*2), by = 2),
             HCR = HCRs[i],
             sim = as.factor(j))
    ind.sims <- rbind(ind.sims, single.sim)
  }
}
colnames(ind.sims)[1:2] <- c("spawners", "catch")

#wrangle fwd.sim into summary array of [yr, states(+CIs), HCR] for plotting...
fwd.sim <- NULL
yrs <- as.character(seq(from = last.yr,
                        to = last.yr+((sim.gens-1)*2), by=2)) #final yr of sim + fwd sims

for(i in 1:length(HCRs)){
  sub <- fwd.states[,,,i]
  HCR <- HCRs[i]
  for(j in 1:sim.gens){
    sub_sub <- as.data.frame(sub[,,j]) |>
      reframe(across(1:5, quantile_df, .unpack = TRUE)) |> #could use a better fun.?
      mutate(quant = c("lwr", "med", "upr")) |>
      pivot_wider(values_from = 1:5, names_from = quant) |>
      mutate(HCR = HCRs[i],
             year = as.numeric(yrs[j])) |>
      relocate(c(year, HCR), 1)
    fwd.sim <- rbind(fwd.sim, sub_sub)
  }
}

fwd.sim <- fwd.sim |>
  mutate(scenario = case_when(grepl("low", HCR) ~ "low productivity",
                              grepl("recent", HCR) ~ "recent, baseline",
                              TRUE ~ "baseline"),
         HCR = gsub("low_a_|recent_", "", HCR)) |>
  relocate(scenario, .after=2)

colnames(fwd.sim) <- c("year", "HCR", "scenario", "R_lwr", "R", "R_upr", "S_lwr",  "S",
                       "S_upr", "C_lwr", "C", "C_upr", "U_lwr", "U", "U_upr", "lnresid_lwr",
                       "lnresid", "lnresid_upr")

# summarise performance metrics------------------------------------------------------
perf.metrics <- NULL

#get relative catch index
rel.catch.index <- filter(data, year >= 2000) |>
  arrange(desc(harvest)) |>
  slice(1:3) |>
  summarise(mean(harvest)) |>
  pull()

for(i in 1:length(HCRs)){
  sub.data <- fwd.states[,,2:(sim.gens),i] #slice HCR to NOT include final year of observed data.
  Cs <- sub.data[,3,] #catch sims

  #getting total count of points that dipped above/below lines
  below.Sgen <- (length(which(sub.data[,6,]==1))/length(sub.data[,6,]))*100
  above.Smsy.8 <- (length(which(sub.data[,7,]==1))/length(sub.data[,7,]))*100

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

  perf.metrics <- rbind(perf.metrics, data.frame(HCR = rep(HCRs[i],5),
                                                 value = c(below.Sgen, above.Smsy.8,
                                                           catch, catch.stability, catch.index),
                                     metric = c("below.Sgen", "above.Smsy.8", "median annual catch",
                                                "catch.stability", "catch index")))
}

perf.metrics <- perf.metrics |>
  mutate(scenario = case_when(grepl("low", HCR) ~ "low productivity",
                              grepl("recent", HCR) ~ "recent baseline",
                              TRUE ~ "baseline"),
         HCR = gsub("low_a_|recent_", "", HCR)) |>
  mutate(HCR = ifelse(HCR == "PA_alt", "PA alternate", HCR)) |>
  pivot_wider(names_from = metric, values_from = value)

rm(beta,ln_a, ln_alpha, C, Cs, catch, catch.stability, fwd.states, bench, bench.quant,
   HCR, HCRs, i,j,k,last.lnresid,last.S, last.yr, sub.data, low_a_rows, n.sims,
   phi,post_HCR, pred.R, r, R, S, sigma_R_corr, sim.gens, states,sub_sub, below.Sgen,
   ref.pts, sub, sub.pars,yrs, U, last.yr.ind, above.Smsy.8, over.Smsy.8,
   sub.Sgen, under.Sgen, par.quants, par.summary, pars, single.sim)
