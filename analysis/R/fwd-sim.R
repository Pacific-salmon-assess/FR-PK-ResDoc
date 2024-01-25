#fwd sim and some wrangling
library(tidyverse)
library(here)
library(gsl)
source(here("analysis/R/functions.R"))
set.seed(123)

# read in data ---------------------------------------------------------------------------
stan.fit <- readRDS(here("analysis/data/generated/SS-SR_AR1.stan.fit.rds"))
model.pars <- rstan::extract(stan.fit)

# get the last yr of rec & parameter posteriors -------------------------------------
last.yr <- 2021 #final yr of model fit
sim.gens <- 1+5 #final state in model + nyrs (i.e. gens for pinks) to fwd sim
n.sims <- 1000
states <- c("R", "S", "C", "U", "lnresid")
last.yr.ind <- ncol(model.pars$lnR) #index the last year of data
HCRs <- c("current", "proposed", "low_a_current", "low_a_proposed")
OU.CV <- 0.1 #outcome uncertainty

# calc forecast error---------------------------------------------------------------------
  #using mean absolute percent error approach (MAPE)
for.error <- read.csv(here("analysis/data/raw/PinkSalmonPrediction&ObservationDataFile(from Merran).csv")) |>
  mutate(error = abs(PreSeasonForecast - FinalRunSize)/FinalRunSize)

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

fwd.states <- array(NA, dim = c(n.sims, length(states), sim.gens, length(HCRs)))
ref.pts <- array(NA, dim = c(n.sims, 2, sim.gens-1, length(HCRs)))

if(FALSE){ #turn off so it doesn't render in doc
  #prior checks
  median(model.pars$lnalpha_c) - median(model.pars$lnalpha)
  hist(model.pars$lnalpha)
  hist(model.pars$lnalpha_c)
  hist(model.pars$sigma_R_corr)
}

for(i in 1:length(HCRs)){
  HCR <- HCRs[i]
  sub.pars <- model.pars
  if(grepl("low_a", HCR)){ #overwrite posterior & final state to low productivity if "low_a"
    low_a_rows <- which(sub.pars$lnalpha_c <= quantile(sub.pars$lnalpha_c, probs = 0.1))
    #subset posterior to only draw parms and states from low productivity draws
    sub.pars$lnalpha_c <- sub.pars$lnalpha_c[low_a_rows]
    sub.pars$beta <- sub.pars$beta[low_a_rows]
    sub.pars$sigma_R_corr <- sub.pars$sigma_R_corr[low_a_rows]
    model.pars$phi <- sub.pars$phi[low_a_rows]
    #draw sub.pars states to start from
    sub.pars$R <- sub.pars$R[low_a_rows,]
    sub.pars$S <- sub.pars$S[low_a_rows,]
    sub.pars$C <- sub.pars$C[low_a_rows,]
    sub.pars$U <- sub.pars$U[low_a_rows,]
  }
  for(j in 1:n.sims){
    r <- sample(length(sub.pars$lnalpha_c), 1, replace = TRUE)
    #draw parms for the sim
    lnalpha_c <- sub.pars$lnalpha_c[r]
    beta <- sub.pars$beta[r]
    sigma_R_corr <- sub.pars$sigma_R_corr[r]
    phi <- sub.pars$phi[r]
    #draw final states from model to start fwd sim from
    R <- sub.pars$R[r, last.yr.ind]
    S <- sub.pars$S[r, last.yr.ind]
    C <- sub.pars$C[r, last.yr.ind]
    U <- sub.pars$U[r, last.yr.ind]
    last.lnresid <- sub.pars$lnresid[r, last.yr.ind]
    fwd.states[j,,1,i] <- c(R,S,C,U,last.lnresid) #write final state from model to array
    for(k in 2:sim.gens){
      #get previous state
      last.S <- fwd.states[j,2,k-1,i] # COULD just call 'S' instead of renaming
      last.lnresid <- fwd.states[j,5,k-1,i]
      #calc R and pred R
      R <- exp(lnalpha_c)*last.S*exp(-beta*last.S+phi*last.lnresid+log(sigma_R_corr))
      R <- R*rlnorm(1, 0, for.error) #add forecast error
      pred.R <- exp(lnalpha_c)*last.S*exp(-beta*last.S+phi*last.lnresid)
      last.lnresid <- log(R)-log(pred.R)
      #apply HCR
      if(grepl("current", HCR)){
        post_HCR <- current_HCR(R, OU = 1+rnorm(1, 0, OU.CV))
      }else{
        post_HCR <- alt_HCR(R, OU = 1+rnorm(1, 0, OU.CV))
      }
      fwd.states[j,,k,i] <- c(R,post_HCR,last.lnresid) #populate the next year based on the HCR

      Smsy <- get_Smsy(lnalpha_c, beta)
      Sgen <- get_Sgen(exp(lnalpha_c), beta, -1, 1/beta*2, Smsy)
      ref.pts[j, ,k-1,i] <- c(post_HCR[1]/Smsy, post_HCR[1]/Sgen)
    }
  }
}

#wrangle fwd.sim into summary array of [yr, states(+CIs), HCR]
fwd.states.summary <- NULL
yrs <- as.character(seq(from = last.yr,
                        to = last.yr+((sim.gens-1)*2), by=2)) #final yr of sim + fwd sims

for(i in 1:length(HCRs)){
  sub <- fwd.states[,,,i]
  HCR <- HCRs[i]
  for(j in 1:sim.gens){
    sub_sub <- as.data.frame(sub[,,j]) |>
      reframe(across(1:5, quantile_df, .unpack = TRUE)) |> #could use a better fun.?
      mutate(quant = c("lwr", "lwr_med", "med", "upr_med", "upr")) |>
      pivot_wider(values_from = 1:5, names_from = quant)
    fwd.states.summary <- rbind(fwd.states.summary, sub_sub)
  }
}

fwd.states.summary <- fwd.states.summary |>
  mutate(sim = c(rep(HCRs, each = sim.gens)),
         year = rep(yrs, length(HCRs))) |>
  relocate(c(year, sim), 1) |>
  arrange(year, sim)

colnames(fwd.states.summary) <- c("year", "sim", "R_lwr", "R_mid_lwr", "R", "R_mid_upr",
                                  "R_upr", "S_lwr", "S_mid_lwr", "S", "S_mid_upr", "S_upr",
                                  "C_lwr", "C_mid_lwr", "C", "C_mid_upr", "C_upr", "U_lwr",
                                  "U_mid_lwr", "U", "U_mid_upr", "U_upr", "lnresid_lwr",
                                  "lnresid_mid_lwr", "lnresid", "lnresid_mid_upr", "lnresid_upr")

# summarise some performance metrics------------------------------------------------------
perf.summary <- NULL
perf.total <- NULL
OCP.total <- NULL
for(i in 1:length(HCRs)){
  sub.refs <- ref.pts[,,,i]
  sub.data <- fwd.states[,,,i]
  below.Smsy.all <- nrow(filter(as.data.frame(sub.refs[,1,]), V1<1&V2<1&V3<1))/nrow(sub.refs[,1,])*100
  below.Sgen.all <- nrow(filter(as.data.frame(sub.refs[,2,]), V1<1&V2<1&V3<1))/nrow(sub.refs[,2,])*100
  HCR <- HCRs[i]

  perf.total <- rbind(perf.total, data.frame(sim = HCR,
                                             below.Smsy.all = below.Smsy.all,
                                             below.Sgen.all = below.Sgen.all))

  lwr.OCP <- round(length(which(sub.data[,1,2:4] <= 7.059))/length(sub.data[,1,2:4]), 4)*100 #all yrs of R (except start)
  mid.OCP <- round(length(which(sub.data[,1,2:4] >= 7.059 & sub.data[,1,2:4] <= 20))/
                     length(sub.data[,1,2:4]), 4)*100
  upr.OCP <- round(length(which(sub.data[,1,2:4] >= 20))/length(sub.data[,1,2:4]), 4)*100

  OCP.total <- rbind(OCP.total, data.frame(sim = HCR, lwr.OCP, mid.OCP, upr.OCP))

  for(j in 1:(sim.gens-1)){
    below.Smsy <- (length(which(sub.refs[,1,j]<1))/nrow(sub.refs))*100
    below.Sgen <- (length(which(sub.refs[,2,j]<1))/nrow(sub.refs))*100

    perf.summary <- rbind(perf.summary,
                          data.frame(sim = HCR, year = yrs[j+1], below.Smsy,below.Sgen))
  }
}

#add catch to performance summaries
#could add stability too but need to think about how to do by year
perf.summary <- fwd.states.summary |>
  filter(year != last.yr) |>
  select(year, sim, C) |>
  mutate(catch = round(C, 2)) |>
  select(-C) |>
  left_join(perf.summary, by = c("year", "sim")) |>
  arrange(sim, year)
colnames(perf.summary) <- c("Year", "Simulaiton", "Median catch (millions)", "Below Smsy (%)",
                            "Below Sgen (%)")

perf.total <- fwd.states.summary |>
  filter(year != last.yr) |>
  group_by(sim) |>
  summarise(total.catch = round(sum(C), 2)) |>
  left_join(perf.total, by = "sim") #could add catch stability
colnames(perf.total) <- c("Simulaiton", "Total median catch (millions)", "All yrs. below Smsy (%)",
                          "All yrs. below Sgen (%)")
colnames(OCP.total) <- c("Simulation", "% below lower OCP", "% between OCPs",
                         "% above upper OCP")

rm(beta, C, fwd.states, HCR, HCRs, i,j,k,last.lnresid,last.S, last.yr, lnalpha_c, lwr.OCP,
   mid.OCP, upr.OCP, sub.data, sub.refs, low_a_rows, n.sims, phi, post_HCR, pred.R, r, R,
   S, sigma_R_corr, sim.gens, states,sub_sub, below.Sgen, below.Sgen.all, below.Smsy,
   below.Smsy.all, ref.pts, Sgen,Smsy, sub, sub.pars, yrs, U, last.yr.ind)
