library(tidyverse)
library(here)
library(gsl)
library(cowplot)
library(scales) #for pretty_breaks() on fig axes
source(here("analysis/R/fwd-sim.R"))

# read in data and fit -------------------------------------------------------------------
HCRs <- read.csv(here("analysis/data/raw/HCRs.csv")) |>
  filter(HCR != "alt.TAM.upper")

avg_mass <- read.csv(here("analysis/data/raw/bio/HistoricalPinkWeightDownload.csv"))
colnames(avg_mass) <- c("year", "avg.weight", "reference")

compitetors <- read.csv(here("analysis/data/raw/bio/Ruggerone_Irvine_2018_TS21.csv"))

#WRANGLING -------------------------------------------------------------------------------
#latent states of spawners and recruits---
spwn.quant <- apply(model.pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:32]
rec.quant <- apply(model.pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,2:33]

brood_t <- as.data.frame(cbind(data$year[1:32],t(spwn.quant), t(rec.quant))) |>
  round(2)
colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")

#SR relationship---
spw <- seq(0,max(brood_t$S_upr),length.out=100)
max_samples <- dim(model.pars$ln_alpha)
SR_pred <- matrix(NA,100,max_samples)

for(i in 1:max_samples){
  a <- model.pars$ln_alpha[i] #corrected here too? - it makes it go crazy
  b <- model.pars$beta[i]
  SR_pred[,i] <- (exp(a)*spw*exp(-b*spw))
}

SR_pred <- as.data.frame(cbind(spw,t(apply(SR_pred,c(1),quantile,probs=c(0.1,0.5,0.9),na.rm=T))))|>
  round(2)
colnames(SR_pred) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr")

#create spawner & esc df --- #should map() these...
spwn.quant <- apply(model.pars$S, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
spwn_df <- as.data.frame(cbind(data$year, t(spwn.quant)))
colnames(spwn_df) <- c("year","lwr","mid_lwr","mid","mid_upr","upr")

er.quant <- apply(model.pars$U, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
er_df <- as.data.frame(cbind(data$year, t(er.quant)))
colnames(er_df) <- c("year","lwr","mid_lwr","mid","mid_upr","upr")

R.quant <- apply(model.pars$R, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
R_df <- as.data.frame(cbind(data$year, t(R.quant)))
colnames(R_df) <- c("year","lwr","mid_lwr","mid","mid_upr","upr")

C.quant <- apply(model.pars$C, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
C_df <- as.data.frame(cbind(data$year, t(C.quant)))
colnames(C_df) <- c("year","lwr","mid_lwr","mid","mid_upr","upr")

kobe_df <- data.frame(S = spwn_df$mid,
                      S_LCI = spwn_df$lwr,
                      S_UCI = spwn_df$upr,
                      U = er_df$mid,
                      U_LCI = er_df$lwr,
                      U_UCI = er_df$upr,
                      year = spwn_df$year) |>
  #index benchmarks table
  mutate(U_Umsy = round(U/Umsy ,2),
         U_Umsy_LCI = round(U_LCI/Umsy ,2),
         U_Umsy_UCI = round(U_UCI/Umsy ,2),
         S_Smsy = round(S/(Smsy.8*1.25), 2), #correct from 80%Smsy to Smsy
         S_Smsy_LCI = round(S_LCI/(Smsy.8*1.25), 2),
         S_Smsy_UCI = round(S_UCI/(Smsy.8*1.25), 2))

rm(er.quant, rec.quant, spwn.quant, max_samples, percentiles, a,b, spw)

# INFERENCE ------------------------------------------------------------------------------
#get DATA (not latent state) of avg run sizes for context/background section
rec_prop <- data |>
  mutate(rec = harvest + spawn,
         ER = (harvest/rec)*100)

# FIGS & TABLES --------------------------------------------------------------------------
theme_set(theme_bw(base_size = 14))

#wrapper function to set dims before saving, so things scale proper each time
  #ASSUMED this is the correct scaling based on index.Rmd
my.ggsave <- function(filename = default_name(plot), width= 9, height = 5.562, dpi= 180){
  ggsave(filename=filename, width=width, height=height, dpi=dpi)
}

# plot catch & escapement ---
catch_esc <- data |>
  select(year, harvest, spawn) |>
  pivot_longer(!year, names_to = "type", values_to = "n")

p1 <- ggplot(catch_esc, aes(x = year, y = n)) +
  geom_bar(position="stack", stat="identity", aes(fill = type)) +
  geom_line(data = data, aes(x=year, y = ER*30), color = "red", lwd = 1) +   #ADD ER to plot
  scale_fill_manual(values = c("darkgrey", "black"), name = "Return type") +
  scale_y_continuous(sec.axis = sec_axis(~(./30)*100, name = "Exploitation Rate (%)")) +
  theme(legend.position = "bottom") +
  labs(x = "Return year",
       y = "Total return (millions of fish)")

p1
my.ggsave(here("figure/catch-esc.png"))

# HCR and realized harvest ---
p2 <- ggplot(filter(HCRs, HCR=="current")) +
  geom_line(aes(x=run_size, y = ER)) +
  geom_point(data = filter(data, year >= 1987),
             aes(x=(harvest+spawn), y = harvest/(harvest+spawn), color = year)) +
  scale_color_viridis_c() +
  labs(x = "run size", y = "target exploitation rate") +
  theme(legend.position = "bottom")
p <- plot_grid(p1, p2, ncol = 1)
p
my.ggsave(here("figure/catch-esc-HCR.png"))

#plot avg mass through time---
if(FALSE){
  #can't figure out how to scale these using line and bar graph.
  #Can't set limits easily with 2 axes. I get we want to spread the line but not sure how.
  rescale <- function(x) (x-0)/(max(x) - 0) * max(compitetors$n_pink)
  rescale(avg_mass$avg.weight)

  ggplot(avg_mass, aes(year, rescale(avg.weight))) +
    geom_col(data = compitetors, aes(x = return_yr, y = n_pink), fill = "pink", color = "pink") +
    geom_line() +
    geom_point() +
    #scale_y_continuous(sec.axis = sec_axis(~.*200, name = "N. Pacific Pink abundance (millions)")) +
    labs(x = "Year", y = "Average body mass (kg)")
}

ggplot(avg_mass, aes(year, avg.weight)) +
  geom_col(data = compitetors, aes(x = return_yr, y = n_pink/200), fill = "pink", color = "pink") +
  geom_line() +
  geom_point() +
  scale_y_continuous(sec.axis = sec_axis(~.*200, name = "N. Pacific Pink abundance (millions)")) +
  labs(x = "Year", y = "Average body mass (kg)")

my.ggsave(here("figure/avg-mass.png"))

# plot HCRs ---
p1 <- ggplot(filter(HCRs, HCR!="alt.TAM.lower"), aes(x=run_size, y=ER, color = HCR)) +
  geom_line(size=1.1, alpha = 0.7) +
  geom_vline(xintercept = R.Smsy.8) +
  annotate("text", x = R.Smsy.8+1.5, y = .7,
           label = expression(italic(R[paste("80%",S)[MSY]]))) +
  geom_vline(xintercept = Sgen) +
  annotate("text", x = Sgen+1, y = .7,
           label = "italic(S[gen])", parse = TRUE) +
  scale_color_viridis_d() +
  ylim(c(0,1)) +
  labs(x = NULL,
       y = "Target ER") +
  guides(color = "none") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- ggplot(filter(HCRs, HCR!="alt.TAM.lower"), aes(x=run_size, y=esc_goal, color = HCR)) +
  geom_line(size=1.1, alpha = 0.7) +
  scale_color_viridis_d() +
  geom_vline(xintercept = R.Smsy.8) +
  geom_vline(xintercept = Sgen) +
  labs(x = "Run size (millions)",
       y = "Target escapement") +
  theme(legend.position = "bottom")

p <- plot_grid(p1, p2, nrow = 2)
p
my.ggsave(here("figure/HCRs.png"))

# plot SR relationship ---
ggplot() +
  geom_ribbon(data = SR_pred, aes(x = Spawn, ymin = Rec_lwr, ymax = Rec_upr),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_line(data = SR_pred, aes(x = Spawn, y = Rec_med), color="black", size = 1) +
  geom_errorbar(data = brood_t, aes(x= S_med, y = R_med, ymin = R_lwr, ymax = R_upr),
                colour="grey", width=0, size=0.3) +
  geom_errorbarh(data = brood_t, aes(y = R_med, xmin = S_lwr, xmax = S_upr),
                 height=0, colour = "grey", size = 0.3) +
  geom_point(data = brood_t,
             aes(x = S_med,
                 y = R_med,
                 color=BroodYear),
             size = 3)+
  coord_cartesian(xlim=c(0, 20), ylim=c(0,max(brood_t[,7])), expand = FALSE) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (millions)",
      y = "Recruits (millions)") +
  theme(legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_text(size=9),
    legend.text = element_text(size=8))+
  geom_abline(intercept = 0, slope = 1,col="dark grey")

my.ggsave(here("figure/SRR.png"))

# plot Kobe ---
ggplot(kobe_df, aes(S_Smsy, U_Umsy)) +
  #draw data and error bars on final year
  geom_point(aes(color = year), size=3) +
  #geom_path() + #if you want to connect the dots
  geom_errorbar(data = filter(kobe_df, year == max(kobe_df$year)),
                aes(x = S_Smsy, ymin = U_Umsy_LCI, ymax = U_Umsy_UCI), width = 0) +
  geom_errorbarh(data = filter(kobe_df, year == max(kobe_df$year)),
                 aes(y = U_Umsy, xmin = S_Smsy_LCI, xmax = S_Smsy_UCI), height = 0) +
  #add "crosshairs"
  geom_vline(xintercept = 1, lty = 2) +
  geom_hline(yintercept = 1, lty = 2) +
  #geom_vline(xintercept = 0.8, lty = 3) +
  #add labels to 80% Smsy, first and last year of data
  #annotate("text", x = 0.8, y = .4, hjust = 1,
  #         label = expression(italic(paste("80%",S)[MSY]))) +
  geom_text(data = filter(kobe_df, year== min(kobe_df$year)|year== max(kobe_df$year)),
            aes(x = S_Smsy, y = U_Umsy, label = c("'59", "'23")), #CHANGE THESE WITH NEW DATA!
            hjust = 0-.2, vjust = 0-.2) +
  scale_colour_viridis_c(name="Year") +
  labs(y="U/Umsy", x= "S/Smsy") +
  theme(legend.position = "bottom")

ggsave(here("figure/kobe.png"), width= 9, height = 9, dpi= 180)

# then residuals---
resid.quant <- apply(model.pars$lnresid, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))[,1:33]

resids <- as.data.frame(cbind(data$year[1:32], t(resid.quant)))
colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr")

ggplot(resids, aes(x=year, y = mid)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
  geom_line(lwd = 1.1) +
  coord_cartesian(ylim=c(-2,2)) +
  labs(x = "Return year",
       y = "Recruitment residuals") +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)

my.ggsave(here("figure/rec-resid.png"))

#plot fwd sims of spawners & catch ---
#how much of the old data do you want to show?
d_start <- 2013
d_end <- 2023

#spawners---
p1 <- ggplot(data = filter(fwd.sim, scenario == "baseline")) +
  #draw the fwd.sim
  geom_ribbon(aes(x = year, ymin = S_lwr, ymax = S_upr, fill = HCR),alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(spwn_df, year >=d_start, year <= d_end), aes(x = year, y = mid),
            lwd = 1.2) +
  geom_ribbon(data = filter(spwn_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  geom_hline(yintercept = benchmarks[1,1]) +
  geom_line(aes(x = year, y = S, color = HCR), lwd = 1.2) +
  geom_line(data = ind.sims, aes(x = year, y = spawners, lty = sim, color = HCR)) +
  annotate("text", x = 2020, y = benchmarks[1,1]+.4,
           label = expression(italic(paste("80%",S)[MSY]))) +
  geom_hline(yintercept = benchmarks[2,1]) +
  annotate("text", x = 2020, y = 2.2,
           label = "italic(S[gen])", parse = TRUE) +
  scale_x_continuous(breaks= pretty_breaks(),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "Return year", y = "Escapement") +
  scale_fill_viridis_d(name = "HCR") +
  scale_color_viridis_d(name = "HCR") +
 # scale_linetype_manual(values=c(1,1,1,1)) + #hack to get lines to stay the same since group arg is broken
  theme(legend.position = "bottom") +
  guides(lty = "none")

#catch ---
p2 <- ggplot(data = filter(fwd.sim, scenario == "baseline")) +
  #draw the fwd.sim
  geom_ribbon(aes(x = year, ymin = C_lwr, ymax = C_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(C_df, year >=d_start, year <= d_end), aes(x = year, y = mid),
            lwd = 1.2) +
  geom_ribbon(data = filter(C_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  geom_line(aes(x = year, y = C, color = HCR), lwd = 1.2) +
  geom_line(data = ind.sims, aes(x = year, y = catch, lty = sim, color = HCR)) +
  geom_hline(yintercept = rel.catch.index, lty = 2) +
  annotate("text", x = d_start + 4, y = rel.catch.index+.2,
           label = "catch index") +
  scale_x_continuous(breaks= pretty_breaks(),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "Return year", y = "Catch (M)") +
  scale_fill_viridis_d(name = "HCR") +
  scale_color_viridis_d(name = "HCR") +
  theme(legend.position = "bottom") +
  guides(lty = "none")

p <- plot_grid(p1, p2) #+ draw_grob(legend) #fix to make a single legend?
p
my.ggsave(here("figure/fwd-SC.png"))

# APPENDIX FIGS? -------------------------------------------------------------------------
if(FALSE){
#and fwd sim harvest rate
ggplot(data = fwd.sim) +
  #draw the fwd.sim
  geom_line(aes(x = year, y = U, color = HCR), lwd = 1, lty = 2) +
  geom_ribbon(aes(x = year, ymin = U_mid_lwr, ymax = U_mid_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(er_df, year >=d_start, year <= d_end), aes(x = year, y = mid)) +
  geom_ribbon(data = filter(er_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(data = filter(er_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  #draw the benchmarks
  geom_hline(yintercept = benchmarks$median[3]) +
  annotate("text", x = d_start + 0.5, y = 0.55,
           label = "italic(U[MSY])", parse = TRUE) +
  scale_x_continuous(breaks= pretty_breaks()) +
  #coord_cartesian(expand = FALSE) +
  labs(x = "Return year", y = "Harvest rate") +
  scale_fill_viridis_d(name = "HCR",
                       breaks = c("current", "alt", "WSP"),
                       labels = c("current", "alternate", "WSP")) +
  scale_color_viridis_d(name = "HCR",
                        breaks = c("current", "alt", "WSP"),
                        labels = c("current", "alternate", "WSP")) +
  facet_grid(prod~.) +
  theme(legend.position = "bottom")

my.ggsave(here("figure/fwd-U.png"))

#and recruits
ggplot(data = fwd.sim) +
  #draw the fwd.sim
  geom_line(aes(x = year, y = R, color = HCR), lwd = 1, lty = 2) +
  #geom_ribbon(aes(x = year, ymin = R_lwr, ymax = R_upr,  fill = HCR), alpha = 0.5) +
  geom_ribbon(aes(x = year, ymin = R_mid_lwr, ymax = R_mid_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(R_df, year >=d_start, year <= d_end), aes(x = year, y = mid)) +
  geom_ribbon(data = filter(R_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(data = filter(R_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  #draw the benchmarks
  #geom_rect(aes(xmin = d_start, xmax = max(fwd.sim$year),
  #              ymin = benchmarks$`lower 95% CI`[2], ymax = benchmarks$`upper 95% CI`[2]),
  #          fill = "darkred", alpha = 0.02) +
  geom_hline(yintercept = benchmarks$median[2], color = "darkred") +
  annotate("text", x = d_start + 2, y = 4,
           label = "italic(S[gen])", parse = TRUE, color = "darkred") +
  scale_x_continuous(breaks= pretty_breaks()) +
  #coord_cartesian(expand = FALSE) +
  labs(x = "return year", y = "recruits") +
  scale_fill_viridis_d(name = "HCR",
                       breaks = c("current", "alt", "WSP"),
                       labels = c("current", "alternate", "WSP")) +
  scale_color_viridis_d(name = "HCR",
                        breaks = c("current", "alt", "WSP"),
                        labels = c("current", "alternate", "WSP")) +
  facet_grid(prod~.) +
  theme(legend.position = "bottom")

my.ggsave(here("figure/fwd-R.png"))
}

#make objs to read into text
model.summary <- as.data.frame(rstan::summary(stan.fit)$summary)
model.summary.93 <- as.data.frame(rstan::summary(stan.fit.93)$summary)

rm(p, p1, p2, C_df, C.quant, R_df, R.quant, resid.quant, resids)
