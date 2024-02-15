library(tidyverse)
library(here)
library(gsl)
library(patchwork)
library(scales) #for pretty_breaks() on fig axes
source(here("analysis/R/functions.R"))
source(here("analysis/R/fwd-sim.R"))

# read in data and fit -------------------------------------------------------------------
current.hcr <- read.csv(here("analysis/data/raw/current_hcr.csv"))

avg_mass <- read.csv(here("analysis/data/raw/bio/HistoricalPinkWeightDownload.csv"))
colnames(avg_mass) <- c("year", "avg.weight", "reference")

compitetors <- read.csv(here("analysis/data/raw/bio/competitor_density_long.csv")) |>
  filter(species == "pink", area == "np")

benchmarks <- as.data.frame(benchmarks) |>
  round(2)

#WRANGLING -------------------------------------------------------------------------------
#latent states of spawners and recruits---
spwn.quant <- apply(model.pars$S, 2, quantile, probs=c(0.025,0.5,0.975))[,1:32]
rec.quant <- apply(model.pars$R, 2, quantile, probs=c(0.025,0.5,0.975))[,2:33]

brood_t <- as.data.frame(cbind(data$year[1:32],t(spwn.quant), t(rec.quant))) |>
  round(2)
colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")

#SR relationship---
spw <- seq(0,max(brood_t$S_upr),length.out=100)
max_samples <- dim(model.pars$lnalpha)
SR_pred <- matrix(NA,100,max_samples)

for(i in 1:max_samples){
  a <- model.pars$lnalpha[i]
  b <- model.pars$beta[i]
  SR_pred[,i] <- (exp(a)*spw*exp(-b*spw))
}

SR_pred <- as.data.frame(cbind(spw,t(apply(SR_pred,c(1),quantile,probs=c(0.025,0.5,0.975),na.rm=T))))|>
  round(2)
colnames(SR_pred) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr")

#create spawner & esc df --- #should map() these...
spwn.quant <- apply(model.pars$S, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
spwn_df <- as.data.frame(cbind(data$year, t(spwn.quant)))
colnames(spwn_df) <- c("year","lwr","mid_lwr","mid","mid_upr","upr")

er.quant <- apply(model.pars$U, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
er_df <- as.data.frame(cbind(data$year, t(er.quant)))
colnames(er_df) <- c("year","lwr","mid_lwr","mid","mid_upr","upr")

R.quant <- apply(model.pars$R, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
R_df <- as.data.frame(cbind(data$year, t(R.quant)))
colnames(R_df) <- c("year","lwr","mid_lwr","mid","mid_upr","upr")

C.quant <- apply(model.pars$C, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))
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
  mutate(U_Umsy = round(U/benchmarks[3,1] ,2), ## index these as names!
         U_Umsy_LCI = round(U_LCI/benchmarks[3,1] ,2),
         U_Umsy_UCI = round(U_UCI/benchmarks[3,1] ,2),
         S_Smsy = round(S/(benchmarks[1,1]*1.25), 2), #correct from 80%Smsy to Smsy
         S_Smsy_LCI = round(S_LCI/(benchmarks[1,1]*1.25), 2),
         S_Smsy_UCI = round(S_UCI/(benchmarks[1,1]*1.25), 2))

rm(er.quant, rec.quant, spwn.quant, max_samples, percentiles, a,b, spw)

# INFERENCE ------------------------------------------------------------------------------
#get DATA (not latent state) of avg run sizes for context/background section
rec_prop <- data |>
  mutate(rec = harvest + spawn,
         ER = (harvest/rec)*100)

# FIGS & TABLES --------------------------------------------------------------------------
theme_set(theme_bw())

#wrapper function to set dims before saving, so things scale proper each time
  #ASSUMED this is the correct scaling based on index.Rmd
my.ggsave <- function(filename = default_name(plot), width= 9, height = 5.562, dpi= 180){
  ggsave(filename=filename, width=width, height=height, dpi=dpi)
}

# plot catch & escapement ---
catch_esc <- data |>
  select(year, harvest, spawn) |>
  pivot_longer(!year, names_to = "type", values_to = "n")

ggplot(catch_esc, aes(x = year, y = n, fill = type)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("darkgrey", "black"), name = "Return type") +
  theme(legend.position = "bottom") +
  labs(x = "Return year",
       y = "Total return (millions of fish)")

my.ggsave(here("figure/catch-esc.png"))

#plot avg mass through time
ggplot(avg_mass, aes(year, scale(avg.weight))) +
  geom_line() +
  geom_point() +
  geom_point(data = compitetors, aes(x = Year, y = scale(value)), color = "red") +
  stat_smooth(data = compitetors, aes(x = Year, y = scale(value)), color = "red",
              se = FALSE) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "N. Pacific Pink abundance (scaled)")) +
  labs(x = "Year", y = "Average body mass (scaled)")

my.ggsave(here("figure/avg-mass.png"))

# plot current HCR -----------------------------------------------------------------------
  # does it make sense to slap ref pts on here? AMH say no - they're OCPs
p1 <- ggplot(current.hcr, aes(x=run_size, y=esc)) +
  geom_line(size=1.1) +
  geom_vline(xintercept = 7.059, lty = 2)+
  geom_vline(xintercept = 20, lty = 2)+
  labs(x = "Run size (millions)",
       y = "Target escapement")

p2 <- ggplot(current.hcr, aes(x=run_size, y=er)) +
  geom_line(size=1.1) +
  geom_vline(xintercept = 7.059, lty = 2)+
  geom_vline(xintercept = 20, lty = 2)+
  ylim(c(0,1)) +
  labs(x = "Run size (millions)",
       y = "Target exploitation rate")

p1/p2

my.ggsave(here("figure/current-HCR.png"))

# plot SR relationship -------------------------------------------------------------------
ggplot() +
  geom_ribbon(data = SR_pred, aes(x = Spawn, ymin = Rec_lwr, ymax = Rec_upr),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_line(data = SR_pred, aes(x = Spawn, y = Rec_med), color="black", size = 1) +
  geom_errorbar(data = brood_t, aes(x= S_med, y = R_med, ymin = R_lwr, ymax = R_upr),
                colour="grey", width=0, size=0.3) +
  geom_errorbarh(data = brood_t, aes(y = R_med, xmin = S_lwr, xmax = S_upr),
                 height=0, colour = "grey", height = 0, size = 0.3) +
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

# then residuals
resid.quant <- apply(model.pars$lnresid, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))[,1:32]

resids <- as.data.frame(cbind(data$year[1:32], t(resid.quant)))
colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr")

ggplot(resids, aes(x=year, y = mid)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
  geom_line(lwd = 1.1) +
  coord_cartesian(ylim=c(-2,2)) +
  labs(x = "Brood year",
       y = "Recruitment residuals") +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)

my.ggsave(here("figure/rec-resid.png"))

#plot fwd sim of escapement --------------------------------------------------------------
#how much of the old data do you want to show?
d_start <- 2017
d_end <- 2023

#spawners
ggplot(data = fwd.sim) +
  #draw the fwd.sim
  geom_line(aes(x = year, y = S, color = HCR), lwd = 1, lty = 2) +
  #geom_ribbon(aes(x = year, ymin = S_lwr, ymax = S_upr,  fill = HCR), alpha = 0.5) +
  geom_ribbon(aes(x = year, ymin = S_mid_lwr, ymax = S_mid_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(spwn_df, year >=d_start, year <= d_end), aes(x = year, y = mid)) +
  geom_ribbon(data = filter(spwn_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(data = filter(spwn_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  #draw the benchmarks
 # geom_rect(aes(xmin = d_start, xmax = max(fwd.sim$year),
 #               ymin = benchmarks$`lower 95% CI`[1], ymax = benchmarks$`upper 95% CI`[1]),
 #           fill = "forestgreen", alpha = 0.02) +
  geom_hline(yintercept = benchmarks$median[1],
             color = "forestgreen") +
  annotate("text", x = d_start + 2, y = 5,
           label = "italic(S[MSY])", parse = TRUE, color = "forestgreen") +
#  geom_rect(aes(xmin = d_start, xmax = max(fwd.sim$year),
#                ymin = benchmarks$`lower 95% CI`[2], ymax = benchmarks$`upper 95% CI`[2]),
#            fill = "darkred", alpha = 0.02) +
  geom_hline(yintercept = benchmarks$median[2], color = "darkred") +
  annotate("text", x = d_start + 2, y = 2,
           label = "italic(S[gen])", parse = TRUE, color = "darkred") +
  scale_x_continuous(breaks= pretty_breaks()) +
  #coord_cartesian(expand = FALSE) +
  labs(x = "return year", y = "escapement") +
  scale_fill_viridis_d(name = "HCR",
                       breaks = c("current", "alt", "WSP"),
                       labels = c("current", "alternate", "WSP")) +
  scale_color_viridis_d(name = "HCR",
                        breaks = c("current", "alt", "WSP"),
                        labels = c("current", "alternate", "WSP")) +
  facet_grid(prod~.) +
  theme(legend.position = "bottom")

my.ggsave(here("figure/fwd-S.png"))

#and fwd sim harvest rate
ggplot(data = fwd.sim) +
  #draw the fwd.sim
  geom_line(aes(x = year, y = U, color = HCR), lwd = 1, lty = 2) +
  #geom_ribbon(aes(x = year, ymin = U_lwr, ymax = U_upr,  fill = HCR), alpha = 0.5) +
  geom_ribbon(aes(x = year, ymin = U_mid_lwr, ymax = U_mid_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(er_df, year >=d_start, year <= d_end), aes(x = year, y = mid)) +
  geom_ribbon(data = filter(er_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(data = filter(er_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  #draw the benchmarks
#  geom_rect(aes(xmin = d_start, xmax = max(fwd.sim$year),
#                ymin = benchmarks$`lower 95% CI`[3], ymax = benchmarks$`upper 95% CI`[3]),
#            fill = "forestgreen", alpha = 0.02) +
  geom_hline(yintercept = benchmarks$median[3],
             color = "forestgreen") +
  annotate("text", x = d_start + 0.5, y = 0.55,
           label = "italic(U[MSY])", parse = TRUE, color = "forestgreen") +
  scale_x_continuous(breaks= pretty_breaks()) +
  #coord_cartesian(expand = FALSE) +
  labs(x = "return year", y = "harvest rate") +
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

#and fwd sim catch
ggplot(data = fwd.sim) +
  #draw the fwd.sim
  geom_line(aes(x = year, y = C, color = HCR), lwd = 1, lty = 2) +
  #geom_ribbon(aes(x = year, ymin = U_lwr, ymax = U_upr,  fill = HCR), alpha = 0.5) +
  geom_ribbon(aes(x = year, ymin = C_mid_lwr, ymax = C_mid_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(C_df, year >=d_start, year <= d_end), aes(x = year, y = mid)) +
  geom_ribbon(data = filter(C_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(data = filter(C_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  scale_x_continuous(breaks= pretty_breaks()) +
  #coord_cartesian(expand = FALSE) +
  labs(x = "return year", y = "catch (M)") +
  scale_fill_viridis_d(name = "HCR",
                       breaks = c("current", "alt", "WSP"),
                       labels = c("current", "alternate", "WSP")) +
  scale_color_viridis_d(name = "HCR",
                        breaks = c("current", "alt", "WSP"),
                        labels = c("current", "alternate", "WSP")) +
  facet_grid(prod~.) +
  theme(legend.position = "bottom")

my.ggsave(here("figure/fwd-C.png"))

# performance metrics---------------------------------------------------------------------

perf.metrics <- perf.metrics |>
  pivot_longer(-c(yrs, metric, prod), names_to = "HCR") |>
  mutate(bin = case_when(value >= 0 & value <=.11 ~ "low",
                         value > .11 & value <=.34 ~ "lowmed",
                         value > .34 & value <= .65 ~ "med",
                         value > .65 ~ "high"),
         yrs = factor(paste(yrs, "generations"), levels = c("3 generations",
                                                            "10 generations")))

ggplot(data = filter(perf.metrics, metric == "above.Smsy.8"), aes(x=prod, y = HCR)) +
  geom_raster(aes(fill= bin)) +
  geom_text(aes(label = round(value,2))) +
  scale_fill_manual(values = c("forestgreen", "white", "pink", "hotpink")) +
  facet_grid(.~yrs) +
  guides(fill = "none") +
  labs(x = "productivity", title = bquote("Percent of years spawners is above S"[MSY]))

my.ggsave(here("figure/perf-Smsy.png"))

ggplot(data = filter(perf.metrics, metric == "below.Sgen"), aes(x=prod, y = HCR)) +
  geom_raster(aes(fill= bin)) +
  geom_text(aes(label = round(value,2))) +
  scale_fill_manual(values = c("forestgreen", "white", "pink", "hotpink")) +
  facet_grid(.~yrs) +
  guides(fill = "none") +
  labs(x = "productivity", title = bquote("Percent of years spawners is below S"[gen]))

my.ggsave(here("figure/perf-Sgen.png"))

ggplot(data = filter(perf.metrics, metric == "catch"), aes(x=prod, y = HCR)) +
  geom_raster(aes(fill= value)) +
  geom_text(aes(label = round(value,2))) +
  scale_fill_viridis_c() +
  facet_grid(.~yrs) +
  labs(x = "productivity", title = "Catch (millions of fish) *bad scales but can fix*")+
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(title="catch (M)"))

my.ggsave(here("figure/perf-catch.png"))

ggplot(data = filter(perf.metrics, metric == "catch.stability"), aes(x=prod, y = HCR)) +
  geom_raster(aes(fill= value)) +
  geom_text(aes(label = round(value,2))) +
  scale_fill_viridis_c() +
  facet_grid(.~yrs) +
  labs(x = "productivity", title = "Catch stability (1/CV)") +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(title="stability"))

my.ggsave(here("figure/perf-catch-stab.png"))
