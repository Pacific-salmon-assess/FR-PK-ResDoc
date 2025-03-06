library(ggsidekick)

#wrapper function to set dims before saving, so things scale proper each time
#ASSUMED this is the correct scaling based on index.Rmd
my.ggsave <- function(filename = default_name(plot), width= 6, height = 3.5, dpi= 180){
  ggsave(filename=filename, width=width, height=height, dpi=dpi)
}

my.ggsave.tall <- function(filename = default_name(plot), width= 5.5, height = 4.5, dpi= 180){
  ggsave(filename=filename, width=width, height=height, dpi=dpi)
}

my.ggsave.long <- function(filename = default_name(plot), width= 11, height = 5.562, dpi= 180){
  ggsave(filename=filename, width=width, height=height, dpi=dpi)
}

# plot catch & escapement ---
catch_esc <- data |>
  select(year, harvest, spawn) |>
  pivot_longer(!year, names_to = "type", values_to = "n")

ggplot(catch_esc, aes(x = year, y = n)) +
  geom_bar(position="stack", stat="identity", aes(fill = type)) +
  scale_fill_manual(values = c("darkgrey", "black"), name = "Return type") +
  scale_x_continuous(breaks = c(1961, 1981, 2001, 2021)) +
  theme(legend.position = "bottom") +
  theme_sleek() +
  geom_vline(xintercept = c(1990, 2000, 2006), lty=2) +
  labs(x = "Return year",
       y = "Total return (millions of fish)")

my.ggsave(here("figure/AFS/catch-esc.png"))


# average body mass and competitors ---
ggplot(avg_mass, aes(year, avg.weight)) +
  geom_col(data = compitetors, aes(x = return_yr, y = n_pink/200), fill = "pink", color = "pink") +
  geom_line() +
  geom_point() +
  ylim(1,3) +
  scale_y_continuous(sec.axis = sec_axis(~.*200, name = "N. Pacific Pink abundance (millions)")) +
  theme_sleek() +
  labs(x = "Year", y = "Average body mass (kg)")

my.ggsave(here("figure/AFS/avg-mass.png"))

# plot HCRs ---

curr.HCR <- HCRs |>
  filter(HCR == "current")

p1 <- ggplot(curr.HCR, aes(x=run_size, y=ER)) +
  geom_line(linewidth=1.1, alpha = 0.7, col = "#E69F00") +
  ylim(c(0,1)) +
  labs(x = NULL,
       y = "Target exp. rate") +
  theme_sleek() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
  
  

p2 <- ggplot(curr.HCR, aes(x=run_size, y=esc_goal)) +
  geom_line(linewidth=1.1, alpha = 0.7, col = "#E69F00") +
  labs(x = "Run size (millions)",
       y = "Target spawners") +
  ylim(0,10) +
  theme_sleek() 

ggarrange(p1, p2, nrow = 2,
          align = "hv", common.legend = TRUE, legend = "bottom")
my.ggsave.tall(here("figure/AFS/curr_HCR.png"))

p1 <- ggplot(HCRs, aes(x=run_size, y=ER, color = HCR)) +
  geom_line(linewidth=1.1, alpha = 0.7) +
  geom_vline(xintercept = R.Smsy.8, lty = 2) +
  annotate("text", x = R.Smsy.8+2, y = .7,
           label = expression(italic(R[paste("80%",S)[MSY]]))) +
  geom_vline(xintercept = Sgen, lty =2) +
  annotate("text", x = Sgen+1, y = .7,
           label = "italic(LRP)", parse = TRUE) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  ylim(c(0,1)) +
  labs(x = NULL,
       y = "Target exp. rate") +
  theme_sleek() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none") 

p2 <- ggplot(HCRs, aes(x=run_size, y=esc_goal, color = HCR)) +
  geom_line(linewidth=1.1, alpha = 0.7) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  geom_vline(xintercept = R.Smsy.8, lty = 2) +
  geom_vline(xintercept = Sgen, lty = 2) +
  labs(x = "Run size (millions)",
       y = "Target spawners") +
  theme_sleek() +
  theme(legend.position="none")

ggarrange(p1, p2, nrow = 2,
          align = "hv")
my.ggsave.tall(here("figure/AFS/HCRs.png"))

# time varying alpha ---
ggplot(as.data.frame(a_yrs)) +
  geom_ribbon(aes(x = brood_year, ymin = lwr, ymax = upr), fill = "darkgrey", alpha = 0.5) +
  geom_line(aes(x = brood_year, y = mid), lwd = 2,  color = "black") +
  scale_x_continuous(breaks = c(1961, 1981, 2001, 2021)) +
  labs(y = "Productivity (max. R/S)", x = "Brood year") +
  theme_sleek() 

my.ggsave(here("figure/AFS/tv-alpha.png"))

# plot the last distribution of spawners and the posteriors of Sgen and Smsy
#just use full posteriors from model pars
ggplot() +
  geom_density(data = data.frame(AR1.model.pars$S[,33]), aes(x=AR1.model.pars.S...33.), fill = "grey", alpha = 0.2, bw=0.5) +
  geom_density(data = data.frame(Smsy.8.post), aes(x=Smsy.8.post), fill = "forestgreen",
               alpha = 0.2, color = "forestgreen", bw=0.5) +
  geom_density(data = data.frame(Sgen.post), aes(x=Sgen.post), fill = "darkred",
               alpha = 0.2, color = "darkred", bw=0.5)  +
  annotate("text", x = 5, y = 0.1,
           label = expression(italic(paste("80%",S)[MSY])), size = 4, color = "forestgreen") +
  annotate("text", x = 2, y = 0.1,
           label = "italic(LRP)", parse = TRUE, size = 4, color = "darkred") +
  annotate("text", x = 9.5, y = 0.1,
           label = "italic(S[23])", parse = TRUE, size = 4) +
  #annotate("text", x = 9.5, y = 0.1,
  #         label = "italic(S[19-23])", parse = TRUE, size = 5) +
  coord_cartesian(xlim=c(0,15)) +
  labs(y = "Posterior density", x = "Spawners (millions)")+
  theme_sleek() 

my.ggsave(here("figure/AFS/recent-status.png"))

#plot fwd sims of spawners & catch ---
#how much of the old data do you want to show?
d_start <- 2015
d_end <- 2023

#spawners---
p1 <- ggplot(data = filter(fwd.sim, scenario == "base", HCR != "no_fishing" )) +
  #draw the fwd.sim
  geom_ribbon(aes(x = year, ymin = S_lwr, ymax = S_upr, fill = HCR),alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(spwn_df, year >=d_start, year <= d_end), aes(x = year, y = mid),
            lwd = 1.2) +
  geom_ribbon(data = filter(spwn_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr), fill = "black", alpha=0.2) +
  geom_hline(yintercept = benchmarks[1,1], lty =2) +
  geom_line(aes(x = year, y = S, color = HCR), lwd = 1.2) +
  annotate("text", x = 2020, y = benchmarks[2,1]+0.6,
           label = expression(italic(paste("80%",S)[MSY])), size = 3) +
  geom_hline(yintercept = benchmarks[2,1], lty = 2) +
  annotate("text", x = 2020, y = 2.3,
           label = "italic(LRP)", parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(2015, 2019, 2023, 2027, 2031),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "", y = "Spawners") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  scale_linetype_manual(values=c(1,1))  +
  theme_sleek() +
  theme(legend.position = "none") +
  guides(lty = "none")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#catch ---
p2 <- ggplot(data = filter(fwd.sim, scenario == "base", HCR != "no_fishing" )) +
  #draw the fwd.sim
  geom_ribbon(aes(x = year, ymin = C_lwr, ymax = C_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(C_df, year >=d_start, year <= d_end), aes(x = year, y = mid),
            lwd = 1.2) +
  geom_ribbon(data = filter(C_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr), fill = "black", alpha=0.2) +
  geom_line(aes(x = year, y = C, color = HCR), lwd = 1.2) +
  scale_x_continuous(breaks = c(2015, 2019, 2023, 2027, 2031),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "Return year", y = "Catch") +
  scale_color_viridis_d(labels = c("current",  "PA alternate")) +
  scale_fill_viridis_d(labels = c("current", "PA alternate")) +
  scale_linetype_manual(values=c(1,1))  +
  theme_sleek() +
  theme(legend.position = "bottom") +
  guides(lty = "none")

p <- plot_grid(p1, p2, nrow = 2) #+ draw_grob(legend) #fix to make a single legend?

p
ggsave(here("figure/AFS/fwd-SC.png"), width= 5, height = 5, dpi= 180)


#spawners---
p1 <- ggplot(data = filter(fwd.sim, scenario == "low_prod", HCR != "no_fishing" )) +
  #draw the fwd.sim
  geom_ribbon(aes(x = year, ymin = S_lwr, ymax = S_upr, fill = HCR),alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(spwn_df, year >=d_start, year <= d_end), aes(x = year, y = mid),
            lwd = 1.2) +
  geom_ribbon(data = filter(spwn_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr), fill = "black", alpha=0.2) +
  geom_hline(yintercept = benchmarks[1,1], lty =2) +
  geom_line(aes(x = year, y = S, color = HCR), lwd = 1.2) +
  annotate("text", x = 2020, y = benchmarks[2,1]+0.6,
           label = expression(italic(paste("80%",S)[MSY])), size = 3) +
  geom_hline(yintercept = benchmarks[2,1], lty = 2) +
  annotate("text", x = 2020, y = 2.3,
           label = "italic(LRP)", parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(2015, 2019, 2023, 2027, 2031),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "", y = "Spawners") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  scale_linetype_manual(values=c(1,1))  +
  theme_sleek() +
  theme(legend.position = "none") +
  guides(lty = "none")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#catch ---
p2 <- ggplot(data = filter(fwd.sim, scenario == "low_prod", HCR != "no_fishing" )) +
  #draw the fwd.sim
  geom_ribbon(aes(x = year, ymin = C_lwr, ymax = C_upr, fill = HCR), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(C_df, year >=d_start, year <= d_end), aes(x = year, y = mid),
            lwd = 1.2) +
  geom_ribbon(data = filter(C_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = lwr, ymax = upr), fill = "black", alpha=0.2) +
  geom_line(aes(x = year, y = C, color = HCR), lwd = 1.2) +
  scale_x_continuous(breaks = c(2015, 2019, 2023, 2027, 2031),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "Return year", y = "Catch") +
  scale_color_viridis_d(labels = c("current",  "PA alternate")) +
  scale_fill_viridis_d(labels = c("current", "PA alternate")) +
  scale_linetype_manual(values=c(1,1))  +
  theme_sleek() +
  theme(legend.position = "bottom") +
  guides(lty = "none")

p <- plot_grid(p1, p2, nrow = 2) #+ draw_grob(legend) #fix to make a single legend?

p
ggsave(here("figure/AFS/fwd-SC_low.png"), width= 5, height = 5, dpi= 180)