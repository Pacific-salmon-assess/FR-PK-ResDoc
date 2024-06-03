source(here("analysis/R/inference-figs.R"))
#might need to source fwd sim to get bigger sim object to pull draws from

#make the 3
# initial ---
ggplot(data = fwd.sim) +
  #draw the fwd.sim
  geom_line(aes(x = year, y = S, color = prod), lwd = 1, lty = 2) +
  geom_ribbon(aes(x = year, ymin = S_mid_lwr, ymax = S_mid_upr, fill = prod), alpha=0.2) +
  #draw the exsiting data
  geom_line(data = filter(spwn_df, year >=d_start, year <= d_end), aes(x = year, y = mid)) +
  geom_ribbon(data = filter(spwn_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  geom_hline(yintercept = benchmarks$median[1]) +
  annotate("text", x = d_start + 2, y = 5,
           label = "italic(S[MSY])", parse = TRUE) +
  geom_hline(yintercept = benchmarks$median[2]) +
  annotate("text", x = d_start + 2, y = 2,
           label = "italic(S[gen])", parse = TRUE) +
  scale_x_continuous(breaks= pretty_breaks(),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "return year", y = "escapement") +
  scale_fill_viridis_d(name = "productivity") +
  scale_color_viridis_d(name = "productivity") +
  facet_grid(HCR~.) +
  theme(legend.position = "bottom")

#1 pane show HCRs ---
ggplot(data = filter(fwd.sim, prod == "baseline")) +
  #draw the fwd.sim
  geom_line(aes(x = year, y = S, color = HCR), lwd = 1, lty = 2) +
  geom_ribbon(aes(x = year, ymin = S_mid_lwr, ymax = S_mid_upr, fill = HCR, color = HCR), alpha=0.2, lty = 2) +
  #draw the exsiting data
  geom_line(data = filter(spwn_df, year >=d_start, year <= d_end), aes(x = year, y = mid)) +
  geom_ribbon(data = filter(spwn_df, year >=d_start, year <= d_end),
              aes(x = year, ymin = mid_lwr, ymax = mid_upr), fill = "black", alpha=0.2) +
  geom_hline(yintercept = benchmarks$median[1]) +
  annotate("text", x = d_start + 2, y = 5,
           label = "italic(S[MSY])", parse = TRUE) + #hard time putting "80%" here
  geom_hline(yintercept = benchmarks$median[2]) +
  annotate("text", x = d_start + 2, y = 2,
           label = "italic(S[gen])", parse = TRUE) +
  scale_x_continuous(breaks= pretty_breaks(),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = "return year", y = "escapement") +
  scale_fill_viridis_d(name = "HCR") +
  scale_color_viridis_d(name = "HCR") +
  theme(legend.position = "bottom")

#3 pane show draws ---