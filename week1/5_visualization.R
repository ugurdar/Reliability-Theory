# Week 1: Visualization
source("week1/4_solutions.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

# Time grids
t_max <- 30
t_seq <- seq(0, t_max, length.out = 1000)
t_pos <- seq(0.001, t_max, length.out = 1000)

# Recompute needed scalars (in case sourced env changes)
median_life <- find_median()
mttf <- calculate_mttf()

# Data
df_h <- data.frame(t = t_pos, h = hazard_rate(t_pos))
df_R <- data.frame(t = t_seq, R = reliability(t_seq))
df_f <- data.frame(t = t_pos, f = pdf_func(t_pos))
df_F <- data.frame(t = t_seq, F = cdf_func(t_seq))
df_m <- data.frame(t = t_seq, m = mean_residual_life(t_seq))

# Plots
p1 <- ggplot(df_h, aes(t, h)) +
  geom_line(color = "#E31A1C", linewidth = 1.2) +
  labs(title = "Hazard (z(t))", x = "t", y = "z(t)") +
  theme_minimal()

p2 <- ggplot(df_R, aes(t, R)) +
  geom_line(color = "#1F78B4", linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = median_life, linetype = "dashed", color = "grey50") +
  labs(title = "Reliability R(t)", subtitle = paste("Median ≈", round(median_life, 3)),
       x = "t", y = "R(t)") +
  theme_minimal()

p3 <- ggplot(df_f, aes(t, f)) +
  geom_line(color = "#33A02C", linewidth = 1.2) +
  geom_vline(xintercept = 0, color = "red") +
  labs(title = "PDF f(t)", subtitle = "Monotonically decreasing; mode at t = 0",
       x = "t", y = "f(t)") +
  theme_minimal()

p4 <- ggplot(df_F, aes(t, F)) +
  geom_line(color = "#FF7F00", linewidth = 1.2) +
  labs(title = "CDF F(t)", x = "t", y = "F(t)") +
  theme_minimal()

p5 <- ggplot(df_m, aes(t, m)) +
  geom_line(color = "#6A3D9A", linewidth = 1.2) +
  labs(title = "Mean Residual Life m(t) = t + 5", x = "t", y = "m(t)") +
  theme_minimal()

grid <- grid.arrange(p1, p2, p3, p4, p5, ncol = 2)
ggsave("week1/all.png", grid, width = 12, height = 12, dpi = 300)

ggsave("week1/plot_hazard.png", p1, width = 8, height = 6, dpi = 300)
ggsave("week1/plot_reliability.png", p2, width = 8, height = 6, dpi = 300)
ggsave("week1/plot_pdf.png", p3, width = 8, height = 6, dpi = 300)
ggsave("week1/plot_cdf.png", p4, width = 8, height = 6, dpi = 300)
ggsave("week1/plot_mrlt.png", p5, width = 8, height = 6, dpi = 300)

cat("Saved plots to week1/*.png\n")