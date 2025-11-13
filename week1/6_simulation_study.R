# Week 1: Simulation Study - 50 Unit Sample
source("week1/4_solutions.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

set.seed(2025)
cat("Simulation Study: 50-Unit Sample Analysis\n")
cat(paste(rep("=", 45), collapse = ""), "\n\n")

# Inverse transform sampler for F(t) = 1 - 1/(0.2 t + 1)^2
# Derivation:
#   U = F(t) = 1 - 1/(0.2 t + 1)^2
#   1 - U = 1/(0.2 t + 1)^2
#   (0.2 t + 1)^2 = 1 / (1 - U)
#   0.2 t + 1 = 1 / sqrt(1 - U)
#   t = (1 / sqrt(1 - U) - 1)/0.2 = 5 * (1 / sqrt(1 - U) - 1)
inverse_transform <- function(u) {
  u <- pmin(u, 0.999999)
  5 * (1 / sqrt(1 - u) - 1)
}

# Sample
n_sample <- 100
u_vals <- runif(n_sample)
sample_data <- inverse_transform(u_vals)
sample_data <- sample_data[is.finite(sample_data) & sample_data >= 0]
n_actual <- length(sample_data)

cat("Generated:", n_actual, "failure times\n")
cat("Range: [", round(min(sample_data), 3), ",", round(max(sample_data), 3), "]\n\n")

# Empirical stats
s_mean <- mean(sample_data)
s_var  <- var(sample_data)
s_sd   <- sd(sample_data)
s_med  <- median(sample_data)
s_min <- min(sample_data)
# Theoretical (finite ones only)
theoretical_mttf <- calculate_mttf()
theoretical_med  <- find_median()
theoretical_var  <- calculate_variance()   # Inf
theoretical_sd   <- calculate_std_dev()    # Inf

# Comparison table via base R
cmp <- data.frame(
  Metric = c("MTTF", "Median", "Variance", "Std_Dev","min"),
  Theoretical = c(theoretical_mttf, theoretical_med, theoretical_var, theoretical_sd,0),
  Sample = c(s_mean, s_med, s_var, s_sd,s_min),
  stringsAsFactors = FALSE
)
print(cmp)
cmp$Absolute_Error <- ifelse(is.finite(cmp$Theoretical), abs(cmp$Sample - cmp$Theoretical), NA_real_)
cmp$Relative_Error_Pct <- ifelse(is.finite(cmp$Theoretical) & cmp$Theoretical != 0,
                                 100 * cmp$Absolute_Error / abs(cmp$Theoretical), NA_real_)

write.csv(cmp, "week1/6_simulation_zcomparison_table.csv", row.names = FALSE)

# KS test
theoretical_cdf <- function(t) cdf_func(t)
ks_res <- ks.test(sample_data, theoretical_cdf)

# Plots (use data.frames to avoid aes vector mapping issues)
t_plot <- seq(0, max(sample_data) * 1.2, length.out = 600)

p_hist <- ggplot(data.frame(x = sample_data)) +
  geom_histogram(aes(x = x, y = ..density..), bins = 15,
                 fill = "steelblue", color = "black", alpha = 0.7) +
  geom_line(data = data.frame(t = t_plot, pdf = pdf_func(t_plot)),
            aes(x = t, y = pdf), color = "red", linewidth = 1.2) +
  geom_vline(xintercept = theoretical_mttf, linetype = "dashed", color = "red") +
  geom_vline(xintercept = s_mean, linetype = "dashed", color = "blue") +
  labs(title = "Histogram + Theoretical PDF", x = "t", y = "density") +
  theme_minimal()

emp_ecdf <- ecdf(sample_data)
t_cdf <- seq(0, max(sample_data), length.out = 600)
p_cdf <- ggplot() +
  geom_step(data = data.frame(t = t_cdf, ecdf = emp_ecdf(t_cdf)),
            aes(x = t, y = ecdf), color = "blue", linewidth = 1.1) +
  geom_line(data = data.frame(t = t_cdf, tcdf = theoretical_cdf(t_cdf)),
            aes(x = t, y = tcdf), color = "red", linewidth = 1.1) +
  labs(title = "Empirical vs Theoretical CDF", x = "t", y = "CDF") +
  theme_minimal()

# QQ plot with distribution quantiles
prob_seq <- (1:n_actual - 0.5) / n_actual
q_theo <- 5 * (1 / sqrt(1 - prob_seq) - 1)
q_samp <- sort(sample_data)
p_qq <- ggplot(data.frame(theoretical = q_theo, sample = q_samp),
               aes(theoretical, sample)) +
  geom_point(color = "darkblue", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()



dash <- grid.arrange(p_hist, p_cdf, p_qq, ncol = 2)
ggsave("week1/6_simulation_zall.png", dash, width = 12, height = 10, dpi = 300)

# Print results
sep_line <- paste(rep("=", 60), collapse = "")
cat("\n", sep_line, "\n", "SIMULATION RESULTS:\n", sep_line, "\n", sep = "")
print(cmp)
cat("\nKS Test:\n")
cat("D-statistic:", round(ks_res$statistic, 4), "\n")
cat("p-value:", round(ks_res$p.value, 4), "\n")
cat("Fit:", ifelse(ks_res$p.value > 0.05, "Good (fail to reject H0)", "Poor (reject H0)"), "\n")

cat("\nFiles saved:\n- comparison_table.csv\n- 6_simulation_dashboard.png\n")