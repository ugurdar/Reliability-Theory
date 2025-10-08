# Week 1: Reliability Theory Analysis Solutions
# Hazard Rate: z(t) = h(t) = 0.4/(0.2*t + 1)

# Core functions
hazard_rate <- function(t) 0.4 / (0.2 * t + 1)
cumulative_hazard <- function(t) 2 * log(0.2 * t + 1)          # Z(t)
reliability <- function(t) 1 / ((0.2 * t + 1)^2)               # R(t)
pdf_func <- function(t) 0.4 / ((0.2 * t + 1)^3)                # f(t) = z(t)R(t)
cdf_func <- function(t) 1 - 1 / ((0.2 * t + 1)^2)              # F(t) = 1 - R(t)

# Metrics (analytical where possible)
calculate_mttf <- function() 5.0                                # ∫ R(t) dt
calculate_variance <- function() Inf                            # Lomax(α=2): Var = ∞
calculate_std_dev <- function() Inf
find_median <- function() 5 * (sqrt(2) - 1)                     # R(t)=0.5
find_mode <- function() 0.0                                     # f'(t)<0 ⇒ mode at 0
calculate_cv <- function() Inf
mean_residual_life <- function(t) t + 5                         # m(t) = t + 5
analyze_hazard_type <- function() "Decreasing Hazard Rate (DHR) - Negative Aging"

# Compute and print
cat("=== Week 1: Reliability Theory Analysis ===\n")
cat("Hazard Rate Function: z(t) = 0.4/(0.2t + 1)\n\n")

mttf <- calculate_mttf()
variance <- calculate_variance()
std_dev <- calculate_std_dev()
median_life <- find_median()
mode_life <- find_mode()
cv <- calculate_cv()

cat("ANALYTICAL SOLUTIONS:\n")
cat("====================\n")
cat("MTTF:", mttf, "\n")
cat("Variance:", variance, "\n")
cat("Standard Deviation:", std_dev, "\n")
cat("Median:", round(median_life, 4), "\n")
cat("Mode:", mode_life, "\n")
cat("Coefficient of Variation:", cv, "\n")
cat("Hazard Type:", analyze_hazard_type(), "\n\n")
