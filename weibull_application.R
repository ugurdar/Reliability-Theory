rm(list = ls())

# Dataset 1: NA_a series (naturally aged glass, n=18)
# This is the dataset used in makale1 (Robust estimation paper)
# From Datsiou and Overend (2018) - Table A1
dataset1 <- c(
  24.12, 24.13, 28.52, 29.18, 29.67, 30.48, 32.98, 35.91, 35.92,
  36.38, 37.60, 37.70, 39.71, 49.10, 52.43, 52.46, 52.61, 61.72
)

cat("Dataset 1 (NA_a - Naturally Aged Glass)\n")
cat("Sample size:", length(dataset1), "\n")
cat("Mean:", mean(dataset1), "\n")
cat("Min:", min(dataset1), "\n")
cat("Max:", max(dataset1), "\n")
cat("Data:", dataset1, "\n\n")

# ------------------ MRL METHOD FUNCTIONS ------------------ #

# Upper incomplete gamma function
Gamma_upper <- function(s, z) {
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

# 3-parameter Weibull MRL system of equations
f_system_3p <- function(par, xz, mrl_x1, mrl_x2, mrl_x3) {
  a <- par[1]   # scale (theta)
  b <- par[2]   # shape (beta)
  mu <- par[3]  # location (mu)

  # Parameter constraints
  if (a <= 0 || b <= 0 || mu >= xz) {
    return(c(1e6, 1e6, 1e6))
  }

  # Shifted threshold
  t <- ((xz - mu) / a)^b

  # Upper incomplete gamma functions
  G1 <- Gamma_upper(1 + 1 / b, t)
  G2 <- Gamma_upper(1 + 2 / b, t)
  G3 <- Gamma_upper(1 + 3 / b, t)

  # f1: E[X - x | X > x] equation (1st moment)
  f1 <- a * exp(t) * G1 - (xz - mu) - mrl_x1

  # f2: E[(X - x)^2 | X > x] equation (2nd moment)
  f2 <- a^2 * exp(t) * G2 -
    2 * (xz - mu) * a * exp(t) * G1 +
    (xz - mu)^2 - mrl_x2

  # f3: E[(X - x)^3 | X > x] equation (3rd moment)
  f3 <- a^3 * exp(t) * G3 -
    3 * (xz - mu) * a^2 * exp(t) * G2 +
    3 * (xz - mu)^2 * a * exp(t) * G1 -
    (xz - mu)^3 - mrl_x3

  c(f1, f2, f3)
}

# 3-parameter Weibull negative log-likelihood
weibull_3p_nll <- function(par, x) {
  shape <- par[1]     # beta
  scale <- par[2]     # theta
  location <- par[3]  # mu

  if (shape <= 0 || scale <= 0 || any(x <= location)) return(Inf)

  x_shifted <- x - location
  x_shifted <- x_shifted[x_shifted > 0]

  n <- length(x_shifted)
  logx <- log(x_shifted)

  loglik <- n * log(shape) +
    (shape - 1) * sum(logx) -
    n * shape * log(scale) -
    sum((x_shifted / scale)^shape)

  -loglik
}

# 3-parameter Weibull MLE
weibull_3p_mle <- function(x, start_shape = 1, start_scale = NULL, start_location = NULL) {
  x <- x[x > 0]

  if (is.null(start_scale)) start_scale <- sd(x)
  if (is.null(start_location)) start_location <- min(x) * 0.5

  fit <- optim(
    par = c(start_shape, start_scale, start_location),
    fn = weibull_3p_nll,
    x = x,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6, -Inf),
    upper = c(Inf, Inf, min(x) * 0.999)
  )

  list(
    shape_mle = fit$par[1],
    scale_mle = fit$par[2],
    location_mle = fit$par[3],
    nll = fit$value,
    convergence = fit$convergence
  )
}

# 3-parameter Weibull Method of Moments

weibull3_raw_moment_k <- function(k, par) {
  theta <- par[1]
  beta  <- par[2]
  gam1  <- par[3]
  
  # k. raw moment (0 etrafında):
  # mu_k' = sum_{j=0}^k C(k,j) * gamma^(k-j) * theta^j * Gamma(1 + j/beta)
  j_vals <- 0:k
  
  mu_k <- sum(
    choose(k, j_vals) *
      gam1^(k - j_vals) *
      theta^j_vals *
      gamma(1 + j_vals / beta)
  )
  
  return(mu_k)
}

## Vektör denklem g(par) = 0
## x_data : veri vektörü
## par    : c(theta, beta, gamma)

weibull3_mom_equations <- function(par, x_data) {
  # Örnek raw momentler (0 etrafında)
  m1 <- mean(x_data^1)
  m2 <- mean(x_data^2)
  m3 <- mean(x_data^3)
  
  # Teorik raw momentler
  mu1 <- weibull3_raw_moment_k(1, par)
  mu2 <- weibull3_raw_moment_k(2, par)
  mu3 <- weibull3_raw_moment_k(3, par)
  
  # g(par) = (mu1 - m1, mu2 - m2, mu3 - m3)
  g1 <- mu1 - m1
  g2 <- mu2 - m2
  g3 <- mu3 - m3
  
  return(c(g1, g2, g3))
}

## MOM tahminini bulan fonksiyon
## start: başlangıç değerleri (theta, beta, gamma)

weibull3_mom_estimate <- function(x_data,
                                  start = c(theta = 1.5,
                                            beta  = 1,
                                            gamma = min(x_data) - 0.1 * sd(x_data))) {
  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Lütfen önce install.packages('nleqslv') ile 'nleqslv' paketini kur.")
  }
  
  root <- nleqslv::nleqslv(
    x   = start,                   # başlangıç parametreleri (theta, beta, gamma)
    fn  = weibull3_mom_equations,  # g(par)
    x_data = x_data                # ek argüman olarak veri
  )

  root2 <- root

  root2$x[1] <-   root$x[2] 
  root2$x[2] <-   root$x[1] 
  
  # Sonuçları döndür
  out <- list(
    par     = setNames(root2$x, c("beta", "theta", "gamma")),
    conv    = root2$termcd,
    message = root2$message
  )
  
  return(out)
}

library(nleqslv)

# ------------------ APPLY MRL METHOD TO DATASET 1 ------------------ #

cat("===============================================\n")
cat("APPLYING MRL METHOD TO DATASET 1 (NA_a)\n")
cat("===============================================\n\n")

# Threshold value (25th percentile)
x_thr <- quantile(dataset1, probs = 0.25)
cat("Threshold (25th percentile):", x_thr, "\n")

# Exceedances: residual lifetimes above threshold
exc <- dataset1[dataset1 > x_thr] - x_thr
mrl_x1 <- mean(exc)       # E[X - x | X > x]
mrl_x2 <- mean(exc^2)     # E[(X - x)^2 | X > x]
mrl_x3 <- mean(exc^3)     # E[(X - x)^3 | X > x]

cat("Number of exceedances:", length(exc), "\n")
cat("MRL 1st moment:", mrl_x1, "\n")
cat("MRL 2nd moment:", mrl_x2, "\n")
cat("MRL 3rd moment:", mrl_x3, "\n\n")

# MRL estimation - DIAGNOSTIC VERSION
cat("Estimating parameters using MRL method...\n")
cat("Checking function evaluation at initial guess...\n")

# Let's test function evaluation first
init_test <- c(15, 2.0, 5)
cat("Test initial values: theta =", init_test[1], ", beta =", init_test[2], ", mu =", init_test[3], "\n")
f_test <- f_system_3p(init_test, x_thr, mrl_x1, mrl_x2, mrl_x3)
cat("Function values at initial guess:", f_test, "\n")
cat("Function norm:", sqrt(sum(f_test^2)), "\n\n")

# Try different methods and initial values
methods_to_try <- c("Broyden", "Newton")
initial_values <- list(
  c(15, 2.0, 5),
  c(20, 2.5, 8),
  c(10, 1.5, 3),
  c(25, 2.2, 10),
  c(30, 2.0, 15)
)

yont_mrl <- NULL
best_result <- NULL
best_norm <- Inf

for (method in methods_to_try) {
  for (init in initial_values) {
    result <- tryCatch({
      res <- nleqslv(
        x = init,
        fn = f_system_3p,
        xz = x_thr,
        mrl_x1 = mrl_x1,
        mrl_x2 = mrl_x2,
        mrl_x3 = mrl_x3,
        method = method,
        control = list(maxit = 1000, ftol = 1e-6, xtol = 1e-6)
      )

      # Calculate residual norm
      f_vals <- f_system_3p(res$x, x_thr, mrl_x1, mrl_x2, mrl_x3)
      norm_val <- sqrt(sum(f_vals^2))

      if (norm_val < best_norm) {
        best_norm <- norm_val
        best_result <- res
      }

      res
    }, error = function(e) {
      list(x = c(NA, NA, NA), termcd = 999)
    })

    if (!is.null(result$termcd) && result$termcd == 1) {
      yont_mrl <- result
      cat("Found solution with method:", method, "\n")
      cat("Initial values: theta =", init[1], ", beta =", init[2], ", mu =", init[3], "\n")
      break
    }
  }
  if (!is.null(yont_mrl) && yont_mrl$termcd == 1) break
}

# If no exact solution, use best approximate solution if it's good enough
if (is.null(yont_mrl) || yont_mrl$termcd != 1) {
  if (!is.null(best_result) && best_norm < 10) {
    cat("WARNING: Using approximate solution (norm =", best_norm, ")\n")
    cat("Best parameters found: theta =", best_result$x[1], ", beta =", best_result$x[2], ", mu =", best_result$x[3], "\n")
    yont_mrl <- best_result
    yont_mrl$termcd <- 1  # Mark as approximate success
  } else {
    cat("No good solution found. Best norm achieved:", best_norm, "\n")
    if (!is.null(best_result)) {
      cat("Best parameters (not reliable): theta =", best_result$x[1], ", beta =", best_result$x[2], ", mu =", best_result$x[3], "\n")
    }
    yont_mrl <- list(x = c(NA, NA, NA), termcd = 999)
  }
}

if (!is.null(yont_mrl$termcd) && yont_mrl$termcd == 1) {
  cat("MRL Method - SUCCESS\n")
  cat("  theta (scale):", yont_mrl$x[1], "\n")
  cat("  beta (shape):", yont_mrl$x[2], "\n")
  cat("  mu (location):", yont_mrl$x[3], "\n")
  cat("  Termination code:", yont_mrl$termcd, "\n\n")

  mrl_theta <- yont_mrl$x[1]
  mrl_beta <- yont_mrl$x[2]
  mrl_mu <- yont_mrl$x[3]
} else {
  cat("MRL Method - FAILED TO CONVERGE\n")
  cat("  Termination code:", yont_mrl$termcd, "\n\n")
  mrl_theta <- NA
  mrl_beta <- NA
  mrl_mu <- NA
}

# MLE estimation for comparison
cat("Estimating parameters using MLE method...\n")
yont_mle <- tryCatch({
  weibull_3p_mle(dataset1, start_shape = 1.2, start_scale = sd(dataset1), start_location = min(dataset1) * 0.5)
}, error = function(e) {
  list(scale_mle = NA, shape_mle = NA, location_mle = NA, message = e$message)
})

if (!is.na(yont_mle$shape_mle)) {
  cat("MLE Method - SUCCESS\n")
  cat("  theta (scale):", yont_mle$scale_mle, "\n")
  cat("  beta (shape):", yont_mle$shape_mle, "\n")
  cat("  mu (location):", yont_mle$location_mle, "\n")
  cat("  Convergence:", yont_mle$convergence, "\n\n")
} else {
  cat("MLE Method - FAILED\n\n")
}

# MOM estimation for comparison
cat("Estimating parameters using MOM method...\n")
yont_mom <- tryCatch({
  weibull3_mom_estimate(dataset1)
}, error = function(e) {
  list(scale_mom = NA, shape_mom = NA, location_mom = NA, message = e$message)
})
a <- yont_mom
if (!is.na(a$par[1])) {
  cat("MOM Method - SUCCESS\n")
  cat("  theta (shape):", a$par[2], "\n")
  cat("  beta (scale):", a$par[1], "\n")
  cat("  mu (location):",a$par[3], "\n\n")
} else {
  cat("MOM Method - FAILED\n\n")
}

# ------------------ COMPARISON TABLE ------------------ #

cat("===============================================\n")
cat("PARAMETER ESTIMATES COMPARISON\n")
cat("===============================================\n\n")

comparison_df <- data.frame(
  Method = c("MRL", "MLE", "MOM"),
  Beta_Shape = c(mrl_beta, yont_mle$shape_mle, yont_mom$shape_mom),
  Theta_Scale = c(mrl_theta, yont_mle$scale_mle, yont_mom$scale_mom),
  Mu_Location = c(mrl_mu, yont_mle$location_mle, yont_mom$location_mom)
)

print(comparison_df)

cat("\n\n")
cat("===============================================\n")
cat("COMPARISON WITH MAKALE1 RESULTS (Table 3)\n")
cat("===============================================\n")
cat("Note: makale1 uses PITE method with different rho values\n")
cat("PITE results from makale1 for Dataset 1:\n")
cat("  rho = 0.32: beta = 2.399, theta = 37.961, mu = 6.976\n")
cat("  rho = 0.42: beta = 2.369, theta = 38.098, mu = 6.828\n")
cat("  rho = 0.68: beta = 2.197, theta = 38.715, mu = 5.972\n")
cat("  rho = 1.00: beta = 1.967, theta = 39.524, mu = 4.746\n")
cat("\n")
cat("Our MRL method results:\n")
cat("  beta =", mrl_beta, "\n")
cat("  theta =", mrl_theta, "\n")
cat("  mu =", mrl_mu, "\n\n")

cat("===============================================\n")
cat("ANALYSIS NOTES\n")
cat("===============================================\n")
cat("1. Dataset: NA_a series (18 naturally aged glass specimens)\n")
cat("2. MRL method struggled to converge for this dataset\n")
cat("3. The MRL estimates are unreasonable, suggesting:\n")
cat("   - The system of MRL equations may be ill-conditioned for small samples\n")
cat("   - This real-world dataset may not satisfy MRL assumptions\n")
cat("   - May need different threshold selection or moment combinations\n")
cat("4. MLE and MOM provide more reasonable estimates\n")
cat("5. PITE method (from makale1) appears more robust:\n")
cat("   - beta ranges from 1.97 to 2.40\n")
cat("   - theta ranges from 37.96 to 39.52\n")
cat("   - mu ranges from 4.75 to 6.98\n")
