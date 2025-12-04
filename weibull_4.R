################################################################################
# ÇOKLU THRESHOLD + ÇOKLU MOMENT MRL
# Her threshold için hem m1 hem m2 kullan
################################################################################

rm(list = ls())

Gamma_upper <- function(s, z) {
  if (z < 0 || is.na(z)) return(gamma(s))
  if (z > 700) return(0)
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

calc_mrl_moments <- function(y0, beta, theta, gamma) {
  x0 <- y0 - gamma
  if (x0 <= 0) return(c(NA, NA, NA))
  
  t <- (x0 / theta)^beta
  if (t > 700) return(c(NA, NA, NA))
  
  exp_t <- exp(t)
  G1 <- Gamma_upper(1 + 1/beta, t)
  G2 <- Gamma_upper(1 + 2/beta, t)
  G3 <- Gamma_upper(1 + 3/beta, t)
  
  EX1 <- theta * exp_t * G1
  EX2 <- theta^2 * exp_t * G2
  EX3 <- theta^3 * exp_t * G3
  
  m1 <- EX1 - x0
  m2 <- EX2 - 2*x0*EX1 + x0^2
  m3 <- EX3 - 3*x0*EX2 + 3*x0^2*EX1 - x0^3
  
  c(m1, m2, m3)
}

# ==============================================================================
# YAKLAŞIM 1: 3 THRESHOLD, HER BİRİNDEN m1
# ==============================================================================

weibull_3p_mrl_3thr_m1 <- function(y, probs = c(0.15, 0.35, 0.55)) {
  n <- length(y)
  y_min <- min(y)
  y0_vec <- quantile(y, probs = probs, names = FALSE)
  
  m1_emp <- sapply(y0_vec, function(y0) {
    exc <- y[y > y0] - y0
    if (length(exc) < 5) return(NA)
    mean(exc)
  })
  
  if (any(is.na(m1_emp))) return(list(beta = NA, theta = NA, gamma = NA))
  
  objective <- function(par) {
    theta <- par[1]; beta <- par[2]; gamma <- par[3]
    if (theta <= 0 || beta <= 0 || gamma >= min(y0_vec)) return(1e10)
    
    err <- 0
    for (i in seq_along(y0_vec)) {
      m <- calc_mrl_moments(y0_vec[i], beta, theta, gamma)
      if (is.na(m[1])) return(1e10)
      err <- err + ((m[1] - m1_emp[i]) / m1_emp[i])^2
    }
    err
  }
  
  mom <- tryCatch(weibull_3p_mom(y), error = function(e) NULL)
  start <- if (!is.null(mom)) c(mom$theta, mom$beta, mom$gamma) else c(sd(y), 1.5, y_min*0.5)
  
  res <- tryCatch({
    optim(start, objective, method = "L-BFGS-B",
          lower = c(1e-6, 0.1, -Inf), upper = c(Inf, 20, min(y0_vec) - 1e-6),
          control = list(maxit = 2000))
  }, error = function(e) NULL)
  
  if (is.null(res)) return(list(beta = NA, theta = NA, gamma = NA))
  list(theta = res$par[1], beta = res$par[2], gamma = res$par[3], obj = res$value)
}

# ==============================================================================
# YAKLAŞIM 2: 2 THRESHOLD, HER BİRİNDEN m1 + m2 (4 denklem, 3 bilinmeyen)
# ==============================================================================

weibull_3p_mrl_2thr_m1m2 <- function(y, probs = c(0.20, 0.50)) {
  n <- length(y)
  y_min <- min(y)
  y0_vec <- quantile(y, probs = probs, names = FALSE)
  
  emp <- lapply(y0_vec, function(y0) {
    exc <- y[y > y0] - y0
    if (length(exc) < 5) return(c(NA, NA))
    c(mean(exc), mean(exc^2))
  })
  
  if (any(sapply(emp, function(x) any(is.na(x))))) {
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  objective <- function(par) {
    theta <- par[1]; beta <- par[2]; gamma <- par[3]
    if (theta <= 0 || beta <= 0 || gamma >= min(y0_vec)) return(1e10)
    
    err <- 0
    for (i in seq_along(y0_vec)) {
      m <- calc_mrl_moments(y0_vec[i], beta, theta, gamma)
      if (any(is.na(m[1:2]))) return(1e10)
      
      err <- err + ((m[1] - emp[[i]][1]) / emp[[i]][1])^2
      err <- err + ((m[2] - emp[[i]][2]) / emp[[i]][2])^2
    }
    err
  }
  
  mom <- tryCatch(weibull_3p_mom(y), error = function(e) NULL)
  
  starts <- list(
    if (!is.null(mom)) c(mom$theta, mom$beta, mom$gamma) else c(sd(y), 1.5, y_min*0.5),
    c(sd(y), 1.5, y_min * 0.3),
    c(sd(y) * 0.8, 2.0, 0),
    c(sd(y) * 1.2, 1.0, y_min * 0.7)
  )
  
  best <- NULL; best_val <- Inf
  for (start in starts) {
    res <- tryCatch({
      optim(start, objective, method = "L-BFGS-B",
            lower = c(1e-6, 0.1, -Inf), upper = c(Inf, 20, min(y0_vec) - 1e-6),
            control = list(maxit = 2000))
    }, error = function(e) NULL)
    if (!is.null(res) && res$value < best_val) {
      best_val <- res$value; best <- res
    }
  }
  
  if (is.null(best)) return(list(beta = NA, theta = NA, gamma = NA))
  list(theta = best$par[1], beta = best$par[2], gamma = best$par[3], obj = best_val)
}

# ==============================================================================
# YAKLAŞIM 3: TEK THRESHOLD + MOM GAMMA (hibrit)
# ==============================================================================

weibull_3p_mrl_hybrid <- function(y, threshold_prob = 0.25) {
  n <- length(y)
  y_min <- min(y)
  
  # MOM'dan gamma al
  mom <- weibull_3p_mom(y)
  gamma_fixed <- mom$gamma
  if (gamma_fixed >= y_min) gamma_fixed <- y_min * 0.9
  
  y0 <- quantile(y, probs = threshold_prob, names = FALSE)
  exc <- y[y > y0] - y0
  if (length(exc) < 5) return(list(beta = NA, theta = NA, gamma = NA))
  
  m1_emp <- mean(exc)
  m2_emp <- mean(exc^2)
  
  # Sadece theta ve beta optimize et
  objective <- function(par) {
    theta <- par[1]; beta <- par[2]
    if (theta <= 0 || beta <= 0) return(1e10)
    
    m <- calc_mrl_moments(y0, beta, theta, gamma_fixed)
    if (any(is.na(m[1:2]))) return(1e10)
    
    ((m[1] - m1_emp) / m1_emp)^2 + ((m[2] - m2_emp) / m2_emp)^2
  }
  
  res <- tryCatch({
    optim(c(mom$theta, mom$beta), objective, method = "L-BFGS-B",
          lower = c(1e-6, 0.1), upper = c(Inf, 20))
  }, error = function(e) NULL)
  
  if (is.null(res)) return(list(beta = NA, theta = NA, gamma = NA))
  list(theta = res$par[1], beta = res$par[2], gamma = gamma_fixed, obj = res$value)
}

# ==============================================================================
# MOM ve MLE
# ==============================================================================

weibull_3p_mom <- function(y) {
  n <- length(y)
  y_bar <- mean(y)
  s2 <- var(y) * (n-1) / n
  m2c <- mean((y - y_bar)^2)
  m3c <- mean((y - y_bar)^3)
  skew <- m3c / (m2c^(3/2))
  
  weibull_skew <- function(b) {
    g1 <- gamma(1 + 1/b); g2 <- gamma(1 + 2/b); g3 <- gamma(1 + 3/b)
    (g3 - 3*g1*g2 + 2*g1^3) / (g2 - g1^2)^(3/2)
  }
  
  opt <- optimize(function(b) (weibull_skew(b) - skew)^2, c(0.1, 20))
  beta_hat <- opt$minimum
  g1 <- gamma(1 + 1/beta_hat); g2 <- gamma(1 + 2/beta_hat)
  theta_hat <- sqrt(s2 / (g2 - g1^2))
  gamma_hat <- y_bar - theta_hat * g1
  
  list(beta = beta_hat, theta = theta_hat, gamma = gamma_hat)
}

weibull_3p_mle <- function(y) {
  n <- length(y); y_min <- min(y)
  
  nll <- function(par) {
    theta <- par[1]; beta <- par[2]; gamma <- par[3]
    if (theta <= 0 || beta <= 0 || gamma >= y_min) return(Inf)
    x <- y - gamma
    if (any(x <= 0)) return(Inf)
    -sum(dweibull(x, shape = beta, scale = theta, log = TRUE))
  }
  
  mom <- tryCatch(weibull_3p_mom(y), error = function(e) NULL)
  start <- if (!is.null(mom)) c(mom$theta, mom$beta, mom$gamma) else c(sd(y), 1.5, y_min*0.5)
  
  res <- tryCatch({
    optim(start, nll, method = "L-BFGS-B",
          lower = c(1e-6, 1e-6, -Inf), upper = c(Inf, Inf, y_min - 1e-6))
  }, error = function(e) NULL)
  
  if (is.null(res)) return(list(beta = NA, theta = NA, gamma = NA))
  list(theta = res$par[1], beta = res$par[2], gamma = res$par[3])
}

# ==============================================================================
# SİMÜLASYON
# ==============================================================================

cat("=" , rep("=", 70), "\n", sep="")
cat("SİMÜLASYON: FARKLI MRL YAKLAŞIMLARI (n=200, 500 tekrar)\n")
cat("=" , rep("=", 70), "\n\n", sep="")

set.seed(123)
beta_true <- 1.5; theta_true <- 1.0; gamma_true <- 0.5
n_sim <- 500; n <- 200

mom_res <- mle_res <- mrl_3thr_res <- mrl_2thr_res <- mrl_hybrid_res <- matrix(NA, n_sim, 3)

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (i in 1:n_sim) {
  X <- rweibull(n, shape = beta_true, scale = theta_true)
  Y <- gamma_true + X
  
  tryCatch({ fit <- weibull_3p_mom(Y); mom_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  tryCatch({ fit <- weibull_3p_mle(Y); if (!is.na(fit$beta)) mle_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  tryCatch({ fit <- weibull_3p_mrl_3thr_m1(Y); if (!is.na(fit$beta)) mrl_3thr_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  tryCatch({ fit <- weibull_3p_mrl_2thr_m1m2(Y); if (!is.na(fit$beta)) mrl_2thr_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  tryCatch({ fit <- weibull_3p_mrl_hybrid(Y); if (!is.na(fit$beta)) mrl_hybrid_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  
  setTxtProgressBar(pb, i)
}
close(pb)

calc_stats <- function(res, name, true_vals) {
  valid <- complete.cases(res)
  n_valid <- sum(valid)
  if (n_valid < 10) return(data.frame(Method = name, beta_mean = NA, beta_mse = NA,
                                      theta_mean = NA, theta_mse = NA,
                                      gamma_mean = NA, gamma_mse = NA, n_valid = n_valid))
  res_clean <- res[valid, ]
  means <- colMeans(res_clean)
  vars <- apply(res_clean, 2, var)
  bias <- means - true_vals
  mse <- bias^2 + vars
  
  data.frame(Method = name,
             beta_mean = means[1], beta_bias = bias[1], beta_mse = mse[1],
             theta_mean = means[2], theta_bias = bias[2], theta_mse = mse[2],
             gamma_mean = means[3], gamma_bias = bias[3], gamma_mse = mse[3],
             n_valid = n_valid)
}

true_vals <- c(beta_true, theta_true, gamma_true)

results <- rbind(
  calc_stats(mom_res, "MOM", true_vals),
  calc_stats(mle_res, "MLE", true_vals),
  calc_stats(mrl_3thr_res, "MRL_3thr_m1", true_vals),
  calc_stats(mrl_2thr_res, "MRL_2thr_m1m2", true_vals),
  calc_stats(mrl_hybrid_res, "MRL_hybrid", true_vals)
)

cat("\n\n")
cat("Sonuçlar:\n")
cat("-" , rep("-", 90), "\n", sep="")
print(results[, c("Method", "beta_mean", "beta_mse", "theta_mean", "theta_mse", 
                  "gamma_mean", "gamma_mse", "n_valid")], digits = 4)

cat("\nGerçek: beta =", beta_true, ", theta =", theta_true, ", gamma =", gamma_true, "\n")

# Toplam MSE
results$total_mse <- results$beta_mse + results$theta_mse + results$gamma_mse
cat("\n\nToplam MSE sıralaması:\n")
ord <- order(results$total_mse)
for (i in ord) {
  cat("  ", results$Method[i], ": ", round(results$total_mse[i], 5), "\n", sep="")
}

# ==============================================================================
# BÜYÜK ÖRNEKLEM
# ==============================================================================

cat("\n\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("BÜYÜK ÖRNEKLEM: n=500, 300 tekrar\n")
cat("=" , rep("=", 70), "\n\n", sep="")

set.seed(456)
n_sim <- 300; n <- 500

mom_res <- mle_res <- mrl_2thr_res <- mrl_hybrid_res <- matrix(NA, n_sim, 3)

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (i in 1:n_sim) {
  X <- rweibull(n, shape = beta_true, scale = theta_true)
  Y <- gamma_true + X
  
  tryCatch({ fit <- weibull_3p_mom(Y); mom_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  tryCatch({ fit <- weibull_3p_mle(Y); if (!is.na(fit$beta)) mle_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  tryCatch({ fit <- weibull_3p_mrl_2thr_m1m2(Y); if (!is.na(fit$beta)) mrl_2thr_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  tryCatch({ fit <- weibull_3p_mrl_hybrid(Y); if (!is.na(fit$beta)) mrl_hybrid_res[i, ] <- c(fit$beta, fit$theta, fit$gamma) }, error = function(e) NULL)
  
  setTxtProgressBar(pb, i)
}
close(pb)

results2 <- rbind(
  calc_stats(mom_res, "MOM", true_vals),
  calc_stats(mle_res, "MLE", true_vals),
  calc_stats(mrl_2thr_res, "MRL_2thr_m1m2", true_vals),
  calc_stats(mrl_hybrid_res, "MRL_hybrid", true_vals)
)

cat("\n\n")
print(results2[, c("Method", "beta_mean", "beta_mse", "theta_mean", "theta_mse",
                   "gamma_mean", "gamma_mse", "n_valid")], digits = 4)

# ==============================================================================
# MAKALE VERİSİ
# ==============================================================================

cat("\n\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("MAKALE VERİSİ (NA_a)\n")
cat("=" , rep("=", 70), "\n\n", sep="")

dataset1 <- c(24.12, 24.13, 28.52, 29.18, 29.67, 30.48, 32.98, 35.91, 35.92,
              36.38, 37.60, 37.70, 39.71, 49.10, 52.43, 52.46, 52.61, 61.72)

cat("n =", length(dataset1), ", mean =", round(mean(dataset1), 2), "\n\n")

cat("MOM:\n"); fit <- weibull_3p_mom(dataset1)
cat("  beta =", round(fit$beta, 4), ", theta =", round(fit$theta, 4), ", gamma =", round(fit$gamma, 4), "\n\n")

cat("MLE:\n"); fit <- weibull_3p_mle(dataset1)
cat("  beta =", round(fit$beta, 4), ", theta =", round(fit$theta, 4), ", gamma =", round(fit$gamma, 4), "\n\n")

cat("MRL_2thr_m1m2:\n"); fit <- weibull_3p_mrl_2thr_m1m2(dataset1)
cat("  beta =", round(fit$beta, 4), ", theta =", round(fit$theta, 4), ", gamma =", round(fit$gamma, 4), "\n\n")

cat("MRL_hybrid:\n"); fit <- weibull_3p_mrl_hybrid(dataset1)
cat("  beta =", round(fit$beta, 4), ", theta =", round(fit$theta, 4), ", gamma =", round(fit$gamma, 4), "\n\n")

cat("Makaledeki PITE:\n")
cat("  rho=0.32: beta=2.399, theta=37.961, gamma=6.976\n")
cat("  rho=1.00: beta=1.967, theta=39.524, gamma=4.746\n")