################################################################################
# GELİŞTİRİLMİŞ MRL YÖNTEMİ - ÇOKLU THRESHOLD
# 
# Problem: Tek threshold ile 3 moment denklemi ill-conditioned
# Çözüm: Birden fazla threshold kullan, her birinden sadece 1. moment (MRL) al
################################################################################

rm(list = ls())

# Üst eksik gamma fonksiyonu
gamma_upper <- function(s, z) {
  if (z < 0) return(gamma(s))
  if (z > 700) return(0)
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

# ==============================================================================
# YENİ YAKLAŞIM 1: ÇOKLU THRESHOLD, SADECE 1. MOMENT
# ==============================================================================
# 
# Her threshold y0_i için:
#   m(y0_i) = E[Y - y0_i | Y > y0_i] = theta * exp(t_i) * Gamma(1+1/beta, t_i) - z0_i
# 
# 3 farklı threshold → 3 denklem, daha iyi koşullu olabilir
# ==============================================================================

# Teorik MRL (sadece 1. moment)
calc_mrl_m1 <- function(y0, beta, theta, gamma) {
  z0 <- y0 - gamma
  if (z0 <= 0) return(NA)
  
  t <- (z0 / theta)^beta
  if (t > 700) return(NA)
  
  exp_t <- exp(t)
  G1 <- gamma_upper(1 + 1/beta, t)
  
  theta * exp_t * G1 - z0
}

# Çoklu threshold MRL yöntemi
weibull_3p_mrl_multi <- function(y, threshold_probs = c(0.10, 0.30, 0.50)) {
  n <- length(y)
  y_min <- min(y)
  
  # Threshold değerlerini hesapla
  y0_vec <- quantile(y, probs = threshold_probs, names = FALSE)
  
  # Her threshold için empirik MRL
  mrl_emp <- sapply(y0_vec, function(y0) {
    exc <- y[y > y0] - y0
    if (length(exc) < 5) return(NA)
    mean(exc)
  })
  
  if (any(is.na(mrl_emp))) {
    warning("Bazı threshold'larda yetersiz gözlem")
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  # Objective: Her threshold için (teorik - empirik)^2 / empirik^2
  objective <- function(par) {
    theta <- par[1]
    beta <- par[2]
    gamma <- par[3]
    
    if (theta <= 0 || beta <= 0 || gamma >= min(y0_vec)) {
      return(1e10)
    }
    
    err_sum <- 0
    for (i in seq_along(y0_vec)) {
      m1_theo <- calc_mrl_m1(y0_vec[i], beta, theta, gamma)
      if (is.na(m1_theo)) return(1e10)
      
      err_sum <- err_sum + ((m1_theo - mrl_emp[i]) / mrl_emp[i])^2
    }
    
    err_sum
  }
  
  # Başlangıç değerleri
  mom <- tryCatch({
    y_bar <- mean(y)
    s2 <- var(y) * (n-1) / n
    m2 <- mean((y - y_bar)^2)
    m3 <- mean((y - y_bar)^3)
    skew <- m3 / (m2^(3/2))
    
    weibull_skew <- function(b) {
      g1 <- gamma(1 + 1/b)
      g2 <- gamma(1 + 2/b)
      g3 <- gamma(1 + 3/b)
      (g3 - 3*g1*g2 + 2*g1^3) / (g2 - g1^2)^(3/2)
    }
    
    opt <- optimize(function(b) (weibull_skew(b) - skew)^2, c(0.5, 10))
    beta_hat <- opt$minimum
    
    g1 <- gamma(1 + 1/beta_hat)
    g2 <- gamma(1 + 2/beta_hat)
    theta_hat <- sqrt(s2 / (g2 - g1^2))
    gamma_hat <- y_bar - theta_hat * g1
    
    c(theta_hat, beta_hat, gamma_hat)
  }, error = function(e) c(sd(y), 1.5, y_min * 0.5))
  
  starts <- list(
    mom,
    c(mom[1] * 0.8, mom[2] * 1.2, mom[3]),
    c(mom[1] * 1.2, mom[2] * 0.8, mom[3]),
    c(mom[1], mom[2], mom[3] * 0.5),
    c(mom[1], mom[2], 0)
  )
  
  best <- NULL
  best_val <- Inf
  
  for (start in starts) {
    res <- tryCatch({
      optim(
        par = start,
        fn = objective,
        method = "L-BFGS-B",
        lower = c(1e-6, 0.3, -Inf),
        upper = c(Inf, 15, min(y0_vec) - 1e-6),
        control = list(maxit = 2000)
      )
    }, error = function(e) NULL)
    
    if (!is.null(res) && res$value < best_val) {
      best_val <- res$value
      best <- res
    }
  }
  
  if (is.null(best)) {
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  list(
    theta = best$par[1],
    beta = best$par[2],
    gamma = best$par[3],
    objective = best_val,
    thresholds = threshold_probs
  )
}

# ==============================================================================
# YENİ YAKLAŞIM 2: İKİ AŞAMALI TAHMİN
# ==============================================================================
#
# Aşama 1: gamma'yı MOM'dan al veya sabit tut
# Aşama 2: Sadece theta ve beta'yı MRL ile tahmin et (2 denklem, 2 bilinmeyen)
# ==============================================================================

weibull_3p_mrl_twostage <- function(y, threshold_prob = 0.25) {
  n <- length(y)
  y_min <- min(y)
  
  # Aşama 1: gamma'yı MOM ile tahmin et
  y_bar <- mean(y)
  s2 <- var(y) * (n-1) / n
  m2 <- mean((y - y_bar)^2)
  m3 <- mean((y - y_bar)^3)
  skew <- m3 / (m2^(3/2))
  
  weibull_skew <- function(b) {
    g1 <- gamma(1 + 1/b)
    g2 <- gamma(1 + 2/b)
    g3 <- gamma(1 + 3/b)
    (g3 - 3*g1*g2 + 2*g1^3) / (g2 - g1^2)^(3/2)
  }
  
  opt <- tryCatch({
    optimize(function(b) (weibull_skew(b) - skew)^2, c(0.5, 10))
  }, error = function(e) list(minimum = 1.5))
  
  beta_init <- opt$minimum
  g1 <- gamma(1 + 1/beta_init)
  g2 <- gamma(1 + 2/beta_init)
  theta_init <- sqrt(s2 / (g2 - g1^2))
  gamma_fixed <- y_bar - theta_init * g1
  
  # gamma'nın y_min'den küçük olduğundan emin ol
  if (gamma_fixed >= y_min) {
    gamma_fixed <- y_min * 0.9
  }
  
  # Aşama 2: gamma sabit, theta ve beta'yı MRL ile tahmin et
  y0 <- quantile(y, probs = threshold_prob, names = FALSE)
  exc <- y[y > y0] - y0
  
  if (length(exc) < 10) {
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  m1_emp <- mean(exc)
  m2_emp <- mean(exc^2)
  
  objective_2p <- function(par) {
    theta <- par[1]
    beta <- par[2]
    
    if (theta <= 0 || beta <= 0) return(1e10)
    
    z0 <- y0 - gamma_fixed
    if (z0 <= 0) return(1e10)
    
    t <- (z0 / theta)^beta
    if (t > 700) return(1e10)
    
    exp_t <- exp(t)
    G1 <- gamma_upper(1 + 1/beta, t)
    G2 <- gamma_upper(1 + 2/beta, t)
    
    EZ_cond <- theta * exp_t * G1
    EZ2_cond <- theta^2 * exp_t * G2
    
    m1_theo <- EZ_cond - z0
    m2_theo <- EZ2_cond - 2*z0*EZ_cond + z0^2
    
    err1 <- ((m1_theo - m1_emp) / m1_emp)^2
    err2 <- ((m2_theo - m2_emp) / m2_emp)^2
    
    err1 + err2
  }
  
  starts <- list(
    c(theta_init, beta_init),
    c(theta_init * 0.8, beta_init * 1.2),
    c(theta_init * 1.2, beta_init * 0.8),
    c(sd(y), 1.5),
    c(sd(y) * 0.5, 2.0)
  )
  
  best <- NULL
  best_val <- Inf
  
  for (start in starts) {
    res <- tryCatch({
      optim(
        par = start,
        fn = objective_2p,
        method = "L-BFGS-B",
        lower = c(1e-6, 0.3),
        upper = c(Inf, 15),
        control = list(maxit = 2000)
      )
    }, error = function(e) NULL)
    
    if (!is.null(res) && res$value < best_val) {
      best_val <- res$value
      best <- res
    }
  }
  
  if (is.null(best)) {
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  list(
    theta = best$par[1],
    beta = best$par[2],
    gamma = gamma_fixed,
    objective = best_val,
    method = "two-stage"
  )
}

# ==============================================================================
# YENİ YAKLAŞIM 3: MRL ORANI (ratio-based)
# ==============================================================================
#
# İki farklı threshold için MRL oranı:
# m(y0_1) / m(y0_2) = [theta*exp(t1)*G1(t1) - z0_1] / [theta*exp(t2)*G1(t2) - z0_2]
#
# Bu oran theta'dan bağımsızlaşabilir, önce beta'yı buluruz
# ==============================================================================

weibull_3p_mrl_ratio <- function(y, threshold_probs = c(0.25, 0.50)) {
  n <- length(y)
  y_min <- min(y)
  
  y0_1 <- quantile(y, probs = threshold_probs[1], names = FALSE)
  y0_2 <- quantile(y, probs = threshold_probs[2], names = FALSE)
  
  exc1 <- y[y > y0_1] - y0_1
  exc2 <- y[y > y0_2] - y0_2
  
  if (length(exc1) < 10 || length(exc2) < 10) {
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  m1_emp <- mean(exc1)
  m2_emp <- mean(exc2)
  ratio_emp <- m1_emp / m2_emp
  
  # MOM'dan gamma tahmini
  y_bar <- mean(y)
  s2 <- var(y) * (n-1) / n
  m2c <- mean((y - y_bar)^2)
  m3c <- mean((y - y_bar)^3)
  skew <- m3c / (m2c^(3/2))
  
  weibull_skew <- function(b) {
    g1 <- gamma(1 + 1/b)
    g2 <- gamma(1 + 2/b)
    g3 <- gamma(1 + 3/b)
    (g3 - 3*g1*g2 + 2*g1^3) / (g2 - g1^2)^(3/2)
  }
  
  opt <- tryCatch({
    optimize(function(b) (weibull_skew(b) - skew)^2, c(0.5, 10))
  }, error = function(e) list(minimum = 1.5))
  
  beta_init <- opt$minimum
  g1 <- gamma(1 + 1/beta_init)
  g2 <- gamma(1 + 2/beta_init)
  theta_init <- sqrt(s2 / (g2 - g1^2))
  gamma_init <- y_bar - theta_init * g1
  
  if (gamma_init >= y_min) gamma_init <- y_min * 0.9
  
  # Objective
  objective <- function(par) {
    theta <- par[1]
    beta <- par[2]
    gamma <- par[3]
    
    if (theta <= 0 || beta <= 0 || gamma >= y0_1) return(1e10)
    
    m1_theo <- calc_mrl_m1(y0_1, beta, theta, gamma)
    m2_theo <- calc_mrl_m1(y0_2, beta, theta, gamma)
    
    if (is.na(m1_theo) || is.na(m2_theo) || m2_theo <= 0) return(1e10)
    
    ratio_theo <- m1_theo / m2_theo
    
    # Hem ratio'yu hem de mutlak değerleri kullan
    err_ratio <- (ratio_theo - ratio_emp)^2
    err1 <- ((m1_theo - m1_emp) / m1_emp)^2
    err2 <- ((m2_theo - m2_emp) / m2_emp)^2
    
    err_ratio + err1 + err2
  }
  
  starts <- list(
    c(theta_init, beta_init, gamma_init),
    c(theta_init * 0.8, beta_init * 1.2, gamma_init),
    c(theta_init * 1.2, beta_init * 0.8, gamma_init),
    c(theta_init, beta_init, gamma_init * 0.5),
    c(sd(y), 1.5, 0)
  )
  
  best <- NULL
  best_val <- Inf
  
  for (start in starts) {
    res <- tryCatch({
      optim(
        par = start,
        fn = objective,
        method = "L-BFGS-B",
        lower = c(1e-6, 0.3, -Inf),
        upper = c(Inf, 15, y0_1 - 1e-6),
        control = list(maxit = 2000)
      )
    }, error = function(e) NULL)
    
    if (!is.null(res) && res$value < best_val) {
      best_val <- res$value
      best <- res
    }
  }
  
  if (is.null(best)) {
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  list(
    theta = best$par[1],
    beta = best$par[2],
    gamma = best$par[3],
    objective = best_val
  )
}

# ==============================================================================
# MOM FONKSİYONU
# ==============================================================================

weibull_3p_mom <- function(y) {
  n <- length(y)
  y_bar <- mean(y)
  s2 <- var(y) * (n - 1) / n
  
  m2 <- mean((y - y_bar)^2)
  m3 <- mean((y - y_bar)^3)
  skew_data <- m3 / (m2^(3/2))
  
  weibull_skewness <- function(beta) {
    if (beta <= 0) return(NA)
    g1 <- gamma(1 + 1/beta)
    g2 <- gamma(1 + 2/beta)
    g3 <- gamma(1 + 3/beta)
    (g3 - 3*g1*g2 + 2*g1^3) / (g2 - g1^2)^(3/2)
  }
  
  opt_beta <- optimize(function(b) (weibull_skewness(b) - skew_data)^2, c(0.1, 20))
  beta_hat <- opt_beta$minimum
  
  g1 <- gamma(1 + 1/beta_hat)
  g2 <- gamma(1 + 2/beta_hat)
  
  theta_hat <- sqrt(s2 / (g2 - g1^2))
  gamma_hat <- y_bar - theta_hat * g1
  
  list(beta = beta_hat, theta = theta_hat, gamma = gamma_hat)
}

# ==============================================================================
# SİMÜLASYON
# ==============================================================================

cat("=" , rep("=", 70), "\n", sep = "")
cat("YENİ MRL YÖNTEMLERİ SİMÜLASYONU (n=100, 500 tekrar)\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

set.seed(123)

beta_true <- 1.5
theta_true <- 1.0
gamma_true <- 0.5

n_sim <- 500
n <- 100

mom_res <- matrix(NA, n_sim, 3)
mrl_multi_res <- matrix(NA, n_sim, 3)
mrl_twostage_res <- matrix(NA, n_sim, 3)
mrl_ratio_res <- matrix(NA, n_sim, 3)

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (i in 1:n_sim) {
  x <- rweibull(n, shape = beta_true, scale = theta_true)
  y <- gamma_true + x
  
  # MOM
  tryCatch({
    fit <- weibull_3p_mom(y)
    mom_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
  }, error = function(e) NULL)
  
  # MRL Multi-threshold
  tryCatch({
    fit <- weibull_3p_mrl_multi(y, threshold_probs = c(0.15, 0.35, 0.55))
    if (!is.na(fit$beta)) {
      mrl_multi_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
    }
  }, error = function(e) NULL)
  
  # MRL Two-stage
  tryCatch({
    fit <- weibull_3p_mrl_twostage(y, threshold_prob = 0.25)
    if (!is.na(fit$beta)) {
      mrl_twostage_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
    }
  }, error = function(e) NULL)
  
  # MRL Ratio
  tryCatch({
    fit <- weibull_3p_mrl_ratio(y, threshold_probs = c(0.20, 0.50))
    if (!is.na(fit$beta)) {
      mrl_ratio_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
    }
  }, error = function(e) NULL)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# İstatistikleri hesapla
calc_stats <- function(res, true_vals, name) {
  valid <- complete.cases(res)
  n_valid <- sum(valid)
  
  if (n_valid < 10) {
    return(data.frame(
      method = name,
      beta_mean = NA, beta_bias = NA, beta_mse = NA,
      theta_mean = NA, theta_bias = NA, theta_mse = NA,
      gamma_mean = NA, gamma_bias = NA, gamma_mse = NA,
      n_valid = n_valid
    ))
  }
  
  res_clean <- res[valid, , drop = FALSE]
  
  means <- colMeans(res_clean)
  vars <- apply(res_clean, 2, var)
  bias <- means - true_vals
  mse <- bias^2 + vars
  
  data.frame(
    method = name,
    beta_mean = means[1], beta_bias = bias[1], beta_mse = mse[1],
    theta_mean = means[2], theta_bias = bias[2], theta_mse = mse[2],
    gamma_mean = means[3], gamma_bias = bias[3], gamma_mse = mse[3],
    n_valid = n_valid
  )
}

true_vals <- c(beta_true, theta_true, gamma_true)

results <- rbind(
  calc_stats(mom_res, true_vals, "MOM"),
  calc_stats(mrl_multi_res, true_vals, "MRL_multi"),
  calc_stats(mrl_twostage_res, true_vals, "MRL_2stage"),
  calc_stats(mrl_ratio_res, true_vals, "MRL_ratio")
)

cat("\n\nSonuçlar:\n")
cat("-" , rep("-", 90), "\n", sep = "")
print(results, digits = 4)

cat("\n\nGerçek değerler: beta =", beta_true, ", theta =", theta_true, ", gamma =", gamma_true, "\n")

# MSE karşılaştırması
cat("\n\nMSE Karşılaştırması:\n")
cat("-" , rep("-", 50), "\n", sep = "")
mse_table <- results[, c("method", "beta_mse", "theta_mse", "gamma_mse")]
print(mse_table, digits = 4)

# En iyi yöntem
total_mse <- rowSums(results[, c("beta_mse", "theta_mse", "gamma_mse")], na.rm = TRUE)
cat("\n\nToplam MSE sıralaması:\n")
ord <- order(total_mse)
for (i in ord) {
  cat("  ", results$method[i], ": ", round(total_mse[i], 4), "\n", sep = "")
}

# ==============================================================================
# DAHA BÜYÜK ÖRNEKLEM
# ==============================================================================

cat("\n\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat("BÜYÜK ÖRNEKLEM TESTİ (n=500, 200 tekrar)\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

n_sim2 <- 200
n2 <- 500

mom_res2 <- matrix(NA, n_sim2, 3)
mrl_multi_res2 <- matrix(NA, n_sim2, 3)
mrl_twostage_res2 <- matrix(NA, n_sim2, 3)
mrl_ratio_res2 <- matrix(NA, n_sim2, 3)

pb <- txtProgressBar(min = 0, max = n_sim2, style = 3)

for (i in 1:n_sim2) {
  x <- rweibull(n2, shape = beta_true, scale = theta_true)
  y <- gamma_true + x
  
  tryCatch({
    fit <- weibull_3p_mom(y)
    mom_res2[i, ] <- c(fit$beta, fit$theta, fit$gamma)
  }, error = function(e) NULL)
  
  tryCatch({
    fit <- weibull_3p_mrl_multi(y, threshold_probs = c(0.15, 0.35, 0.55))
    if (!is.na(fit$beta)) mrl_multi_res2[i, ] <- c(fit$beta, fit$theta, fit$gamma)
  }, error = function(e) NULL)
  
  tryCatch({
    fit <- weibull_3p_mrl_twostage(y, threshold_prob = 0.25)
    if (!is.na(fit$beta)) mrl_twostage_res2[i, ] <- c(fit$beta, fit$theta, fit$gamma)
  }, error = function(e) NULL)
  
  tryCatch({
    fit <- weibull_3p_mrl_ratio(y, threshold_probs = c(0.20, 0.50))
    if (!is.na(fit$beta)) mrl_ratio_res2[i, ] <- c(fit$beta, fit$theta, fit$gamma)
  }, error = function(e) NULL)
  
  setTxtProgressBar(pb, i)
}
close(pb)

results2 <- rbind(
  calc_stats(mom_res2, true_vals, "MOM"),
  calc_stats(mrl_multi_res2, true_vals, "MRL_multi"),
  calc_stats(mrl_twostage_res2, true_vals, "MRL_2stage"),
  calc_stats(mrl_ratio_res2, true_vals, "MRL_ratio")
)

cat("\n\nSonuçlar (n=500):\n")
print(results2, digits = 4)

cat("\n\nMSE Karşılaştırması (n=500):\n")
mse_table2 <- results2[, c("method", "beta_mse", "theta_mse", "gamma_mse")]
print(mse_table2, digits = 4)