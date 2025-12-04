################################################################################
# 3 PARAMETRELİ WEIBULL DAĞILIMI İÇİN PARAMETRE TAHMİNİ
# Y ~ Weibull(beta, theta, gamma)  
# beta : şekil (shape), 
# theta: ölçek (scale), 
# gamma: konum (location)
#
# Y = X + gamma,  X ~ Weibull(beta, theta) (2 parametreli standart)
################################################################################

rm(list = ls())

# ==============================================================================
# BÖLÜM 1: TEORİK ALTYAPI
# ==============================================================================
# 
# 2 parametreli Weibull için:
#   E(X)   = theta * Gamma(1 + 1/beta)
#   E(X^2) = theta^2 * Gamma(1 + 2/beta)
#   E(X^3) = theta^3 * Gamma(1 + 3/beta)
#
# 3 parametreli Weibull için (Y = X + gamma):
#   E(Y) = gamma + theta * Gamma(1 + 1/beta)
#   E(Y^2) = E[(X+gamma)^2] 
#          = E(X^2) + 2*gamma*E(X) + gamma^2
#   E(Y^3) = E[(X+gamma)^3] 
#          = E(X^3) + 3*gamma*E(X^2) + 3*gamma^2*E(X) + gamma^3
#
# Önemli ilişkiler:
#   Var(Y) = Var(X) 
#          = theta^2 * [Gamma(1+2/beta) - Gamma(1+1/beta)^2]
#
#   Çarpıklık (Skewness): sadece beta'ya bağlı (theta ve gamma'dan bağımsız)
#
# STRATEJİ:
# 1. Çarpıklıktan beta'yı bul (çarpıklık sadece beta'ya bağlı)
# 2. Varyanstan theta'yı bul  
# 3. Ortalamadan gamma'yı bul
# ==============================================================================


# ==============================================================================
# BÖLÜM 2: MOM (Momentler Yöntemi) - SADE VE TUTARLI VERSİYON
# ==============================================================================

weibull_3p_mom <- function(y) {
  y <- as.numeric(y)
  n <- length(y)
  
  # Örnek istatistikler
  y_bar <- mean(y)                        # E(Y) tahmini
  m2    <- mean((y - y_bar)^2)           # 2. merkezi moment (Var(Y), n ile bölünmüş)
  m3    <- mean((y - y_bar)^3)           # 3. merkezi moment
  
  # Çarpıklık katsayısı (skewness)
  skew_data <- m3 / (m2^(3/2))
  
  # -----------------------------------------------------
  # ADIM 1: Beta'yı bul (çarpıklıktan)
  # -----------------------------------------------------
  # Weibull çarpıklık formülü (sadece beta'ya bağlı):
  # skew = [Γ(1+3/β) - 3*Γ(1+1/β)*Γ(1+2/β) + 2*Γ(1+1/β)^3] / 
  #        [Γ(1+2/β) - Γ(1+1/β)^2]^(3/2)
  
  weibull_skewness <- function(beta) {
    if (beta <= 0) return(NA_real_)
    g1 <- gamma(1 + 1/beta)
    g2 <- gamma(1 + 2/beta)
    g3 <- gamma(1 + 3/beta)
    
    numerator   <- g3 - 3*g1*g2 + 2*g1^3
    denominator <- (g2 - g1^2)^(3/2)
    
    numerator / denominator
  }
  
  objective_beta <- function(beta) {
    s_th <- weibull_skewness(beta)
    if (!is.finite(s_th)) return(1e9)
    (s_th - skew_data)^2
  }
  
  # Beta için optimizasyon
  opt_beta <- optimize(objective_beta, interval = c(0.1, 20))
  beta_hat <- opt_beta$minimum
  
  # -----------------------------------------------------
  # ADIM 2: Theta'yı bul (varyanstan)
  # -----------------------------------------------------
  # Var(Y) = Var(X) = theta^2 * [Γ(1+2/β) - Γ(1+1/β)^2]
  # => theta = sqrt(Var(Y) / [Γ(1+2/β) - Γ(1+1/β)^2])
  
  g1 <- gamma(1 + 1/beta_hat)
  g2 <- gamma(1 + 2/beta_hat)
  
  var_factor <- g2 - g1^2
  theta_hat  <- sqrt(m2 / var_factor)
  
  # -----------------------------------------------------
  # ADIM 3: Gamma'yı bul (ortalamadan)
  # -----------------------------------------------------
  # E(Y) = gamma + theta * Γ(1+1/β)
  # => gamma = E(Y) - theta * Γ(1+1/β)
  
  gamma_hat <- y_bar - theta_hat * g1
  
  # Sonuçları döndür
  list(
    beta            = beta_hat,
    theta           = theta_hat,
    gamma           = gamma_hat,
    skewness_data   = skew_data,
    skewness_fitted = weibull_skewness(beta_hat)
  )
}


# ==============================================================================
# BÖLÜM 3: MLE (Maximum Likelihood Estimation)
# ==============================================================================

weibull_3p_mle <- function(y, start = NULL) {
  y     <- as.numeric(y)
  n     <- length(y)
  y_min <- min(y)
  
  # Negatif log-likelihood
  neg_loglik <- function(par) {
    beta  <- par[1]   # shape
    theta <- par[2]   # scale
    gamma <- par[3]   # location
    
    # Kısıtlamalar
    if (beta <= 0 || theta <= 0 || gamma >= y_min) return(Inf)
    
    # Shifted değerler
    x <- y - gamma
    if (any(x <= 0)) return(Inf)
    
    # Log-likelihood (Weibull(beta, theta), x > 0)
    ll <- n * log(beta) - n * beta * log(theta) +
      (beta - 1) * sum(log(x)) -
      sum((x/theta)^beta)
    
    -ll
  }
  
  # Başlangıç değerleri (MOM'dan al)
  if (is.null(start)) {
    mom_est <- tryCatch(weibull_3p_mom(y), error = function(e) NULL)
    if (!is.null(mom_est)) {
      start <- c(mom_est$beta, mom_est$theta, mom_est$gamma)
    } else {
      start <- c(1.5, sd(y), y_min * 0.5)
    }
  }
  
  fit <- optim(
    par    = start,
    fn     = neg_loglik,
    method = "L-BFGS-B",
    lower  = c(1e-6, 1e-6, -Inf),
    upper  = c(Inf,  Inf,  y_min - 1e-6)
  )
  
  list(
    beta        = fit$par[1],
    theta       = fit$par[2],
    gamma       = fit$par[3],
    neg_loglik  = fit$value,
    convergence = fit$convergence
  )
}


# ==============================================================================
# BÖLÜM 4: MRL (Mean Residual Life) YÖNTEMİ - DÜZELTİLMİŞ VERSİYON
# ==============================================================================
#
# MRL: E[Y - y0 | Y > y0] = mean residual life at y0
#
# 3 parametreli Weibull için (Y = X + gamma):
# E[Y - y0 | Y > y0] = theta * exp(t) * Γ_upper(1 + 1/β, t) - (y0 - γ)
#  burada t = ((y0 - γ) / θ)^β
#
# Benzer şekilde 2. ve 3. momentler için:
#
# E[(Y - y0)^2 | Y > y0]  = theta^2 * e^t * G2 
#                           - 2 (y0 - γ) theta e^t G1 + (y0 - γ)^2
#
# E[(Y - y0)^3 | Y > y0]  = theta^3 * e^t * G3
#                           - 3 (y0 - γ) theta^2 e^t G2
#                           + 3 (y0 - γ)^2 theta e^t G1
#                           - (y0 - γ)^3
#
# Burada Gk = Γ_upper(1 + k/β, t)
# ==============================================================================

# Üst eksik gamma fonksiyonu: Γ(s, z) = integral_z^inf t^(s-1)*e^(-t) dt
gamma_upper <- function(s, z) {
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

weibull_3p_mrl <- function(y, threshold_prob = 0.25) {
  y     <- as.numeric(y)
  n     <- length(y)
  y_min <- min(y)
  
  # Threshold değeri
  y0 <- quantile(y, probs = threshold_prob, names = FALSE)
  
  # Threshold üzerindeki gözlemler
  exc   <- y[y > y0] - y0
  n_exc <- length(exc)
  
  if (n_exc < 10) {
    warning("Threshold üzerinde yeterli gözlem yok")
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  # Empirik MRL momentleri
  mrl_1 <- mean(exc)        # E[Y - y0 | Y > y0]
  mrl_2 <- mean(exc^2)      # E[(Y - y0)^2 | Y > y0]
  mrl_3 <- mean(exc^3)      # E[(Y - y0)^3 | Y > y0]
  
  # -----------------------------------------------------
  # MRL denklem sistemi - kare toplamını minimize et
  # Parametre sırası: par = c(beta, theta, gamma)
  # -----------------------------------------------------
  mrl_objective <- function(par) {
    beta  <- par[1]  # shape
    theta <- par[2]  # scale
    gamma <- par[3]  # location
    
    # Kısıtlamalar
    if (beta <= 0 || theta <= 0 || gamma >= y0) {
      return(1e10)
    }
    
    # z = y0 - gamma > 0 olmalı
    z <- y0 - gamma
    if (z <= 0) return(1e10)
    
    # t = ((y0 - gamma) / theta)^beta
    t <- (z / theta)^beta
    if (!is.finite(t) || t <= 0 || t > 700) return(1e10)
    
    exp_t <- exp(t)
    
    # Üst eksik gamma değerleri
    G1 <- gamma_upper(1 + 1/beta, t)
    G2 <- gamma_upper(1 + 2/beta, t)
    G3 <- gamma_upper(1 + 3/beta, t)
    
    # Teorik MRL momentleri
    theo_1 <- theta * exp_t * G1 - z
    theo_2 <- theta^2 * exp_t * G2 -
      2 * z * theta * exp_t * G1 + z^2
    theo_3 <- theta^3 * exp_t * G3 -
      3 * z * theta^2 * exp_t * G2 +
      3 * z^2 * theta * exp_t * G1 -
      z^3
    
    # Normalize edilmiş kare farklar toplamı
    err1 <- (theo_1 - mrl_1)^2 / (mrl_1^2 + 1e-10)
    err2 <- (theo_2 - mrl_2)^2 / (mrl_2^2 + 1e-10)
    err3 <- (theo_3 - mrl_3)^2 / (mrl_3^2 + 1e-10)
    
    err1 + err2 + err3
  }
  
  # Başlangıç değerleri için MOM kullan
  mom_est <- tryCatch(weibull_3p_mom(y), error = function(e) NULL)
  
  if (!is.null(mom_est) && !is.na(mom_est$theta)) {
    # MOM tahminleri: beta, theta, gamma
    start_values <- list(
      c(mom_est$beta,       mom_est$theta,       mom_est$gamma),
      c(mom_est$beta * 1.2, mom_est$theta * 0.8, mom_est$gamma * 0.9),
      c(mom_est$beta * 0.8, mom_est$theta * 1.2, mom_est$gamma * 1.1),
      c(1.5,                sd(y),              y_min * 0.5),
      c(2.0,                sd(y) * 1.5,        y_min * 0.3)
    )
  } else {
    start_values <- list(
      c(1.5, sd(y),      y_min * 0.5),
      c(2.0, sd(y) * 1.5,y_min * 0.3),
      c(1.0, sd(y) * 0.7,y_min * 0.7),
      c(2.0, mean(y),    0),
      c(3.0, sd(y),      y_min * 0.8)
    )
  }
  
  best_result <- NULL
  best_value  <- Inf
  
  for (start in start_values) {
    result <- tryCatch({
      optim(
        par    = start,
        fn     = mrl_objective,
        method = "L-BFGS-B",
        lower  = c(1e-6, 1e-6, -Inf),
        upper  = c(Inf,  Inf,  y0 - 1e-6),
        control = list(maxit = 1000)
      )
    }, error = function(e) NULL)
    
    if (!is.null(result) && is.finite(result$value) && result$value < best_value) {
      best_value  <- result$value
      best_result <- result
    }
  }
  
  if (is.null(best_result) || !is.finite(best_value) || best_value > 1) {
    warning(paste("MRL yakınsamadı, en iyi değer:", round(best_value, 4)))
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  list(
    beta            = best_result$par[1],
    theta           = best_result$par[2],
    gamma           = best_result$par[3],
    objective_value = best_value,
    convergence     = best_result$convergence
  )
}


# ==============================================================================
# BÖLÜM 5: SİMÜLASYON ÇALIŞMASI
# ==============================================================================

run_simulation <- function(n = 50, 
                           beta_true = 1.5, 
                           theta_true = 1.0, 
                           gamma_true = 0.5,
                           n_sim = 100,
                           add_outliers = FALSE) {
  
  cat("Simülasyon Parametreleri:\n")
  cat("  n =", n, ", beta =", beta_true, ", theta =", theta_true, ", gamma =", gamma_true, "\n")
  cat("  Tekrar sayısı:", n_sim, "\n")
  cat("  Aykırı değer:", ifelse(add_outliers, "Evet", "Hayır"), "\n\n")
  
  # Sonuç matrisleri
  results_mom <- matrix(NA, nrow = n_sim, ncol = 3)
  results_mle <- matrix(NA, nrow = n_sim, ncol = 3)
  results_mrl <- matrix(NA, nrow = n_sim, ncol = 3)
  colnames(results_mom) <- colnames(results_mle) <- colnames(results_mrl) <- c("beta", "theta", "gamma")
  
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  
  for (i in 1:n_sim) {
    # Veri üret: Y = gamma + X, X ~ Weibull(beta, theta)
    x <- rweibull(n, shape = beta_true, scale = theta_true)
    y <- gamma_true + x
    
    # Aykırı değer ekle
    if (add_outliers) {
      n_out <- floor(n / 10)
      idx   <- sample(1:n, n_out)
      y[idx] <- y[idx] * 2
    }
    
    # MOM
    tryCatch({
      mom_fit <- weibull_3p_mom(y)
      results_mom[i, ] <- c(mom_fit$beta, mom_fit$theta, mom_fit$gamma)
    }, error = function(e) NULL)
    
    # MLE
    tryCatch({
      mle_fit <- weibull_3p_mle(y)
      if (mle_fit$convergence == 0) {
        results_mle[i, ] <- c(mle_fit$beta, mle_fit$theta, mle_fit$gamma)
      }
    }, error = function(e) NULL)
    
    # MRL
    tryCatch({
      mrl_fit <- weibull_3p_mrl(y, threshold_prob = 0.25)
      if (!is.na(mrl_fit$beta)) {
        results_mrl[i, ] <- c(mrl_fit$beta, mrl_fit$theta, mrl_fit$gamma)
      }
    }, error = function(e) NULL)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Özet istatistikler
  calc_stats <- function(est, true_val, param_name) {
    est_clean <- est[!is.na(est)]
    n_valid   <- length(est_clean)
    
    if (n_valid < 5) {
      return(c(mean = NA, var = NA, bias = NA, mse = NA, n_valid = n_valid))
    }
    
    mean_est <- mean(est_clean)
    var_est  <- var(est_clean)
    bias     <- mean_est - true_val
    mse      <- bias^2 + var_est
    
    c(mean = mean_est, var = var_est, bias = bias, mse = mse, n_valid = n_valid)
  }
  
  true_vals <- c(beta = beta_true, theta = theta_true, gamma = gamma_true)
  
  summary_df <- data.frame(
    Method    = rep(c("MOM", "MLE", "MRL"), each = 3),
    Parameter = rep(c("beta", "theta", "gamma"), 3),
    True      = rep(true_vals, 3)
  )
  
  # MOM stats
  for (j in 1:3) {
    stats <- calc_stats(results_mom[, j], true_vals[j], names(true_vals)[j])
    summary_df[j, c("Mean", "Var", "Bias", "MSE", "N_valid")] <- stats
  }
  
  # MLE stats
  for (j in 1:3) {
    stats <- calc_stats(results_mle[, j], true_vals[j], names(true_vals)[j])
    summary_df[3 + j, c("Mean", "Var", "Bias", "MSE", "N_valid")] <- stats
  }
  
  # MRL stats
  for (j in 1:3) {
    stats <- calc_stats(results_mrl[, j], true_vals[j], names(true_vals)[j])
    summary_df[6 + j, c("Mean", "Var", "Bias", "MSE", "N_valid")] <- stats
  }
  
  summary_df
}


# ==============================================================================
# BÖLÜM 6: TEK ÖRNEK TESTİ
# ==============================================================================

cat("=====================================\n")
cat("TEK ÖRNEK TESTİ\n")
cat("=====================================\n\n")

set.seed(42)
# Gerçek parametreler
beta_true  <- 1.5
theta_true <- 1.0
gamma_true <- 0.5

# Veri üret
n <- 100
x <- rweibull(n, shape = beta_true, scale = theta_true)
y <- gamma_true + x

cat("Gerçek parametreler:\n")
cat("  beta  =", beta_true,  "\n")
cat("  theta =", theta_true, "\n")
cat("  gamma =", gamma_true, "\n\n")

cat("Örnek istatistikler:\n")
cat("  n         =", length(y),          "\n")
cat("  Ortalama  =", round(mean(y), 4),  "\n")
cat("  Std Sapma =", round(sd(y), 4),    "\n")
cat("  Min       =", round(min(y), 4),   "\n")
cat("  Max       =", round(max(y), 4),   "\n\n")

# MOM
cat("MOM Tahminleri:\n")
mom_fit <- weibull_3p_mom(y)
cat("  beta  =", round(mom_fit$beta,  4), "\n")
cat("  theta =", round(mom_fit$theta, 4), "\n")
cat("  gamma =", round(mom_fit$gamma, 4), "\n\n")

# MLE
cat("MLE Tahminleri:\n")
mle_fit <- weibull_3p_mle(y)
cat("  beta  =", round(mle_fit$beta,  4), "\n")
cat("  theta =", round(mle_fit$theta, 4), "\n")
cat("  gamma =", round(mle_fit$gamma, 4), "\n\n")

# MRL
cat("MRL Tahminleri:\n")
mrl_fit <- weibull_3p_mrl(y)
if (!is.na(mrl_fit$beta)) {
  cat("  beta  =", round(mrl_fit$beta,  4), "\n")
  cat("  theta =", round(mrl_fit$theta, 4), "\n")
  cat("  gamma =", round(mrl_fit$gamma, 4), "\n")
  cat("  Objective value =", round(mrl_fit$objective_value, 6), "\n\n")
} else {
  cat("  MRL yakınsamadı\n\n")
}


# ==============================================================================
# BÖLÜM 7: SİMÜLASYON SONUÇLARI
# ==============================================================================

cat("=====================================\n")
cat("SİMÜLASYON SONUÇLARI (AYKIRISIZ)\n")
cat("=====================================\n\n")

sim_results <- run_simulation(
  n          = 100, 
  beta_true  = 1.5, 
  theta_true = 1.0, 
  gamma_true = 0.5,
  n_sim      = 200,
  add_outliers = FALSE
)

print(sim_results)

cat("\n=====================================\n")
cat("AYKIRI DEĞERLİ SİMÜLASYON\n")
cat("=====================================\n\n")

sim_results_outlier <- run_simulation(
  n          = 100, 
  beta_true  = 1.5, 
  theta_true = 1.0, 
  gamma_true = 0.5,
  n_sim      = 200,
  add_outliers = TRUE
)

print(sim_results_outlier)


# ==============================================================================
# BÖLÜM 8: MAKALE VERİSİ (NA_a - Naturally Aged Glass)
# ==============================================================================

cat("\n=====================================\n")
cat("MAKALE VERİSİ (NA_a - Naturally Aged Glass)\n")
cat("=====================================\n\n")

dataset1 <- c(
  24.12, 24.13, 28.52, 29.18, 29.67, 30.48, 32.98, 35.91, 35.92,
  36.38, 37.60, 37.70, 39.71, 49.10, 52.43, 52.46, 52.61, 61.72
)

cat("Veri özeti:\n")
cat("  n         =", length(dataset1),          "\n")
cat("  Ortalama  =", round(mean(dataset1), 2),  "\n")
cat("  Std Sapma =", round(sd(dataset1),   2),  "\n")
cat("  Min       =", min(dataset1),            "\n")
cat("  Max       =", max(dataset1),            "\n\n")

cat("MOM Tahminleri:\n")
mom_real <- weibull_3p_mom(dataset1)
cat("  beta  =", round(mom_real$beta,  4), "\n")
cat("  theta =", round(mom_real$theta, 4), "\n")
cat("  gamma =", round(mom_real$gamma, 4), "\n\n")

cat("MLE Tahminleri:\n")
mle_real <- weibull_3p_mle(dataset1)
cat("  beta  =", round(mle_real$beta,  4), "\n")
cat("  theta =", round(mle_real$theta, 4), "\n")
cat("  gamma =", round(mle_real$gamma, 4), "\n\n")

cat("MRL Tahminleri:\n")
mrl_real <- weibull_3p_mrl(dataset1, threshold_prob = 0.25)
if (!is.na(mrl_real$beta)) {
  cat("  beta  =", round(mrl_real$beta,  4), "\n")
  cat("  theta =", round(mrl_real$theta, 4), "\n")
  cat("  gamma =", round(mrl_real$gamma, 4), "\n")
} else {
  cat("  MRL yakınsamadı\n")
}

cat("\n")
cat("Makaledeki PITE sonuçları (karşılaştırma için):\n")
cat("  rho = 0.32: beta = 2.399, theta = 37.961, mu = 6.976\n")
cat("  rho = 0.42: beta = 2.369, theta = 38.098, mu = 6.828\n")
cat("  rho = 0.68: beta = 2.197, theta = 38.715, mu = 5.972\n")
cat("  rho = 1.00: beta = 1.967, theta = 39.524, mu = 4.746\n")
