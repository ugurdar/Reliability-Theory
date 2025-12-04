################################################################################
# MRL YÖNTEMİNİN DETAYLI ANALİZİ
# Problem: MRL neden kötü sonuç veriyor?
################################################################################

rm(list = ls())

# ==============================================================================
# BÖLÜM 1: TEORİK KONTROL
# ==============================================================================
# 
# 3 parametreli Weibull: Y = gamma + X, X ~ Weibull(beta, theta)
# 
# CDF: F(y) = 1 - exp(-((y-gamma)/theta)^beta),  y > gamma
# PDF: f(y) = (beta/theta) * ((y-gamma)/theta)^(beta-1) * exp(-((y-gamma)/theta)^beta)
# 
# Mean Residual Life (MRL) at threshold y0:
# m(y0) = E[Y - y0 | Y > y0]
#
# Genel formül:
# m(y0) = (1/S(y0)) * integral_{y0}^{inf} S(y) dy
# burada S(y) = 1 - F(y) = exp(-((y-gamma)/theta)^beta)
#
# Değişken dönüşümü: u = ((y-gamma)/theta)^beta
# y = gamma + theta * u^(1/beta)
# dy = (theta/beta) * u^(1/beta - 1) du
#
# y = y0 için: u0 = ((y0-gamma)/theta)^beta = t diyelim
#
# MRL formülü:
# m(y0) = theta * exp(t) * Gamma_upper(1 + 1/beta, t) - (y0 - gamma)
#
# Bu formül DOĞRU görünüyor. Peki 2. ve 3. momentler?
# ==============================================================================

# Üst eksik gamma fonksiyonu
gamma_upper <- function(s, z) {
  if (z < 0) return(gamma(s))
  if (z > 700) return(0)  # Overflow önleme
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

# ==============================================================================
# BÖLÜM 2: TEORİK MRL MOMENTLER - TÜRETME
# ==============================================================================
#
# E[(Y - y0)^k | Y > y0] için genel formül:
#
# Substitution: Z = Y - gamma (standart 2-param Weibull)
# z0 = y0 - gamma
#
# E[(Y - y0)^k | Y > y0] = E[(Z - z0)^k | Z > z0]
#
# 1. MOMENT (k=1):
# E[Z - z0 | Z > z0] = theta * exp(t) * Gamma(1+1/beta, t) - z0
# burada t = (z0/theta)^beta
#
# 2. MOMENT (k=2):
# E[(Z - z0)^2 | Z > z0] = E[Z^2 | Z > z0] - 2*z0*E[Z | Z > z0] + z0^2
#
# E[Z^2 | Z > z0] = (1/S(z0)) * integral_{z0}^{inf} z^2 * f(z) dz
#                 = theta^2 * exp(t) * Gamma(1 + 2/beta, t)
#
# E[Z | Z > z0] = theta * exp(t) * Gamma(1 + 1/beta, t)
#
# Dolayısıyla:
# E[(Z-z0)^2 | Z>z0] = theta^2*exp(t)*G(1+2/b,t) - 2*z0*theta*exp(t)*G(1+1/b,t) + z0^2
#
# 3. MOMENT (k=3):
# E[(Z-z0)^3 | Z>z0] = E[Z^3|Z>z0] - 3*z0*E[Z^2|Z>z0] + 3*z0^2*E[Z|Z>z0] - z0^3
#
# E[Z^3 | Z > z0] = theta^3 * exp(t) * Gamma(1 + 3/beta, t)
# ==============================================================================

# Teorik MRL momentleri hesapla (DOĞRULANMIŞ FORMÜLLER)
calc_theoretical_mrl <- function(y0, beta, theta, gamma) {
  z0 <- y0 - gamma
  if (z0 <= 0) return(list(m1 = NA, m2 = NA, m3 = NA))
  
  t <- (z0 / theta)^beta
  exp_t <- exp(t)
  
  G1 <- gamma_upper(1 + 1/beta, t)
  G2 <- gamma_upper(1 + 2/beta, t)
  G3 <- gamma_upper(1 + 3/beta, t)
  
  # Conditional expectations E[Z^k | Z > z0]
  EZ_cond <- theta * exp_t * G1
  EZ2_cond <- theta^2 * exp_t * G2
  EZ3_cond <- theta^3 * exp_t * G3
  
  # MRL momentleri E[(Z - z0)^k | Z > z0]
  m1 <- EZ_cond - z0
  m2 <- EZ2_cond - 2*z0*EZ_cond + z0^2
  m3 <- EZ3_cond - 3*z0*EZ2_cond + 3*z0^2*EZ_cond - z0^3
  
  list(m1 = m1, m2 = m2, m3 = m3, t = t, z0 = z0)
}

# ==============================================================================
# BÖLÜM 3: EMPİRİK MRL - n vs n-1 KONTROLÜ
# ==============================================================================

# Empirik MRL momentleri hesapla - FARKLI VERSİYONLAR
calc_empirical_mrl <- function(y, y0, version = "biased") {
  # Threshold üzerindeki değerler
  exc <- y[y > y0] - y0
  n_exc <- length(exc)
  
  if (n_exc < 3) {
    return(list(m1 = NA, m2 = NA, m3 = NA, n_exc = n_exc))
  }
  
  if (version == "biased") {
    # Basit ortalamalar (MLE-tipi, biased)
    m1 <- mean(exc)
    m2 <- mean(exc^2)
    m3 <- mean(exc^3)
  } else if (version == "unbiased_var") {
    # 2. moment için unbiased varyans düzeltmesi
    m1 <- mean(exc)
    # E[(Y-y0)^2] = Var(Y-y0) + [E(Y-y0)]^2
    # Unbiased varyans: sum((x-mean)^2)/(n-1)
    var_exc <- var(exc)  # Bu zaten n-1 ile böler
    m2 <- var_exc + m1^2  # Var + mean^2 = E[X^2]
    # Ama bu yanlış! Biz E[(Y-y0)^2] istiyoruz, E[(Y-y0-mean)^2] + mean^2 değil
    # Aslında: E[(Y-y0)^2] = mean(exc^2), bu doğru
    m2 <- mean(exc^2)
    m3 <- mean(exc^3)
  } else if (version == "kernel") {
    # Kernel density estimation ile
    m1 <- mean(exc)
    m2 <- mean(exc^2)
    m3 <- mean(exc^3)
  }
  
  list(m1 = m1, m2 = m2, m3 = m3, n_exc = n_exc)
}

# ==============================================================================
# BÖLÜM 4: TEORİK DEĞERLERLE KARŞILAŞTIRMA
# ==============================================================================

cat("=" , rep("=", 60), "\n", sep = "")
cat("TEORİK VS EMPİRİK MRL KARŞILAŞTIRMASI\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

set.seed(123)

# Gerçek parametreler
beta_true <- 1.5
theta_true <- 1.0
gamma_true <- 0.5

# Büyük örneklem (teorik değerlere yakınsamalı)
n <- 10000
x <- rweibull(n, shape = beta_true, scale = theta_true)
y <- gamma_true + x

# Farklı threshold değerleri dene
probs <- c(0.10, 0.25, 0.50, 0.75)

cat("Gerçek parametreler: beta =", beta_true, ", theta =", theta_true, ", gamma =", gamma_true, "\n\n")

for (p in probs) {
  y0 <- quantile(y, probs = p, names = FALSE)
  
  # Teorik değerler
  theo <- calc_theoretical_mrl(y0, beta_true, theta_true, gamma_true)
  
  # Empirik değerler
  emp <- calc_empirical_mrl(y, y0)
  
  cat("Threshold quantile:", p, " (y0 =", round(y0, 4), ")\n")
  cat("  n_exc =", emp$n_exc, "\n")
  cat("  Teorik m1 =", round(theo$m1, 6), ", Empirik m1 =", round(emp$m1, 6), 
      ", Fark =", round(emp$m1 - theo$m1, 6), "\n")
  cat("  Teorik m2 =", round(theo$m2, 6), ", Empirik m2 =", round(emp$m2, 6),
      ", Fark =", round(emp$m2 - theo$m2, 6), "\n")
  cat("  Teorik m3 =", round(theo$m3, 6), ", Empirik m3 =", round(emp$m3, 6),
      ", Fark =", round(emp$m3 - theo$m3, 6), "\n\n")
}

# ==============================================================================
# BÖLÜM 5: MRL DENKLEM SİSTEMİNİN CONDITION NUMBER ANALİZİ
# ==============================================================================

cat("=" , rep("=", 60), "\n", sep = "")
cat("MRL DENKLEM SİSTEMİ JACOBIAN ANALİZİ\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

# MRL denklemleri
mrl_equations <- function(par, y0, m1_emp, m2_emp, m3_emp) {
  theta <- par[1]
  beta <- par[2]
  gamma <- par[3]
  
  if (theta <= 0 || beta <= 0 || gamma >= y0) {
    return(c(1e10, 1e10, 1e10))
  }
  
  theo <- calc_theoretical_mrl(y0, beta, theta, gamma)
  
  c(theo$m1 - m1_emp, theo$m2 - m2_emp, theo$m3 - m3_emp)
}

# Sayısal Jacobian hesapla
numerical_jacobian <- function(par, y0, m1_emp, m2_emp, m3_emp, eps = 1e-6) {
  n_par <- length(par)
  f0 <- mrl_equations(par, y0, m1_emp, m2_emp, m3_emp)
  n_eq <- length(f0)
  
  J <- matrix(0, n_eq, n_par)
  
  for (j in 1:n_par) {
    par_plus <- par
    par_plus[j] <- par[j] + eps
    f_plus <- mrl_equations(par_plus, y0, m1_emp, m2_emp, m3_emp)
    J[, j] <- (f_plus - f0) / eps
  }
  
  J
}

# Gerçek parametrelerde Jacobian'ı kontrol et
y0_test <- quantile(y, 0.25, names = FALSE)
emp_test <- calc_empirical_mrl(y, y0_test)

J <- numerical_jacobian(
  c(theta_true, beta_true, gamma_true),
  y0_test, emp_test$m1, emp_test$m2, emp_test$m3
)

cat("Jacobian at true parameters:\n")
print(round(J, 4))

# Condition number
svd_J <- svd(J)
cond_num <- max(svd_J$d) / min(svd_J$d)
cat("\nCondition number:", round(cond_num, 2), "\n")
cat("Singular values:", round(svd_J$d, 6), "\n\n")

# ==============================================================================
# BÖLÜM 6: ÖLÇEKLEME PROBLEMİ
# ==============================================================================

cat("=" , rep("=", 60), "\n", sep = "")
cat("ÖLÇEKLEME PROBLEMİ ANALİZİ\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

# m1, m2, m3 çok farklı ölçeklerde!
cat("Moment büyüklükleri (n=10000 örnek):\n")
cat("  m1 =", round(emp_test$m1, 4), "\n")
cat("  m2 =", round(emp_test$m2, 4), "\n")
cat("  m3 =", round(emp_test$m3, 4), "\n")
cat("  m3/m1 oranı =", round(emp_test$m3/emp_test$m1, 2), "\n\n")

cat("Bu demek ki:\n")
cat("  - 3. moment denklemi domine ediyor\n")
cat("  - Optimizasyon m3'ü minimize etmeye odaklanıyor\n")
cat("  - m1 ve m2 ihmal ediliyor\n\n")

# ==============================================================================
# BÖLÜM 7: DÜZELTİLMİŞ MRL - ÖLÇEKLENMİŞ DENKLEMLER
# ==============================================================================

cat("=" , rep("=", 60), "\n", sep = "")
cat("DÜZELTİLMİŞ MRL YÖNTEMİ\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

# Ölçeklenmiş objective function
mrl_objective_scaled <- function(par, y0, m1_emp, m2_emp, m3_emp) {
  theta <- par[1]
  beta <- par[2]
  gamma <- par[3]
  
  if (theta <= 0 || beta <= 0 || gamma >= y0) {
    return(1e10)
  }
  
  z0 <- y0 - gamma
  if (z0 <= 0) return(1e10)
  
  t <- (z0 / theta)^beta
  if (t > 700) return(1e10)
  
  exp_t <- exp(t)
  G1 <- gamma_upper(1 + 1/beta, t)
  G2 <- gamma_upper(1 + 2/beta, t)
  G3 <- gamma_upper(1 + 3/beta, t)
  
  EZ_cond <- theta * exp_t * G1
  EZ2_cond <- theta^2 * exp_t * G2
  EZ3_cond <- theta^3 * exp_t * G3
  
  m1_theo <- EZ_cond - z0
  m2_theo <- EZ2_cond - 2*z0*EZ_cond + z0^2
  m3_theo <- EZ3_cond - 3*z0*EZ2_cond + 3*z0^2*EZ_cond - z0^3
  
  # ÖLÇEKLENMİŞ HATALAR (her momenti kendi büyüklüğüne böl)
  err1 <- ((m1_theo - m1_emp) / m1_emp)^2
  err2 <- ((m2_theo - m2_emp) / m2_emp)^2
  err3 <- ((m3_theo - m3_emp) / m3_emp)^2
  
  # Eşit ağırlık
  return(err1 + err2 + err3)
}

# Alternatif: Log-moment kullan
mrl_objective_log <- function(par, y0, m1_emp, m2_emp, m3_emp) {
  theta <- par[1]
  beta <- par[2]
  gamma <- par[3]
  
  if (theta <= 0 || beta <= 0 || gamma >= y0) {
    return(1e10)
  }
  
  z0 <- y0 - gamma
  if (z0 <= 0) return(1e10)
  
  t <- (z0 / theta)^beta
  if (t > 700) return(1e10)
  
  exp_t <- exp(t)
  G1 <- gamma_upper(1 + 1/beta, t)
  G2 <- gamma_upper(1 + 2/beta, t)
  G3 <- gamma_upper(1 + 3/beta, t)
  
  EZ_cond <- theta * exp_t * G1
  EZ2_cond <- theta^2 * exp_t * G2
  EZ3_cond <- theta^3 * exp_t * G3
  
  m1_theo <- EZ_cond - z0
  m2_theo <- EZ2_cond - 2*z0*EZ_cond + z0^2
  m3_theo <- EZ3_cond - 3*z0*EZ2_cond + 3*z0^2*EZ_cond - z0^3
  
  if (m1_theo <= 0 || m2_theo <= 0 || m3_theo <= 0) return(1e10)
  
  # LOG DÖNÜŞÜMÜ ile ölçekleme
  err1 <- (log(m1_theo) - log(m1_emp))^2
  err2 <- (log(m2_theo) - log(m2_emp))^2
  err3 <- (log(m3_theo) - log(m3_emp))^2
  
  return(err1 + err2 + err3)
}

# MOM ile başlangıç değeri
weibull_3p_mom <- function(y) {
  y_bar <- mean(y)
  n <- length(y)
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

# DÜZELTİLMİŞ MRL FONKSİYONU
weibull_3p_mrl_v2 <- function(y, threshold_prob = 0.25, method = "scaled") {
  n <- length(y)
  y_min <- min(y)
  y0 <- quantile(y, probs = threshold_prob, names = FALSE)
  
  exc <- y[y > y0] - y0
  n_exc <- length(exc)
  
  if (n_exc < 10) {
    warning("Yetersiz gözlem")
    return(list(beta = NA, theta = NA, gamma = NA))
  }
  
  m1_emp <- mean(exc)
  m2_emp <- mean(exc^2)
  m3_emp <- mean(exc^3)
  
  # Objective seç
  if (method == "scaled") {
    obj_fn <- function(par) mrl_objective_scaled(par, y0, m1_emp, m2_emp, m3_emp)
  } else {
    obj_fn <- function(par) mrl_objective_log(par, y0, m1_emp, m2_emp, m3_emp)
  }
  
  # MOM'dan başlangıç
  mom <- tryCatch(weibull_3p_mom(y), error = function(e) NULL)
  
  if (!is.null(mom) && !is.na(mom$theta) && mom$theta > 0 && mom$beta > 0) {
    starts <- list(
      c(mom$theta, mom$beta, mom$gamma),
      c(mom$theta * 0.8, mom$beta * 1.2, mom$gamma),
      c(mom$theta * 1.2, mom$beta * 0.8, mom$gamma),
      c(mom$theta, mom$beta, mom$gamma * 0.5),
      c(mom$theta, mom$beta, 0)
    )
  } else {
    starts <- list(
      c(sd(y), 1.5, y_min * 0.5),
      c(sd(y), 2.0, 0),
      c(sd(y) * 0.5, 1.0, y_min * 0.8)
    )
  }
  
  best <- NULL
  best_val <- Inf
  
  for (start in starts) {
    if (any(is.na(start))) next
    
    res <- tryCatch({
      optim(
        par = start,
        fn = obj_fn,
        method = "L-BFGS-B",
        lower = c(1e-6, 0.1, -Inf),
        upper = c(Inf, 20, y0 - 1e-6),
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
    convergence = best$convergence
  )
}

# ==============================================================================
# BÖLÜM 8: KÜÇÜK ÖRNEKLEM TESTİ
# ==============================================================================

cat("\nKüçük örneklem testi (n=50):\n")
cat("-" , rep("-", 40), "\n", sep = "")

set.seed(456)
n_test <- 50
x_test <- rweibull(n_test, shape = beta_true, scale = theta_true)
y_test <- gamma_true + x_test

mom_test <- weibull_3p_mom(y_test)
mrl_scaled <- weibull_3p_mrl_v2(y_test, method = "scaled")
mrl_log <- weibull_3p_mrl_v2(y_test, method = "log")

cat("\nGerçek: beta =", beta_true, ", theta =", theta_true, ", gamma =", gamma_true, "\n\n")

cat("MOM:\n")
cat("  beta =", round(mom_test$beta, 4), "\n")
cat("  theta =", round(mom_test$theta, 4), "\n")
cat("  gamma =", round(mom_test$gamma, 4), "\n\n")

cat("MRL (scaled):\n")
cat("  beta =", round(mrl_scaled$beta, 4), "\n")
cat("  theta =", round(mrl_scaled$theta, 4), "\n")
cat("  gamma =", round(mrl_scaled$gamma, 4), "\n")
cat("  objective =", round(mrl_scaled$objective, 6), "\n\n")

cat("MRL (log):\n")
cat("  beta =", round(mrl_log$beta, 4), "\n")
cat("  theta =", round(mrl_log$theta, 4), "\n")
cat("  gamma =", round(mrl_log$gamma, 4), "\n")
cat("  objective =", round(mrl_log$objective, 6), "\n\n")

# ==============================================================================
# BÖLÜM 9: SİMÜLASYON KARŞILAŞTIRMASI
# ==============================================================================

cat("=" , rep("=", 60), "\n", sep = "")
cat("SİMÜLASYON KARŞILAŞTIRMASI (n=100, 500 tekrar)\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

n_sim <- 500
n <- 100

results <- data.frame(
  method = character(),
  beta_mean = numeric(),
  beta_mse = numeric(),
  theta_mean = numeric(),
  theta_mse = numeric(),
  gamma_mean = numeric(),
  gamma_mse = numeric(),
  n_valid = numeric()
)

mom_res <- matrix(NA, n_sim, 3)
mrl_scaled_res <- matrix(NA, n_sim, 3)
mrl_log_res <- matrix(NA, n_sim, 3)

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (i in 1:n_sim) {
  x <- rweibull(n, shape = beta_true, scale = theta_true)
  y <- gamma_true + x
  
  # MOM
  tryCatch({
    fit <- weibull_3p_mom(y)
    mom_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
  }, error = function(e) NULL)
  
  # MRL scaled
  tryCatch({
    fit <- weibull_3p_mrl_v2(y, method = "scaled")
    if (!is.na(fit$beta)) {
      mrl_scaled_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
    }
  }, error = function(e) NULL)
  
  # MRL log
  tryCatch({
    fit <- weibull_3p_mrl_v2(y, method = "log")
    if (!is.na(fit$beta)) {
      mrl_log_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
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
      beta_mean = NA, beta_mse = NA,
      theta_mean = NA, theta_mse = NA,
      gamma_mean = NA, gamma_mse = NA,
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
    beta_mean = means[1], beta_mse = mse[1],
    theta_mean = means[2], theta_mse = mse[2],
    gamma_mean = means[3], gamma_mse = mse[3],
    n_valid = n_valid
  )
}

true_vals <- c(beta_true, theta_true, gamma_true)

results <- rbind(
  calc_stats(mom_res, true_vals, "MOM"),
  calc_stats(mrl_scaled_res, true_vals, "MRL_scaled"),
  calc_stats(mrl_log_res, true_vals, "MRL_log")
)

cat("\n")
print(results, digits = 4)

cat("\n\nGerçek değerler: beta =", beta_true, ", theta =", theta_true, ", gamma =", gamma_true, "\n")

# ==============================================================================
# BÖLÜM 10: THRESHOLD SEÇİMİNİN ETKİSİ
# ==============================================================================

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("THRESHOLD SEÇİMİNİN ETKİSİ\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

thresholds <- c(0.10, 0.20, 0.25, 0.30, 0.40, 0.50)
n_sim_thr <- 200

thr_results <- data.frame()

for (thr in thresholds) {
  mrl_res <- matrix(NA, n_sim_thr, 3)
  
  for (i in 1:n_sim_thr) {
    x <- rweibull(n, shape = beta_true, scale = theta_true)
    y <- gamma_true + x
    
    tryCatch({
      fit <- weibull_3p_mrl_v2(y, threshold_prob = thr, method = "scaled")
      if (!is.na(fit$beta)) {
        mrl_res[i, ] <- c(fit$beta, fit$theta, fit$gamma)
      }
    }, error = function(e) NULL)
  }
  
  stats <- calc_stats(mrl_res, true_vals, paste0("thr_", thr))
  stats$threshold <- thr
  thr_results <- rbind(thr_results, stats)
}

cat("Threshold etkisi:\n")
print(thr_results[, c("threshold", "beta_mse", "theta_mse", "gamma_mse", "n_valid")], digits = 4)