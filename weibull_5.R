################################################################################
# MRL HİBRİT YÖNTEMİ - 3 PARAMETRELİ WEİBULL TAHMİNİ
# 
# Notasyon (makale ile uyumlu):
#   θ (theta) : shape (şekil) parametresi
#   β (beta)  : scale (ölçek) parametresi  
#   γ (gamma) : location (konum) parametresi
#
# Yazar: [Sizin adınız]
# Tarih: 2025
################################################################################

rm(list = ls())

# ==============================================================================
# BÖLÜM 1: TEORİK ÇERÇEVE
# ==============================================================================
#
# 3-Parametreli Weibull Dağılımı: X ~ Weibull(θ, β, γ)
# 
# PDF:
#   f(x; θ, β, γ) = (θ/β) * ((x-γ)/β)^(θ-1) * exp(-((x-γ)/β)^θ),  x > γ
#
# CDF:
#   F(x; θ, β, γ) = 1 - exp(-((x-γ)/β)^θ)
#
# Survival:
#   S(x; θ, β, γ) = exp(-((x-γ)/β)^θ)
#
# Dönüşüm: Y = X - γ olduğunda, Y ~ Weibull(θ, β) (2-parametreli standart form)
#
# ==============================================================================
# UNCONDITIONAL MOMENTLER (Tüm veri üzerinden):
# ==============================================================================
#
# Y = X - γ için ham momentler:
#   E[Y^k] = β^k · Γ(1 + k/θ)
#
# X için:
#   E[X]   = γ + β · Γ(1 + 1/θ)
#   Var(X) = β² · [Γ(1 + 2/θ) - Γ(1 + 1/θ)²]
#
# ÖNEMLİ: Weibull'un çarpıklığı (skewness) SADECE θ'ya bağlıdır!
#
#   Skewness = E[(X-μ)³] / Var(X)^(3/2)
#            = [Γ(1+3/θ) - 3·Γ(1+1/θ)·Γ(1+2/θ) + 2·Γ(1+1/θ)³] / [Γ(1+2/θ) - Γ(1+1/θ)²]^(3/2)
#
# Bu özellik sayesinde:
#   1) Empirik skewness'tan θ tahmin edilir (β ve γ'dan bağımsız)
#   2) θ bilinince, varyans'tan β hesaplanır
#   3) θ ve β bilinince, ortalamadan γ hesaplanır
#
# ==============================================================================
# MRL (MEAN RESIDUAL LIFE) MOMENTLERİ:
# ==============================================================================
#
# Tanım: Threshold x₀ verildiğinde, koşullu momentler:
#   m_k(x₀) = E[(X - x₀)^k | X > x₀]
#
# Hesaplama için:
#   y₀ = x₀ - γ     (shifted threshold)
#   t  = (y₀/β)^θ   (normalized threshold)
#
# Üst eksik gamma fonksiyonu:
#   Γ_upper(s, z) = ∫_z^∞ t^(s-1) · e^(-t) dt = Γ(s) · P(X > z) where X ~ Gamma(s,1)
#
# Koşullu beklenen değerler (Y = X - γ için):
#   E[Y^k | Y > y₀] = β^k · exp(t) · Γ_upper(1 + k/θ, t)
#
# MRL'nin 1. momenti (Mean Residual Life):
#   m₁(x₀) = E[X - x₀ | X > x₀] 
#          = E[Y | Y > y₀] - y₀
#          = β · exp(t) · Γ_upper(1 + 1/θ, t) - y₀
#
# MRL'nin 2. momenti:
#   m₂(x₀) = E[(X - x₀)² | X > x₀]
#          = E[Y² | Y > y₀] - 2·y₀·E[Y | Y > y₀] + y₀²
#          = β² · exp(t) · Γ_upper(1 + 2/θ, t) - 2·y₀·β·exp(t)·Γ_upper(1 + 1/θ, t) + y₀²
#
# ==============================================================================
# HİBRİT YAKLAŞIM:
# ==============================================================================
#
# Problem: 3 parametre (θ, β, γ) için 3 MRL momenti kullanmak kötü koşullu
#          bir sistem oluşturuyor (condition number ~ 10^6).
#
# Çözüm: İki aşamalı hibrit yaklaşım
#
# AŞAMA 1 - MOM ile γ tahmini:
#   1a) Empirik skewness hesapla: skew_emp = m₃ / m₂^(3/2)
#   1b) Weibull skewness fonksiyonunu tersle: θ_MOM = skew⁻¹(skew_emp)
#   1c) Varyans'tan β tahmin et: β_MOM = √(Var(X) / [Γ(1+2/θ) - Γ(1+1/θ)²])
#   1d) Ortalama'dan γ tahmin et: γ_MOM = X̄ - β_MOM · Γ(1 + 1/θ_MOM)
#
# AŞAMA 2 - MRL ile θ ve β tahmini (γ sabit tutularak):
#   2a) x₀ = quantile(X, 0.25) seç (threshold)
#   2b) Empirik MRL momentlerini hesapla: m₁_emp, m₂_emp
#   2c) γ = γ_MOM sabit tutarak, (θ, β) için optimize et:
#       minimize [(m₁_theo - m₁_emp)/m₁_emp]² + [(m₂_theo - m₂_emp)/m₂_emp]²
#
# Bu yaklaşım çalışır çünkü:
#   - γ tüm veriden tahmin edilir (daha stabil)
#   - 2 parametre + 2 denklem sistemi iyi koşulludur
#   - MRL kuyruk bilgisini kullanarak θ ve β'yı refine eder
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Üst Eksik Gamma Fonksiyonu
# ------------------------------------------------------------------------------
# Γ_upper(s, z) = ∫_z^∞ t^(s-1) · e^(-t) dt
#              = Γ(s) · P(W > z) where W ~ Gamma(s, 1)
# ------------------------------------------------------------------------------
Gamma_upper <- function(s, z) {
  if (z < 0 || is.na(z) || is.nan(z)) return(gamma(s))
  if (z > 700) return(0)  # Overflow koruması
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

# ------------------------------------------------------------------------------
# Weibull Skewness Fonksiyonu
# ------------------------------------------------------------------------------
# Skewness(θ) = [Γ(1+3/θ) - 3·Γ(1+1/θ)·Γ(1+2/θ) + 2·Γ(1+1/θ)³] / [Γ(1+2/θ) - Γ(1+1/θ)²]^(3/2)
#
# NOT: Bu fonksiyon SADECE θ'ya bağlıdır, β ve γ'dan tamamen bağımsızdır!
# ------------------------------------------------------------------------------
weibull_skewness <- function(theta) {
  if (theta <= 0) return(NA)
  
  g1 <- gamma(1 + 1/theta)
  g2 <- gamma(1 + 2/theta)
  g3 <- gamma(1 + 3/theta)
  
  # 2. ve 3. merkezi momentler (β=1 için, β sadeleşir)
  mu2 <- g2 - g1^2
  mu3 <- g3 - 3*g1*g2 + 2*g1^3
  
  # Skewness
  mu3 / (mu2^(3/2))
}

# ------------------------------------------------------------------------------
# Teorik MRL Momentleri
# ------------------------------------------------------------------------------
# Girdi: x₀ (threshold), θ (shape), β (scale), γ (location)
# Çıktı: c(m₁, m₂) - MRL'nin 1. ve 2. momentleri
#
# Formüller:
#   y₀ = x₀ - γ
#   t  = (y₀/β)^θ
#   m₁ = β · exp(t) · Γ_upper(1+1/θ, t) - y₀
#   m₂ = β² · exp(t) · Γ_upper(1+2/θ, t) - 2·y₀·β·exp(t)·Γ_upper(1+1/θ, t) + y₀²
# ------------------------------------------------------------------------------
calc_mrl_moments <- function(x0, theta, beta, gamma) {
  y0 <- x0 - gamma
  if (y0 <= 0) return(c(NA, NA))
  
  t <- (y0 / beta)^theta
  if (t > 700 || is.na(t) || is.nan(t)) return(c(NA, NA))
  
  exp_t <- exp(t)
  G1 <- Gamma_upper(1 + 1/theta, t)
  G2 <- Gamma_upper(1 + 2/theta, t)
  
  # E[Y | Y > y₀] ve E[Y² | Y > y₀]
  EY1 <- beta * exp_t * G1
  EY2 <- beta^2 * exp_t * G2
  
  # MRL momentleri
  m1 <- EY1 - y0
  m2 <- EY2 - 2*y0*EY1 + y0^2
  
  c(m1, m2)
}

# ------------------------------------------------------------------------------
# Empirik MRL Momentleri
# ------------------------------------------------------------------------------
# Girdi: x (veri vektörü), x₀ (threshold)
# Çıktı: list(m1, m2, n) - empirik momentler ve örneklem büyüklüğü
# ------------------------------------------------------------------------------
calc_emp_mrl <- function(x, x0) {
  exc <- x[x > x0] - x0  # Exceedances (aşımlar)
  n_exc <- length(exc)
  
  if (n_exc < 5) return(list(m1 = NA, m2 = NA, n = 0))
  
  list(
    m1 = mean(exc),       # Örnek ortalaması
    m2 = mean(exc^2),     # Örnek 2. ham momenti
    n = n_exc
  )
}

# ==============================================================================
# ANA FONKSİYON: MRL HİBRİT TAHMİN EDİCİ
# ==============================================================================
#
# Algoritma:
#
# AŞAMA 1 (MOM - γ tahmini için):
#   1. Empirik skewness hesapla: skew = m₃/m₂^(3/2)
#   2. θ'yı bul: weibull_skewness(θ) = skew denklemini çöz
#   3. β'yı hesapla: β = √(s² / [Γ(1+2/θ) - Γ(1+1/θ)²])
#   4. γ'yı hesapla: γ = x̄ - β·Γ(1+1/θ)
#
# AŞAMA 2 (MRL - θ ve β refinement):
#   1. x₀ = Q₀.₂₅(X) threshold'u seç
#   2. Empirik m₁, m₂ hesapla
#   3. γ_MOM sabit tutarak, (θ, β) optimize et
#
# ==============================================================================
weibull_3p_mrl_hybrid <- function(x, threshold_prob = 0.25) {
  
  n <- length(x)
  x_bar <- mean(x)
  x_min <- min(x)
  
  # ============================================================================
  # AŞAMA 1: MOM ile başlangıç tahminleri
  # ============================================================================
  
  # 1a) Empirik merkezi momentler
  m2 <- mean((x - x_bar)^2)  # ≈ Var(X)
  m3 <- mean((x - x_bar)^3)  # 3. merkezi moment
  
  # 1b) Empirik skewness
  skew_emp <- m3 / (m2^(3/2))
  
  # 1c) θ tahmini: weibull_skewness(θ) = skew_emp denklemini çöz
  obj_theta <- function(theta) {
    if (theta <= 0.1) return(1e10)
    (weibull_skewness(theta) - skew_emp)^2
  }
  opt_theta <- optimize(obj_theta, c(0.1, 50))
  theta_mom <- opt_theta$minimum
  
  # 1d) β tahmini: Var(X) = β² · [Γ(1+2/θ) - Γ(1+1/θ)²]
  g1 <- gamma(1 + 1/theta_mom)
  g2 <- gamma(1 + 2/theta_mom)
  var_factor <- g2 - g1^2
  s2 <- var(x) * (n-1)/n  # Population variance
  beta_mom <- sqrt(s2 / var_factor)
  
  # 1e) γ tahmini: E[X] = γ + β·Γ(1+1/θ)
  gamma_mom <- x_bar - beta_mom * g1
  
  # γ'nın geçerliliğini kontrol et
  gamma_fixed <- gamma_mom
  if (is.na(gamma_fixed) || gamma_fixed >= x_min) {
    gamma_fixed <- x_min * 0.9
  }
  
  # ============================================================================
  # AŞAMA 2: MRL ile θ ve β optimizasyonu (γ sabit)
  # ============================================================================
  
  # 2a) Threshold seç
  x0 <- quantile(x, probs = threshold_prob, names = FALSE)
  
  # 2b) Empirik MRL momentleri
  emp <- calc_emp_mrl(x, x0)
  
  if (is.na(emp$m1) || is.na(emp$m2)) {
    # MRL hesaplanamıyorsa MOM sonuçlarını döndür
    return(list(
      theta = theta_mom,
      beta = beta_mom,
      gamma = gamma_mom,
      method = "MOM_only"
    ))
  }
  
  # 2c) Objective fonksiyon: normalize edilmiş kare hata
  #     L(θ,β) = [(m₁_theo - m₁_emp)/m₁_emp]² + [(m₂_theo - m₂_emp)/m₂_emp]²
  objective <- function(par) {
    theta <- par[1]
    beta <- par[2]
    
    if (theta <= 0 || beta <= 0) return(1e10)
    
    m_theo <- calc_mrl_moments(x0, theta, beta, gamma_fixed)
    if (any(is.na(m_theo))) return(1e10)
    
    err1 <- ((m_theo[1] - emp$m1) / emp$m1)^2
    err2 <- ((m_theo[2] - emp$m2) / emp$m2)^2
    
    err1 + err2
  }
  
  # 2d) Çoklu başlangıç noktası ile optimizasyon
  starts <- list(
    c(theta_mom, beta_mom),
    c(theta_mom * 0.8, beta_mom * 1.2),
    c(theta_mom * 1.2, beta_mom * 0.8),
    c(2, sd(x))
  )
  
  best <- NULL
  best_val <- Inf
  
  for (start in starts) {
    res <- tryCatch({
      optim(start, objective, method = "L-BFGS-B",
            lower = c(0.1, 1e-6), 
            upper = c(50, Inf),
            control = list(maxit = 2000))
    }, error = function(e) NULL)
    
    if (!is.null(res) && res$value < best_val) {
      best_val <- res$value
      best <- res
    }
  }
  
  # Sonuç
  if (is.null(best)) {
    return(list(
      theta = theta_mom,
      beta = beta_mom,
      gamma = gamma_mom,
      method = "MOM_fallback"
    ))
  }
  
  list(
    theta = best$par[1],
    beta = best$par[2],
    gamma = gamma_fixed,
    objective = best_val,
    method = "MRL_hybrid"
  )
}

# ==============================================================================
# BÖLÜM 2: YARDIMCI FONKSİYONLAR
# ==============================================================================

# MOM Tahmini (karşılaştırma için)
weibull_3p_mom <- function(x) {
  n <- length(x)
  x_bar <- mean(x)
  
  m2 <- mean((x - x_bar)^2)
  m3 <- mean((x - x_bar)^3)
  skew_emp <- m3 / (m2^(3/2))
  
  obj_theta <- function(theta) {
    if (theta <= 0.1) return(1e10)
    (weibull_skewness(theta) - skew_emp)^2
  }
  opt_theta <- optimize(obj_theta, c(0.1, 50))
  theta_hat <- opt_theta$minimum
  
  g1 <- gamma(1 + 1/theta_hat)
  g2 <- gamma(1 + 2/theta_hat)
  var_factor <- g2 - g1^2
  s2 <- var(x) * (n-1)/n
  beta_hat <- sqrt(s2 / var_factor)
  gamma_hat <- x_bar - beta_hat * g1
  
  list(theta = theta_hat, beta = beta_hat, gamma = gamma_hat)
}

# MLE Tahmini (karşılaştırma için)
weibull_3p_mle <- function(x) {
  n <- length(x)
  x_min <- min(x)
  
  nll <- function(par) {
    theta <- par[1]
    beta <- par[2]
    gamma <- par[3]
    
    if (theta <= 0 || beta <= 0 || gamma >= x_min) return(Inf)
    
    y <- x - gamma
    if (any(y <= 0)) return(Inf)
    
    ll <- n * log(theta) - n * theta * log(beta) + 
      (theta - 1) * sum(log(y)) - sum((y/beta)^theta)
    -ll
  }
  
  mom <- tryCatch(weibull_3p_mom(x), error = function(e) NULL)
  start <- if (!is.null(mom)) c(mom$theta, mom$beta, mom$gamma) else c(2, sd(x), x_min * 0.5)
  
  res <- tryCatch({
    optim(start, nll, method = "L-BFGS-B",
          lower = c(1e-6, 1e-6, -Inf),
          upper = c(Inf, Inf, x_min - 1e-6),
          control = list(maxit = 2000))
  }, error = function(e) NULL)
  
  if (is.null(res) || res$convergence != 0) {
    return(list(theta = NA, beta = NA, gamma = NA))
  }
  
  list(theta = res$par[1], beta = res$par[2], gamma = res$par[3])
}

# ==============================================================================
# BÖLÜM 3: TEST VE KARŞILAŞTIRMA
# ==============================================================================

cat("=" , rep("=", 70), "\n", sep="")
cat("MRL HİBRİT YÖNTEMİ - 3 PARAMETRELİ WEİBULL TAHMİNİ\n")
cat("=" , rep("=", 70), "\n\n", sep="")

# Dataset 1: Naturally Aged Glass (Datsiou & Overend, 2018)
dataset1 <- c(
  24.12, 24.13, 28.52, 29.18, 29.67, 30.48, 32.98, 35.91, 35.92,
  36.38, 37.60, 37.70, 39.71, 49.10, 52.43, 52.46, 52.61, 61.72
)

cat("Dataset 1: Naturally Aged Glass Strength (MPa)\n")
cat("-" , rep("-", 50), "\n", sep="")
cat("n =", length(dataset1), "\n")
cat("Ortalama =", round(mean(dataset1), 3), "\n")
cat("Std. Sapma =", round(sd(dataset1), 3), "\n")
cat("Min =", min(dataset1), ", Max =", max(dataset1), "\n\n")

# MOM
cat("MOM Tahminleri:\n")
mom_fit <- weibull_3p_mom(dataset1)
cat("  θ (shape)    =", round(mom_fit$theta, 4), "\n")
cat("  β (scale)    =", round(mom_fit$beta, 4), "\n")
cat("  γ (location) =", round(mom_fit$gamma, 4), "\n\n")

# MLE
cat("MLE Tahminleri:\n")
mle_fit <- weibull_3p_mle(dataset1)
cat("  θ (shape)    =", round(mle_fit$theta, 4), "\n")
cat("  β (scale)    =", round(mle_fit$beta, 4), "\n")
cat("  γ (location) =", round(mle_fit$gamma, 4), "\n\n")

# MRL Hybrid
cat("MRL HİBRİT Tahminleri:\n")
mrl_fit <- weibull_3p_mrl_hybrid(dataset1, threshold_prob = 0.25)
cat("  θ (shape)    =", round(mrl_fit$theta, 4), "\n")
cat("  β (scale)    =", round(mrl_fit$beta, 4), "\n")
cat("  γ (location) =", round(mrl_fit$gamma, 4), "\n")
cat("  Objective    =", round(mrl_fit$objective, 6), "\n")
cat("  Method       =", mrl_fit$method, "\n\n")

# Makale sonuçları
cat("-" , rep("-", 50), "\n", sep="")
cat("MAKALEDEKİ PITE SONUÇLARI (Table 3):\n")
cat("  ρ = 1.00: θ = 1.376, β = 17.404, γ = 22.794\n")
cat("  MLE:      θ = 1.279, β = 15.907, γ = 23.524\n")

# ==============================================================================
# BÖLÜM 4: SİMÜLASYON ÇALIŞMASI
# ==============================================================================

cat("\n\n")
cat("=" , rep("=", 70), "\n", sep="")
cat("SİMÜLASYON ÇALIŞMASI\n")
cat("Gerçek parametreler: θ=2, β=3, γ=1\n")
cat("=" , rep("=", 70), "\n\n", sep="")

set.seed(123)
theta_true <- 2
beta_true <- 3
gamma_true <- 1

n_sim <- 500
n <- 100

mom_res <- mle_res <- mrl_res <- matrix(NA, n_sim, 3)

cat("Simülasyon çalışıyor (n=100, 500 tekrar)...\n")
pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (i in 1:n_sim) {
  # Veri üret: X = γ + Y, Y ~ Weibull(θ, β)
  y <- rweibull(n, shape = theta_true, scale = beta_true)
  x <- gamma_true + y
  
  tryCatch({
    fit <- weibull_3p_mom(x)
    mom_res[i, ] <- c(fit$theta, fit$beta, fit$gamma)
  }, error = function(e) NULL)
  
  tryCatch({
    fit <- weibull_3p_mle(x)
    if (!is.na(fit$theta)) mle_res[i, ] <- c(fit$theta, fit$beta, fit$gamma)
  }, error = function(e) NULL)
  
  tryCatch({
    fit <- weibull_3p_mrl_hybrid(x)
    if (!is.na(fit$theta)) mrl_res[i, ] <- c(fit$theta, fit$beta, fit$gamma)
  }, error = function(e) NULL)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# İstatistikler
calc_stats <- function(res, name, true_vals) {
  valid <- complete.cases(res)
  res_clean <- res[valid, ]
  
  means <- colMeans(res_clean)
  vars <- apply(res_clean, 2, var)
  bias <- means - true_vals
  mse <- bias^2 + vars
  
  data.frame(
    Method = name,
    theta_mean = means[1], theta_mse = mse[1],
    beta_mean = means[2], beta_mse = mse[2],
    gamma_mean = means[3], gamma_mse = mse[3],
    total_mse = sum(mse),
    n_valid = sum(valid)
  )
}

true_vals <- c(theta_true, beta_true, gamma_true)

results <- rbind(
  calc_stats(mom_res, "MOM", true_vals),
  calc_stats(mle_res, "MLE", true_vals),
  calc_stats(mrl_res, "MRL_hybrid", true_vals)
)

cat("\n\nSONUÇLAR (Gerçek: θ=2, β=3, γ=1):\n")
cat("-" , rep("-", 80), "\n", sep="")
print(results, digits = 4, row.names = FALSE)

cat("\n\nSIRALAMA (Toplam MSE):\n")
ord <- order(results$total_mse)
for (i in ord) {
  cat("  ", i, ". ", results$Method[i], ": MSE = ", round(results$total_mse[i], 4), "\n", sep="")
}