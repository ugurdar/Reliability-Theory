rm(list = ls())
n <- 50
beta <- 1.5
theta <- 1
mu <- 0.5  # 3. parametre: konum parametresi (threshold)
tekrar <- 1000

# ------------------ Fonksiyonları döngü dışına al ------------------ #

# Üst eksik gamma: Γ(s, z) = upper incomplete gamma
Gamma_upper <- function(s, z) {
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

# 3-parametreli Weibull için MRL denklem sistemi
# sim_2'deki mantık: TEK threshold, ÜÇ moment (k=1,2,3)
f_system_3p <- function(par, xz, mrl_x1, mrl_x2, mrl_x3) {
  a <- par[1]   # scale (theta)
  b <- par[2]   # shape (beta)
  mu <- par[3]  # location (mu)
  
  # Parametreler uygun aralıkta olmalı
  if (a <= 0 || b <= 0 || mu >= xz) {
    return(c(1e6, 1e6, 1e6))
  }
  
  # 3-parametreli Weibull için shifted threshold
  t <- ((xz - mu) / a)^b
  
  # Üst eksik gamma fonksiyonları
  G1 <- Gamma_upper(1 + 1 / b, t)
  G2 <- Gamma_upper(1 + 2 / b, t)
  G3 <- Gamma_upper(1 + 3 / b, t)
  
  # f1: E[X - x | X > x] denklemi (1. moment)
  f1 <- a * exp(t) * G1 - (xz - mu) - mrl_x1
  
  # f2: E[(X - x)^2 | X > x] denklemi (2. moment)
  f2 <- a^2 * exp(t) * G2 -
    2 * (xz - mu) * a * exp(t) * G1 +
    (xz - mu)^2 - mrl_x2
  
  # f3: E[(X - x)^3 | X > x] denklemi (3. moment)
  f3 <- a^3 * exp(t) * G3 -
    3 * (xz - mu) * a^2 * exp(t) * G2 +
    3 * (xz - mu)^2 * a * exp(t) * G1 -
    (xz - mu)^3 - mrl_x3
  
  
  
  c(f1, f2, f3)
}

# 3-parametreli Weibull için negatif log-likelihood
weibull_3p_nll <- function(par, x) {
  shape <- par[1]  # beta
  scale <- par[2]  # theta
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

# 3-parametreli Weibull MLE fonksiyonu
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

# 3-parametreli Weibull için Momentler yöntemi
weibull_3p_mom <- function(x) {
  x <- x[x > 0]
  
  m1 <- mean(x)
  m2 <- mean(x^2)
  m3 <- mean(x^3)
  
  # Merkezi momentler
  mu2 <- m2 - m1^2
  mu3 <- m3 - 3 * m1 * m2 + 2 * m1^3
  
  # Çarpıklık katsayısı
  skewness <- mu3 / (mu2^(3/2))
  
  # Shape parametresini çarpıklıktan tahmin et
  skewness_theoretical <- function(shape) {
    g1 <- gamma(1 + 1 / shape)
    g2 <- gamma(1 + 2 / shape)
    g3 <- gamma(1 + 3 / shape)
    
    numerator <- g3 - 3 * g1 * g2 + 2 * g1^3
    denominator <- (g2 - g1^2)^(3/2)
    
    numerator / denominator
  }
  
  objective <- function(shape) {
    if (shape <= 0) return(Inf)
    (skewness_theoretical(shape) - skewness)^2
  }
  
  opt <- optimize(objective, interval = c(0.1, 10))
  shape_hat <- opt$minimum
  
  # Scale ve location parametrelerini tahmin et
  g1 <- gamma(1 + 1 / shape_hat)
  g2 <- gamma(1 + 2 / shape_hat)
  
  # Varyans denkleminden scale
  scale_hat <- sqrt(mu2 / (g2 - g1^2))
  
  # Ortalama denkleminden location
  location_hat <- m1 - scale_hat * g1
  
  list(
    shape_mom = shape_hat,
    scale_mom = scale_hat,
    location_mom = location_hat,
    mean_sample = m1,
    var_sample = mu2,
    skewness_sample = skewness
  )
}

library(nleqslv)

# ------------------ SONUÇLARI TOPLAYACAĞIMIZ YAPI ------------------ #

# Her iterasyonun sonucunu koyacağımız bir liste
results_list <- vector("list", tekrar)

# ------------------ SİMÜLASYON DÖNGÜSÜ ------------------ #
# set.seed(123)

for(i in 1:tekrar){
  # 3-parametreli Weibull dağılımından veri üret
  Tvals <- mu + stats::rweibull(n, shape = beta, scale = theta)
  Tvals <- sort(Tvals)
  
  # Aykırı değerler ekle
  index <- 1:floor(n/4)
  Tvals[index] <- Tvals[index] * mean(Tvals) + mean(Tvals)
  
  # MRL yöntemi (3 parametreli)
  # sim_2'deki gibi: TEK threshold, ÜÇ MOMENT
  
  # Threshold değeri
  x_thr <- quantile(Tvals, probs = 0.25)
  
  # Exceedances: threshold üzerindeki değerlerin artık yaşam süreleri
  exc <- Tvals[Tvals > x_thr] - x_thr
  mrl_x1 <- mean(exc)       # E[X - x | X > x]
  mrl_x2 <- mean(exc^2)     # E[(X - x)^2 | X > x]
  mrl_x3 <- mean(exc^3)     # E[(X - x)^3 | X > x]
  
  # 3 parametreyi 3 moment denklemiyle tahmin et
  yont_mrl <- tryCatch({
    nleqslv(
      x = c(sd(Tvals), 1.5, min(Tvals) * 0.5),  # başlangıç: theta, beta, mu
      fn = f_system_3p,
      xz = x_thr,
      mrl_x1 = mrl_x1,
      mrl_x2 = mrl_x2,
      mrl_x3 = mrl_x3,
      method = "Broyden",
      control = list(maxit = 500)
    )
  }, error = function(e) {
    list(x = c(NA, NA, NA), termcd = 999)
  })
  
  # MLE
  yont_mle <- tryCatch({
    weibull_3p_mle(Tvals, start_shape = 1.2, start_scale = 0.4, start_location = 0.1)
  }, error = function(e) {
    list(scale_mle = NA, shape_mle = NA, location_mle = NA)
  })
  
  # MOM
  yont_mom <- tryCatch({
    weibull_3p_mom(Tvals)
  }, error = function(e) {
    list(scale_mom = NA, shape_mom = NA, location_mom = NA)
  })
  
  # MRL sonuçlarını kontrol et (sadece yakınsayan sonuçları kullan)
  if (!is.null(yont_mrl$termcd) && yont_mrl$termcd == 1) {
    # termcd == 1: Başarılı yakınsama
    df_yont_mrl <- data.frame(
      iter = i,
      method = "mrl",
      theta = yont_mrl$x[1],
      beta = yont_mrl$x[2],
      mu = yont_mrl$x[3]
    )
  } else {
    # Yakınsamadı - NA kaydet
    df_yont_mrl <- data.frame(
      iter = i,
      method = "mrl",
      theta = NA,
      beta = NA,
      mu = NA
    )
  }
  
  df_yont_mom <- data.frame(
    iter = i,
    method = "mom",
    theta = yont_mom$scale_mom,
    beta = yont_mom$shape_mom,
    mu = yont_mom$location_mom
  )
  
  df_yont_mle <- data.frame(
    iter = i,
    method = "mle",
    theta = yont_mle$scale_mle,
    beta = yont_mle$shape_mle,
    mu = yont_mle$location_mle
  )
  
  # Bu iterasyonun üç yöntemi:
  results_list[[i]] <- rbind(df_yont_mrl, df_yont_mle, df_yont_mom)
}

# ------------------ HEPSİNİ BİRLEŞTİR VE ORTALAMAYI AL ------------------ #

# Tüm iterasyonları tek data.frame'de birleştir
results_all <- do.call(rbind, results_list)

# Yönteme göre ortalama theta, beta ve mu
library(dplyr)

ortalama_sonuclar <- results_all %>%
  group_by(method) %>%
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    var_beta = var(beta, na.rm = TRUE),
    
    mean_theta = mean(theta, na.rm = TRUE),
    var_theta = var(theta, na.rm = TRUE),
    
    mean_mu = mean(mu, na.rm = TRUE),
    var_mu = var(mu, na.rm = TRUE),
    .groups = "drop"
  )

ortalama_sonuclar <- ortalama_sonuclar |>
  mutate(mse_theta = (mean_theta - theta)^2 + var_theta) |>
  mutate(mse_beta = (mean_beta - beta)^2 + var_beta) |>
  mutate(mse_mu = (mean_mu - mu)^2 + var_mu) |>
  mutate(n = n, tekrar = tekrar, beta_gercek = beta, theta_gercek = theta, mu_gercek = mu)

ortalama_sonuclar |>
  select(Method = method,
         mean_beta, var_beta, mse_beta,
         mean_theta, var_theta, mse_theta,
         mean_mu, var_mu, mse_mu,
         n, tekrar, beta_gercek, theta_gercek, mu_gercek) |>
  as.data.frame()
