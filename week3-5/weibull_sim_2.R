rm(list = ls())
n <- 50
beta <- 1.5
theta <- 1
tekrar <- 1000

# ------------------ Fonksiyonları döngü dışına al ------------------ #

# Üst eksik gamma: Γ(s, z) = upper incomplete gamma
Gamma_upper <- function(s, z) {
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

# f1 ve f2 sistemini tanımla
f_system <- function(par, xz, mrl_x1, mrl_x2) {
  a <- par[1]
  b <- par[2]
  
  # a ve b pozitif olmalı
  if (a <= 0 || b <= 0) {
    return(c(1e6, 1e6))
  }
  
  t <- (xz / a)^b
  
  G1 <- Gamma_upper(1 + 1 / b, t)
  G2 <- Gamma_upper(1 + 2 / b, t)
  
  f1 <- a * exp(t) * G1 - xz - mrl_x1
  
  f2 <- a^2 * exp(t) * G2 -
    2 * xz * a * exp(t) * G1 +
    xz^2 - mrl_x2
  
  c(f1, f2)
}

# Negatif log-likelihood
weibull_nll <- function(par, x) {
  shape <- par[1]  # beta
  scale <- par[2]  # theta
  
  if (shape <= 0 || scale <= 0) return(Inf)
  
  x <- x[x > 0]
  
  n    <- length(x)
  logx <- log(x)
  
  loglik <- n * log(shape) +
    (shape - 1) * sum(logx) -
    n * shape * log(scale) -
    sum((x / scale)^shape)
  
  -loglik
}

# Weibull MLE fonksiyonu
weibull_mle <- function(x, start_shape = 1, start_scale = mean(x)) {
  x <- x[x > 0]
  
  fit <- optim(
    par    = c(start_shape, start_scale),
    fn     = weibull_nll,
    x      = x,
    method = "L-BFGS-B",
    lower  = c(1e-6, 1e-6)
  )
  
  list(
    shape_mle = fit$par[1],
    scale_mle = fit$par[2],
    nll       = fit$value,
    convergence = fit$convergence
  )
}

# Momentler yöntemi
weibull_mom <- function(x) {
  x <- x[x > 0]
  m <- mean(x)
  v <- var(x)
  
  r_hat <- v / m^2
  
  cv2_theoretical <- function(shape) {
    g1 <- gamma(1 + 1 / shape)
    g2 <- gamma(1 + 2 / shape)
    g2 / (g1^2) - 1
  }
  
  objective <- function(shape) {
    (cv2_theoretical(shape) - r_hat)^2
  }
  
  opt <- optimize(objective, interval = c(0.1, 10))
  shape_hat <- opt$minimum
  
  scale_hat <- m / gamma(1 + 1 / shape_hat)
  
  list(
    shape_mom  = shape_hat,
    scale_mom  = scale_hat,
    mean_sample = m,
    var_sample  = v
  )
}

library(nleqslv)

# ------------------ SONUÇLARI TOPLAYACAĞIMIZ YAPI ------------------ #

# Her iterasyonun sonucunu koyacağımız bir liste
results_list <- vector("list", tekrar)

# ------------------ SİMÜLASYON DÖNGÜSÜ ------------------ #
# set.seed(123)

for(i in 1:tekrar){
  Tvals <- stats::rweibull(n, shape = beta, scale = theta)
  Tvals <- sort(Tvals)
#  index <- sample(1:n,floor(n/4))
  index <- 1:floor(n/4)
  Tvals[index] <- Tvals[index] * mean(Tvals) + mean(Tvals)
  x_thr <- quantile(Tvals, probs = 0.25)

  exc <- Tvals[Tvals > x_thr] - x_thr
  mrl_x1 <- mean(exc)
  mrl_x2 <- mean(exc^2)
  
  # MRL yöntemi
  yont_mrl <- nleqslv(
    x = c(0.5, 0.5),
    fn = f_system,
    xz = x_thr,
    mrl_x1 = mrl_x1,
    mrl_x2 = mrl_x2,
    method = "Broyden"
  )
  
  # MLE
  yont_mle <- weibull_mle(Tvals, start_shape = 1.2, start_scale = 0.4)
  
  # MOM
  yont_mom <- weibull_mom(Tvals)
  
  df_yont_mrl <- data.frame(
    iter   = i,
    method = "mrl",
    theta  = yont_mrl$x[1],
    beta   = yont_mrl$x[2]
  )
  
  df_yont_mom <- data.frame(
    iter   = i,
    method = "mom",
    theta  = yont_mom$scale_mom,
    beta   = yont_mom$shape_mom
  )
  
  df_yont_mle <- data.frame(
    iter   = i,
    method = "mle",
    theta  = yont_mle$scale_mle,
    beta   = yont_mle$shape_mle
  )
  
  # Bu iterasyonun üç yöntemi:
  results_list[[i]] <- rbind(df_yont_mrl, df_yont_mle, df_yont_mom)
}

# ------------------ HEPSİNİ BİRLEŞTİR VE ORTALAMAYI AL ------------------ #

# Tüm iterasyonları tek data.frame'de birleştir
results_all <- do.call(rbind, results_list)

# Yönteme göre ortalama theta ve beta
library(dplyr)

ortalama_sonuclar <- results_all %>%
  group_by(method) %>%
  summarise(
    mean_beta  = mean(beta,  na.rm = TRUE),
    var_beta    = var(beta,    na.rm = TRUE),
    
    mean_theta = mean(theta, na.rm = TRUE),
    var_theta   = var(theta,   na.rm = TRUE),
    .groups = "drop"
  )

ortalama_sonuclar <- ortalama_sonuclar |>
  mutate(mse_theta = (mean_theta - theta)^2 + var_theta) |> 
  mutate(mse_beta = (mean_beta - beta)^2 + var_beta) |> 
  mutate(n = n, tekrar = tekrar,beta_gercek = beta, theta_gercek = theta)
ortalama_sonuclar |>
  select(Method= method,mean_beta,var_beta,mse_beta, 
         mean_theta,var_theta, mse_theta,n,tekrar,beta_gercek,theta_gercek) |> 
  as.data.frame()

