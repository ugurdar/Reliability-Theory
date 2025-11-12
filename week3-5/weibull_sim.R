############################################################
## T ~ Weibull(β, θ), shape = β, scale = θ
##
## Amacımız:
##   E[(T - x)^2 | T > x] için
##   (1) teorik formülü kullanarak değerleri hesaplamak,
##   (2) aynı şeyi Monte Carlo simülasyonla ampirik tahmin
##       edip karşılaştırmak.
############################################################

## İstersen temiz başlangıç:
## rm(list = ls())

set.seed(123)

## Weibull parametreleri (GENEL, sen istediğin β, θ'yı koy)
beta  <- 1.5   # shape
theta <- 0.5   # scale

## Simülasyon ayarları
R     <- 10000                   # tekrar sayısı
n_vec <- c(10, 25, 50, 100, 500) # örneklem büyüklükleri

############################################################
## x eşikleri: Weibull'ün teorik Q1, Q2, Q3'ü
##
## F(t) = 1 - exp( - (t/theta)^beta )
##  => Q_p = theta * [ -log(1 - p) ]^(1/beta)
############################################################

p_vec <- c(0.25, 0.50, 0.75)  # Q1, Q2, Q3
x_vec <- theta * (-log(1 - p_vec))^(1 / beta)

############################################################
## Teorik:
##
##   z = (x/theta)^beta
##
##   E[T^r | T > x] = theta^r * exp(z) * Γ(1 + r/β, z)
##
##   E[(T - x)^2 | T > x]
##     = theta^2 * exp(z) * Γ(1 + 2/β, z)
##       - 2 x theta * exp(z) * Γ(1 + 1/β, z)
##       + x^2
##
## R'de upper incomplete gamma:
##   Γ(a, z) = gamma(a) * pgamma(z, shape = a, lower.tail = FALSE)
############################################################

weibull_m2_cond <- function(x, theta, beta) {
  z  <- (x / theta)^beta      # z = (x/theta)^beta
  s1 <- 1 + 1 / beta          # 1 + 1/β
  s2 <- 1 + 2 / beta          # 1 + 2/β
  
  ## Upper incomplete gamma: Γ(s, z)
  G1 <- gamma(s1) * stats::pgamma(z, shape = s1, lower.tail = FALSE)
  G2 <- gamma(s2) * stats::pgamma(z, shape = s2, lower.tail = FALSE)
  
  ## E[T | T > x]  = theta   * exp(z) * Γ(1 + 1/β, z)
  ## E[T^2| T > x] = theta^2 * exp(z) * Γ(1 + 2/β, z)
  ET1 <- theta   * exp(z) * G1
  ET2 <- theta^2 * exp(z) * G2
  
  ## E[(T - x)^2 | T > x] = E[T^2|T>x] - 2x E[T|T>x] + x^2
  m2 <- ET2 - 2 * x * ET1 + x^2
  return(m2)
}

## Her x (Q1, Q2, Q3) için teorik E[(T - x)^2 | T > x] değerleri:
theo_m2_vec <- sapply(x_vec, weibull_m2_cond, theta = theta, beta = beta)
theo_m2_vec
# 1. eleman: Q1 için teorik E[(T - Q1)^2 | T > Q1]
# 2. eleman: Q2 için ...
# 3. eleman: Q3 için ...

############################################################
## Monte Carlo simülasyonu:
##
## Her n ∈ {10, 25, 50, 100, 500} için:
##   R kez:
##     - T_i ~ Weibull(β, θ)
##     - Her x (Q1, Q2, Q3) için ampirik tahmin:
##
##         m2_hat(x) = (1 / #{i : T_i > x}) * Σ_{i : T_i > x} (T_i - x)^2
##
## Sonra m2_hat(x)'lerin ortalamasını ve varyansını hesaplayıp
## teorik değerlerle kıyaslıyoruz.
############################################################

excess2_results <- data.frame()

for (n in n_vec) {
  # Her n için: R x length(x_vec) matris
  m2_excess_vals <- matrix(NA_real_, nrow = R, ncol = length(x_vec))
  colnames(m2_excess_vals) <- paste0("Q", 1:3)  # Q1, Q2, Q3
  
  for (r in 1:R) {
    ## T_i ~ Weibull(β, θ)
    Tvals <- stats::rweibull(n, shape = beta, scale = theta)
    
    for (j in seq_along(x_vec)) {
      x_thr <- x_vec[j]
      
      ## T_i > x koşulunu sağlayanlar için:
      ##   excess = (T_i - x)
      exc <- Tvals[Tvals > x_thr] - x_thr
      
      ## En az bir gözlem varsa conditional ikinci moment:
      ##   m2_hat(x) = ortalama (T_i - x)^2 | T_i > x
      if (length(exc) > 0) {
        m2_excess_vals[r, j] <- mean(exc^2)
      }
    }
  }
  
  ## Her eşik x (Q1, Q2, Q3) için R tekrar üzerinden:
  ##   mean_est = simülasyon ortalaması ≈ teorik E[(T-x)^2 | T>x]
  ##   var_est  = tahmin edicinin simülasyon varyansı
  for (j in seq_along(x_vec)) {
    excess2_results <- rbind(
      excess2_results,
      data.frame(
        n         = n,
        which_Q   = paste0("Q", j),        # "Q1", "Q2", "Q3"
        x_value   = x_vec[j],
        mean_est  = mean(m2_excess_vals[, j], na.rm = TRUE),
        var_est   = stats::var(m2_excess_vals[, j],  na.rm = TRUE),
        theo_val  = theo_m2_vec[j]         # teorik E[(T - x)^2 | T > x]
      )
    )
  }
}

excess2_results
############################################################
## Yorum:
##
## - 'mean_est' sütunu:
##     m2_hat(x) = (1 / #{i : T_i > x}) * Σ_{i : T_i > x} (T_i - x)^2
##   ampirik tahmin edicisinin Monte Carlo ortalamasıdır.
##
## - 'theo_val' sütunu:
##   Genel formüle karşılık gelen teorik değeri gösterir:
##
##   E[(T - x)^2 | T > x]
##     = theta^2 * exp(z) * Γ(1 + 2/beta, z)
##       - 2 * x * theta * exp(z) * Γ(1 + 1/beta, z)
##       + x^2
##   z = (x/theta)^beta.
##
## n büyüdükçe 'mean_est' değerleri 'theo_val'a yaklaşmalı,
## 'var_est' ise küçülmelidir -> tahmin edici tutarlı.
############################################################
