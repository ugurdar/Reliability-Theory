
n <- 1000000
beta <- 1.5
theta <- 1
Tvals <- stats::rweibull(n, shape = beta, scale = theta)
x_thr <- quantile(Tvals, probs = c(0.25))
p_vec <- c(0.25) 
x_vec <- theta * (-log(1 - p_vec))^(1 / beta)
exc <- Tvals[Tvals > x_thr] - x_thr
mrl_x1 <- mean(exc)
mrl_x2 <- mean(exc^2)
# Üst eksik gamma: Γ(s, z) = upper incomplete gamma
Gamma_upper <- function(s, z) {
  gamma(s) * pgamma(z, shape = s, lower.tail = FALSE)
}

# f1 ve f2 sistemini tanımla
f_system <- function(par, x) {
  a <- par[1]
  b <- par[2]
  
  # a ve b pozitif olmalı (Weibull ölçek/şekil gibi düşünüyorsan)
  if (a <= 0 || b <= 0) {
    return(c(1e6, 1e6))  # saçma büyük değer ver, çözücü bu bölgeden kaçar
  }
  
  t <- (x / a)^b
  
  G1 <- Gamma_upper(1 + 1 / b, t)
  G2 <- Gamma_upper(1 + 2 / b, t)
  
  f1 <- a * exp(t) * G1 - x
  
  f2 <- a^2 * exp(t) * G2 -
    2 * x * a * exp(t) * G1 +
    x^2
  
  c(f1, f2)
}

# Çözüm fonksiyonu: verilen x için a ve b'yi bul
solve_ab <- function(x, a0 = 1, b0 = 1.5) {
  # nleqslv paketi:
  # install.packages("nleqslv") gerekebilir
  library(nleqslv)
  
  nleqslv(
    x = c(a0, b0),
    fn = f_system,
    x = x,          # fn içindeki x argümanı
    method = "Broyden"
  )
}

anakutle <- f_system(c(theta,beta),x_vec)
print(anakutle)
orneklem <- c(mrl_x1,mrl_x2)
print(orneklem)
## ÖRNEK KULLANIM
# x_val <- 1  # istediğin bir x değeri
# res <- solve_ab(x_val, a0 = 0.5, b0 = 1.2)
# 
# res$x   # burası (a, b) çözümünü verir
# res$fvec  # kontrol için f1(a,b) ve f2(a,b) değerleri
