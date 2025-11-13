set.seed(123)

theta <- 0.5          # gerçek parametre
R     <- 10000       # tekrar sayısı
n_vec <- c(10, 25, 50, 100, 500)

# Teorik Q1, Q2, Q3 
p_vec <- c(0.25, 0.5, 0.75)
# x_vec <- -theta * log(1 - p_vec)   # −θln(1−p) = F^-1

xbar_results   <- data.frame()
excess_results <- data.frame()

for (n in n_vec) {
  # Her n için 10 000 tane Xbar ve MRL tahmini 
  xbar_vals       <- numeric(R)
  mrl_excess_vals <- matrix(NA_real_, nrow = R, ncol = length(x_vec))
  colnames(mrl_excess_vals) <- paste0("Q", 1:3)  # Q1, Q2, Q3 sütunlarını oluşturuyor
  
  for (r in 1:R) {
    # 1) Örneklem üret: Exp(mean = theta)
    x <- stats::rexp(n, rate = 1/theta)
    x_vec <- quantile(x, probs = c(0.25,0.5,0.75))
    # 2) Parametrik tahmin edici: theta_hat_1 = Xbar
    xbar_vals[r] <- mean(x)
    
    # 3) MRL tabanlı tahmin edici:
    #    MRL_hat(x) = ortalama (X_i - x) | X_i > x
    for (j in seq_along(x_vec)) {
      thr <- x_vec[j]              # eşik: Q1, Q2, Q3
      exc <- x[x > thr] - thr      # t_i - x, t_i > x
      
      if (length(exc) > 0) {
        mrl_excess_vals[r, j] <- mean(exc)
      }
    }
  }
  
  ## Xbar sonuçları (theta_hat_1)
  xbar_results <- rbind(
    xbar_results,
    data.frame(
      n         = n,
      est       = "Xbar",
      mean_est  = mean(xbar_vals),
      var_est   = stats::var(xbar_vals),
      theo_mean = theta,
      theo_var  = theta^2 / n   # Var(Xbar) = theta^2 / n
    )
  )
  
  ## MRL(x) tahmin sonuçları (theta_hat_2(x))
  for (j in seq_along(x_vec)) {
    excess_results <- rbind(
      excess_results,
      data.frame(
        n         = n,
        which_Q   = paste0("Q", j),   # Q1, Q2, Q3
        x_value   = x_vec[j],
        mean_est  = mean(mrl_excess_vals[, j], na.rm = TRUE),
        var_est   = stats::var(mrl_excess_vals[, j], na.rm = TRUE),
        theo_MRL  = theta             # teorik MRL(x) = theta
      )
    )
  }
}

xbar_results$mse <-  xbar_results$var_est + (xbar_results$mean_est - xbar_results$theo_mean)^2 

excess_results$mse <- excess_results$var_est + (excess_results$mean_est - excess_results$theo_MRL)^2

excess_results
xbar_results