my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
  return(u)
}

sd_estimator <- function(x, y, coefs) {
  n <- length(x)
  
  # Compute m^(xi) using the polynomial
  mhat_vals <- m(coefs, x)
  
  # Residuals and sigma²
  residuals_sq <- (y - mhat_vals)^2
  sigma2_hat <- sum(residuals_sq) / (n - 5)
  
  # Return constant function
  function(x_grid) {
    rep(sigma2_hat, length(x_grid))
  }
}

# Local linear regression estimator
est.lin_reg <- function(xdat, ydat, xgrid = xdat, h = 1, alpha = 0.05, sig_2) {
  nu0 <- 1 / (2 * sqrt(pi))
  n <- length(xdat)
  g <- length(xgrid)
  
  # Sort xgrid explicitly
  xgrid_sorted <- sort(xgrid)
  
  m_est <- m_lower <- m_upper <- rep(0, g)
  denom <- n * h
  
  for (j in 1:g) {
    num <- sum(my.kernel((xdat - xgrid_sorted[j]) / h))
    est <- num / denom
    m_est[j] <- sum(ydat * my.kernel((xgrid_sorted[j] - xdat) / h)) / sum(my.kernel((xgrid_sorted[j] - xdat) / h))
    m_upper[j] <- m_est[j] + qnorm(1 - alpha / 2) * sqrt(nu0 * sig_2 / (est * n * h))
    m_lower[j] <- m_est[j] - qnorm(1 - alpha / 2) * sqrt(nu0 * sig_2 / (est * n * h))
  }
  
  return(list(est = m_est, lower = m_lower, upper = m_upper, xgrid = xgrid_sorted))
}

NP_bootstrap_m <- function(y, x, xgrid, B, h, sig_2) {
  G <- length(xgrid)
  
  # Estimate mean
  m <- est.lin_reg(x, y, xgrid, h, sig_2 = sig_2)
  m.est <- m$est
  xgrid_sorted <- m$xgrid  # use sorted grid from inside est.lin_reg
  
  m.hat_estimates <- matrix(0, nrow = B, ncol = G)
  n.boot <- length(y)
  
  for (j in 1:B) {
    indices <- sample(1:n.boot, replace = TRUE)
    x.boot <- x[indices]
    y.boot <- y[indices]
    
    boot_result <- est.lin_reg(x.boot, y.boot, xgrid, h, sig_2 = sig_2)
    m.hat_estimates[j, ] <- boot_result$est
  }
  
  return(list(NP_m_est = m.est,
              NP_bootstrap_est = m.hat_estimates,
              xgrid = xgrid_sorted))
}

# Set up grid for onions density
xgrid <- seq(min(onions.data$dens), max(onions.data$dens), length.out = 100)

poly_fit <- lm(yield ~ poly(dens, 4, raw = TRUE), data = onions.data)
coefs <- coef(poly_fit)

# Estimate conditional variance σ²(x)
Sig.2 <- sd_estimator(onions.data$dens, onions.data$yield, coefs)(0)

# Local linear regression estimate using custom bandwidth
M_onions.m <- est.lin_reg(onions.data$dens, onions.data$yield, xgrid = xgrid, h = 0.424, sig_2 = Sig.2)

# Bootstrap estimates
# SP_B <- SP_bootstrap_m(y = onions.data$yield, x = onions.data$dens, xgrid = xgrid, B = 1000, h = 1.37, sig_2 = Sig.2)

NP_B <- NP_bootstrap_m(y = onions.data$yield, x = onions.data$dens, xgrid = xgrid, B = 1000, h = 1.37, sig_2 = Sig.2)

plot(xgrid, NP_B$NP_m_est, type = "n", xlab = "Density", ylab = "Yield", main = "Non-Parametric")
for (i in 1:1000) {
  lines(xgrid, NP_B$NP_bootstrap_est[i, ], col = "deepskyblue4")
}
lines(xgrid, NP_B$NP_m_est, lwd = 2)

# Confidence Intervals
G <- length(xgrid)

NP_Quantiles <- matrix(0, 2, G)
for (j in 1:G) {
  NP_Quantiles[1, j] <- quantile(NP_B$NP_bootstrap_est[, j], probs = 0.975)
  NP_Quantiles[2, j] <- quantile(NP_B$NP_bootstrap_est[, j], probs = 0.025)
}

# Plot the CI comparisons
plot(xgrid, NP_B$NP_m_est, type = "n", xlab = "Density", ylab = "Yield")
lines(xgrid, NP_Quantiles[1, ], lwd = 2, col = "deepskyblue")
lines(xgrid, NP_Quantiles[2, ], lwd = 2, col = "deepskyblue")
# lines(xgrid, SP_Quantiles[1, ], lwd = 2, col = "deepskyblue4")
# lines(xgrid, SP_Quantiles[2, ], lwd = 2, col = "deepskyblue4")

SP_bootstrap_m <- function(y, x, xgrid, B, h, sig_2) {
  G <- length(xgrid)
  
  m <- est.lin_reg(x, y, xgrid, h = h, sig_2 = sig_2)
  m.est <- m$est
  e.hat <- y - m.est
  n.boot <- length(y)
  
  m.hat_estimates <- matrix(0, nrow = B, ncol = G)
  
  for (j in 1:B) {
    indices <- sample(1:n.boot, replace = TRUE)
    e.boot <- e.hat[indices]  # bootstrap residuals
    y.boot <- m.est + e.boot  # bootstrap y values from fixed-x model
    
    # Estimate m_LL on bootstrap sample
    m.full <- est.lin_reg(x, y.boot, xgrid, h, sig_2)
    m.hat_estimates[j, ] <- m.full$est
  }
  
  return(list(SP_m_est = m.est, SP_bootstrap_est = m.hat_estimates))
}

# Plot both bootstrap distributions
par(mfrow = c(1, 2))

plot(xgrid, SP_B$SP_m_est, type = "n", xlab = "Density", ylab = "Yield", main = "Semi-Parametric")
for (i in 1:1000) {
  lines(xgrid, SP_B$SP_bootstrap_est[i, ], col = "deepskyblue4")
}
lines(xgrid, SP_B$SP_m_est, lwd = 2)


SP_Quantiles <- matrix(0, 2, G)
for (j in 1:G) {
  SP_Quantiles[1, j] <- quantile(SP_B$SP_bootstrap_est[, j], probs = 0.975)
  SP_Quantiles[2, j] <- quantile(SP_B$SP_bootstrap_est[, j], probs = 0.025)
}