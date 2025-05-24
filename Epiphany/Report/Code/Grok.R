# Load required packages
library(KernSmooth)
library(locpol)
library(SemiPar)

my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

# Auxiliary functions for S_{n,j} and v_i(x)
Sn<- function(xdat, x, h, j){
  snj<- sum(my.kernel((xdat-x)/h)*(xdat-x)^j)
  return(snj)
}
vix<-function(xdat, x, h){
  my.kernel((xdat-x)/h)*( Sn(xdat,x,h,2)-(xdat-x)*Sn(xdat,x,h,1) )
}
# produces the entire vector v_1(x),...v_n(x)

# Actual local linear smoother 
my.ll.smoother<- function(xdat, ydat, xgrid=xdat, h){
  G<- length(xgrid)
  est<-rep(0,G)
  for (j in 1:G){
    est[j]<-
      sum(ydat*vix(xdat,xgrid[j],h))/sum(vix(xdat,xgrid[j],h))
  }
  return(as.data.frame(cbind(xgrid, est)))
}

# Helper functions
m <- function(coefs, x) {
  coefs[1] + coefs[2]*x + coefs[3]*x^2 + coefs[4]*x^3 + coefs[5]*x^4
}

m_double_dash <- function(coefs, x) {
  2 * coefs[3] + 6 * coefs[4]*x + 12 * coefs[5]*x^2
}

bias <- function(h, m2, mu2) {
  (1/2) * (h^2) * m2 * mu2
}

sd <- function(n, h, sigma, nu0, f) {
  ((1 / (n * h)) * (sigma^2) * (nu0 / f))^(1/2)
}

est_density <- function(data, xgrid, h = 1, kernel = "gauss") {
  my.kernel <- switch(kernel,
                      "gauss" = function(x) dnorm(x),
                      "epa"   = function(x) { 3/4 * (1 - x^2) * ifelse(abs(x) <= 1, 1, 0) },
                      "uni"   = function(x) { 0.5 * ifelse(abs(x) <= 1, 1, 0) })
  data <- na.omit(data)
  n <- length(data)
  g <- length(xgrid)
  est <- numeric(g)
  denom <- n * h
  for (j in seq_along(xgrid)) {
    num <- sum(my.kernel((data - xgrid[j]) / h))
    est[j] <- num / denom
  }
  return(est)
}

# Simulation function with custom CI
sim.test.v1.ll <- function(data, xgrid, h = 1, kernel = "gauss", alpha = 0.05, B = 300) {
  # Kernel constants
  nu0 <- switch(kernel,
                "gauss" = 1 / (2 * sqrt(pi)),
                "epa"   = 3 / 5,
                "uni"   = 1 / 2)
  mu2 <- switch(kernel,
                "gauss" = 1,
                "epa"   = 1 / 5,
                "uni"   = 1 / 3)
  
  data <- na.omit(data)
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  # Fit original model
  model <- lm(y ~ poly(x, 4, raw = TRUE))
  coefs <- coef(model)
  mhat <- m(coefs, xgrid)
  mddash <- m_double_dash(coefs, xgrid)
  sigma2 <- sum((y - m(coefs, x))^2) / (n - 5)
  f <- est_density(x, xgrid, h, kernel)
  
  # Compute bias and sd
  bias_val <- bias(h, mddash, mu2)
  sd_val <- sapply(f, function(f_i) sd(n, h, sqrt(sigma2), nu0, f_i))
  
  # Bias-corrected estimate
  est <- mhat - bias_val
  
  # Custom CIs: estimate ± (bias + 2 * sd)
  lower_ci <- est - (bias_val + 2 * sd_val)
  upper_ci <- est + (bias_val + 2 * sd_val)
  
  # Sort by xgrid to ensure smooth plotting
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(
    xgrid = xgrid_sorted,
    est = est_sorted,
    lower = lower_ci_sorted,
    upper = upper_ci_sorted
  ))
}

# Run and plot
data(onions)
onions.data <- onions[, -3]  # Remove location variable
xgrid <- seq(15, 185, by = 0.5)

res <- sim.test.v1.ll(onions.data, xgrid = xgrid, h = h_rot, kernel = "gauss", alpha = 0.05, B = 300)

# Plotting
plot(onions.data$dens, onions.data$yield, pch = 19,
     xlab = "Density", ylab = "Yield",
     main = "4th Degree Polynomial Fit with Custom CIs")
lines(res$xgrid, res$est, col = "purple", lwd = 2)
lines(res$xgrid, res$lower, col = "purple", lty = 3)
lines(res$xgrid, res$upper, col = "purple", lty = 3)

# Simulation function with non-parametric paired bootstrap for CI
sim.test.v2.ll <- function(data, xgrid, h = 1, kernel = "gauss", alpha = 0.05, B = 300) {
  # Kernel constants
  nu0 <- switch(kernel,
                "gauss" = 1 / (2 * sqrt(pi)),
                "epa"   = 3 / 5,
                "uni"   = 1 / 2)
  mu2 <- switch(kernel,
                "gauss" = 1,
                "epa"   = 1 / 5,
                "uni"   = 1 / 3)
  
  data <- na.omit(data)
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  # Storage matrices for bootstrap estimates, bias, and sd over xgrid
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  biasMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  sdMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  
  # Non-parametric paired bootstrap
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    model_boot <- lm(y_boot ~ poly(x_boot, 4, raw = TRUE))
    coefs_boot <- coef(model_boot)
    mhat_boot <- m(coefs_boot, xgrid)
    mddash_boot <- m_double_dash(coefs_boot, xgrid)
    sigma2_boot <- sum((y_boot - m(coefs_boot, x_boot))^2) / (n - 5)
    f_boot <- est_density(x_boot, xgrid, h, kernel)
    
    # Compute bias and sd for this bootstrap sample
    bias_val <- bias(h, mddash_boot, mu2)
    sd_val <- sd(n, h, sqrt(sigma2_boot), nu0, f_boot)
    
    # Store the bias-corrected estimate, bias, and sd
    bootMat[b, ] <- mhat_boot - bias_val
    biasMat[b, ] <- bias_val
    sdMat[b, ] <- sd_val
  }
  
  # Compute means from bootstrap
  est_mean <- apply(bootMat, 2, mean)
  bias_mean <- apply(biasMat, 2, mean)
  sd_mean <- apply(sdMat, 2, mean)
  
  z <- qnorm(1-alpha/2)
  # Custom CIs: estimate ± (bias + 2 * sd) using bootstrapped means
  lower_ci <- (est_mean - bias_mean) - z * sd_mean
  upper_ci <- (est_mean - bias_mean) + z * sd_mean
  
  # Sort by xgrid to ensure smooth plotting
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est_mean[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(
    xgrid = xgrid_sorted,
    est = est_sorted,
    lower = lower_ci_sorted,
    upper = upper_ci_sorted
  ))
}
library(locpol)
h_rot <- thumbBw(onions.data$dens, onions.data$yield, deg = 1, kernel = gaussK)
res <- sim.test.v2.ll(onions.data, xgrid = xgrid, h = h_rot_manual, kernel = "gauss", alpha = 0.05, B = 300)

# Plotting
plot(onions.data$dens, onions.data$yield, pch = 19,
     xlab = "Density", ylab = "Yield",
     main = "Simulation Approach")
lines(res$xgrid, res$est, col = "purple", lwd = 2)
lines(res$xgrid, res$lower, col = "purple", lty = 3)
lines(res$xgrid, res$upper, col = "purple", lty = 3)

# Simulation function with non-parametric paired bootstrap for CI
analytical.ll <- function(data, xgrid, h = 1, kernel = "gauss", alpha = 0.05, B = 300) {
  # Kernel constants
  nu0 <- switch(kernel,
                "gauss" = 1 / (2 * sqrt(pi)),
                "epa"   = 3 / 5,
                "uni"   = 1 / 2)
  mu2 <- switch(kernel,
                "gauss" = 1,
                "epa"   = 1 / 5,
                "uni"   = 1 / 3)
  
  data <- na.omit(data)
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  # Storage matrices for bootstrap estimates, bias, and sd over xgrid
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  biasMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  sdMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  
  # Non-parametric paired bootstrap
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    model_boot <- lm(y_boot ~ poly(x_boot, 4, raw = TRUE))
    coefs_boot <- coef(model_boot)
    mhat_boot <- my.ll.smoother(xdat = x_boot, ydat = y_boot, xgrid = xgrid, h = h)$est
    mddash_boot <- m_double_dash(coefs_boot, xgrid)
    sigma2_boot <- sum((y_boot - m(coefs_boot, x_boot))^2) / (n - 5)
    # print(length(f_boot))
    f_boot <- est_density(x_boot, xgrid, h, kernel)
    
    # Compute bias and sd for this bootstrap sample
    bias_val <- bias(h, mddash_boot, mu2)
    sd_val <- sapply(f_boot, function(f_i) sd(n, h, sqrt(sigma2_boot), nu0, f_i))
    
    # Store the bias-corrected estimate, bias, and sd
    # print(mhat_boot)
    # print(bias_val)
    bootMat[b, ] <- mhat_boot - bias_val
    biasMat[b, ] <- bias_val
    sdMat[b, ] <- sd_val
  }
  
  # Compute means from bootstrap
  est_mean <- apply(bootMat, 2, mean)
  bias_mean <- apply(biasMat, 2, mean)
  sd_mean <- apply(sdMat, 2, mean)
  
  # Custom CIs: estimate ± (bias + 2 * sd) using bootstrapped means
  z <- qnorm(1-alpha/2)
  
  lower_ci <- (est_mean - bias_mean) + z * sd_mean
  upper_ci <- (est_mean + bias_mean) + z * sd_mean
  
  # Sort by xgrid to ensure smooth plotting
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est_mean[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(
    xgrid = xgrid_sorted,
    est = est_sorted,
    lower = lower_ci_sorted,
    upper = upper_ci_sorted
  ))
}

h_rot <- thumbBw(onions.data$dens, onions.data$yield, deg = 1, kernel = gaussK)
res <- sim.ll(onions.data, xgrid = xgrid, h = h_rot, kernel = "gauss", alpha = 0.05, B = 300)

# Plotting
plot(onions.data$dens, onions.data$yield, pch = 19,
     xlab = "Density", ylab = "Yield",
     main = "Simulated Approach")
lines(res$xgrid, res$est, col = "purple", lwd = 2)
lines(res$xgrid, res$lower, col = "purple", lty = 3)
lines(res$xgrid, res$upper, col = "purple", lty = 3)

library(SemiPar)
library(locpol)
data(onions)
onions.data <- onions[, -3]
x <- onions.data$dens
y <- onions.data$yield

# Manual computation with intermediate outputs
thumbBw_manual_debug <- function(x, y, kernel = "gauss") {
  n <- length(x)
  nu0 <- switch(kernel, "gauss" = 1 / (2 * sqrt(pi)), "epa" = 3 / 5, "uni" = 1 / 2)
  mu2 <- switch(kernel, "gauss" = 1, "epa" = 1 / 5, "uni" = 1 / 3)
  
  pilot_fit <- lm(y ~ poly(x, degree = 4, raw = TRUE))
  coefs <- coef(pilot_fit)
  
  sigma2 <- sum((y - m(coefs, x))^2) / (n - 5)
  sum_m_double_sq <- sum(m2(coefs, x)^2)
  
  numerator <- nu0 * sigma2
  denominator <- mu2^2 * sum_m_double_sq
  h_rot <- (numerator / denominator)^(1/5) * n^(-1/5)
  
  cat("Manual - n:", n, "\n")
  cat("Manual - nu0:", nu0, "\n")
  cat("Manual - mu2:", mu2, "\n")
  cat("Manual - sigma2:", sigma2, "\n")
  cat("Manual - sum_m_double_sq:", sum_m_double_sq, "\n")
  cat("Manual - numerator:", numerator, "\n")
  cat("Manual - denominator:", denominator, "\n")
  cat("Manual - h_rot:", h_rot, "\n")
  
  return(h_rot)
}

# Run manual
h_rot_manual <- thumbBw_manual_debug(onions.data$dens, onions.data$yield, kernel = "gauss")

# Run thumbBw
h_rot <- thumbBw(onions.data$dens, onions.data$yield, deg = 1, kernel = gaussK)
cat("thumbBw - h_rot:", h_rot, "\n")

thumbBw_manual <- function(x, y, kernel = "gauss") {
  # Input validation
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) {
    stop("x and y must be numeric vectors of the same length")
  }
  n <- length(x)
  if (n < 5) stop("Sample size must be at least 5 for a 4th-degree polynomial fit")
  
  # Scale x
  x_mean <- mean(x)
  x_sd <- sd(x)
  x_scaled <- (x - x_mean) / x_sd
  
  # Kernel constants
  nu0 <- switch(kernel, "gauss" = 1 / (2 * sqrt(pi)), "epa" = 3 / 5, "uni" = 1 / 2)
  mu2 <- switch(kernel, "gauss" = 1, "epa" = 1 / 5, "uni" = 1 / 3)
  
  # Fit polynomial on scaled data
  pilot_fit <- lm(y ~ poly(x_scaled, degree = 4, raw = TRUE))
  coefs <- coef(pilot_fit)
  
  # Compute sigma^2
  fitted_scaled <- m(coefs, x_scaled)
  sigma2 <- sum((y - fitted_scaled)^2) / (n - 5)
  
  # Compute sum of squared second derivatives
  sum_m_double_sq <- sum(m2(coefs, x_scaled)^2)
  
  # Compute h_ROT on scaled data
  numerator <- nu0 * sigma2
  denominator <- mu2^2 * sum_m_double_sq # + 1e-10
  h_rot_scaled <- (numerator / denominator)^(1/5) * n^(-1/5)
  
  # Fit polynomial on scaled data
  pilot_fit <- lm(y ~ poly(x, degree = 4, raw = TRUE))
  coefs <- coef(pilot_fit)
  
  # Compute sigma^2
  fit <- m(coefs, x)
  sigma2 <- sum((y - fit)^2) / (n - 5)
  sum_m_double_sq <- sum(m2(coefs, x)^2)
  
  numerator <- nu0 * sigma2
  denominator <- mu2^2 * sum_m_double_sq # + 1e-10
  h_rot <- (numerator / denominator)^(1/5) * n^(-1/5)
  
  # Debug output
  cat("Manual - n:", n, "\n")
  cat("Manual - nu0:", nu0, "\n")
  cat("Manual - mu2:", mu2, "\n")
  cat("Manual - sigma2:", sigma2, "\n")
  cat("Manual - sum_m_double_sq:", sum_m_double_sq, "\n")
  cat("Manual - numerator:", numerator, "\n")
  cat("Manual - denominator:", denominator, "\n")
  cat("Manual - h_rot (scaled):", h_rot_scaled, "\n")
  cat("Manual - h_rot (rescaled):", h_rot, "\n")
  
  return(h_rot)
}

# Test
library(SemiPar)
library(locpol)
data(onions)
onions.data <- onions[, -3]
h_rot_manual <- thumbBw_manual(onions.data$dens, onions.data$yield, kernel = "gauss")
h_rot <- thumbBw(onions.data$dens, onions.data$yield, deg = 1, kernel = gaussK)
cat("thumbBw - h_rot:", h_rot, "\n")
