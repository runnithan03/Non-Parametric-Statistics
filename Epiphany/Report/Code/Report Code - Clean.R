install.packages(c("KernSmooth", "locpol", "SemiPar"))

library(KernSmooth)
library(locpol)
library(SemiPar)

my.kernel<-function(u){return(dnorm(u))}

Sn <- function(xdat, x, h, j){sum(my.kernel((xdat - x)/h)*(xdat - x)^j)}

vix <- function(xdat, x, h){
  my.kernel((xdat - x)/h)*(Sn(xdat, x, h, 2) - (xdat - x)*Sn(xdat, x, h, 1))}

my.ll.smoother <- function(xdat, ydat, xgrid = xdat, h){
  G <- length(xgrid)
  est <- rep(0, G)
  for (j in 1:G){
    est[j] <- sum(ydat*vix(xdat, xgrid[j], h))/sum(vix(xdat, xgrid[j], h))}
  return(est)}

m <- function(coefs, x) {coefs[1] + coefs[2]*x + coefs[3]*x^2 + coefs[4]*x^3 + coefs[5]*x^4}

m2 <- function(coefs, x) {2 * coefs[3] + 6 * coefs[4]*x + 12 * coefs[5]*x^2}

bias <- function(h, m2, mu2) {(1/2) * (h^2) * m2 * mu2}

var.ll <- function(n, h, sigma, nu0, f) {
  variance <- ((1 / sqrt(n * h)) * (sigma^2) * (nu0 / f))
  return(variance)}

est_density <- function(data, xgrid, h, kernel = "gauss") {
  my.kernel <- switch(kernel,
                      "gauss" = function(x) dnorm(x),
                      "epa"   = function(x) { 3/4 * (1 - x^2) *
                          ifelse(abs(x) <= 1, 1, 0) },
                      "uni"   = function(x) { 0.5 * ifelse(abs(x) <= 1, 1, 0) })
  data <- na.omit(data)
  n <- length(data)
  est <- numeric(length(xgrid))
  denom <- n * h
  for (j in seq_along(xgrid)) {
    est[j] <- sum(my.kernel((data - xgrid[j]) / h)) / denom}
  return(est)}

ana.ll <- function(data, xgrid, h = 1, kernel = "gauss", alpha = 0.05) {
  # Kernel constants
  nu0 <- switch(kernel,"gauss" = 1 / (2 * sqrt(pi)),"epa" = 3 / 5,"uni" = 1 / 2)
  mu2 <- switch(kernel,"gauss" = 1,"epa" = 1 / 5,"uni" = 1 / 3)
  
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  model <- lm(y ~ poly(x, 4, raw = TRUE))   # Fit degree 4 polynomial for derivatives
  coefs <- coef(model)
  mhat <- my.ll.smoother(xdat = x, ydat = y, xgrid = xgrid, h = h)
  
  mddash <- m2(coefs, xgrid)
  sigma2 <- sum((y - m(coefs, x))^2) / (n - 5)
  f <- est_density(x, xgrid, h, kernel)
  
  bias_val <- bias(h, mddash, mu2)
  sd_val <- var.ll(n, h, sqrt(sigma2), nu0, f)^(1/2)
  
  est <- mhat - bias_val # Bias-corrected estimate
  z <- qnorm(1-alpha/2)
  lower_ci <- est - z * sd_val # Compute CI: (estimate - bias) - 2 * sd
  upper_ci <- est + z * sd_val # Compute CI: (estimate - bias) + 2 * sd
  
  # Sort by xgrid to ensure lines() function plots properly
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(xgrid = xgrid_sorted, est = est_sorted,
              lower = lower_ci_sorted, upper = upper_ci_sorted))}

sim.ll <- function(data, xgrid, h, kernel = "gauss", alpha = 0.05, B = 300) {
  nu0 <- switch(kernel,"gauss" = 1 / (2 * sqrt(pi)),"epa" = 3 / 5,"uni" = 1 / 2)
  mu2 <- switch(kernel,"gauss" = 1,"epa"   = 1 / 5,"uni" = 1 / 3)
  
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  # Initialise matrices for bootstrap estimates, bias, and sd over xgrid
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
    mhat_boot <- my.ll.smoother(xdat = x_boot, ydat = y_boot, xgrid = xgrid, h = h)
    mddash_boot <- m2(coefs_boot, xgrid)
    sigma2_boot <- sum((y_boot - m(coefs_boot, x_boot))^2) / (n - 5)
    f_boot <- est_density(x_boot, xgrid, h, kernel)
    
    # Compute bias and sd for this bootstrap sample
    bias_val <- bias(h, mddash_boot, mu2)
    sd_val <- var.ll(n, h, sqrt(sigma2_boot), nu0, f_boot)^(1/2)
    
    # Store the bias-corrected estimate and sd
    bootMat[b, ] <- mhat_boot - bias_val
    sdMat[b, ] <- sd_val}
  
  # Compute means from bootstrap
  est_mean <- apply(bootMat, 2, mean)
  sd_mean <- apply(sdMat, 2, mean)
  
  z <- qnorm(1-alpha/2)
  lower_ci <- est_mean - z * sd_mean # Compute CI: (estimate - bias) +/- 2 * sd
  upper_ci <- est_mean + z * sd_mean # Compute CI: (estimate - bias) +/- 2 * sd
  
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est_mean[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(xgrid = xgrid_sorted, est = est_sorted,
              lower = lower_ci_sorted, upper = upper_ci_sorted))}

thumbBw_manual <- function(x, y, kernel = "gauss") {
  n <- length(x)
  nu0 <- switch(kernel, "gauss" = 1 / (2 * sqrt(pi)), "epa" = 3 / 5, "uni" = 1 / 2)
  mu2 <- switch(kernel, "gauss" = 1, "epa"   = 1 / 5, "uni" = 1 / 3)
  
  pilot_fit <- lm(y ~ poly(x, degree = 4, raw = TRUE))
  coefs <- coef(pilot_fit)
  
  sigma2 <- sum((y - m(coefs, x))^2) / (n - 5)
  sum_m_double_sq <- sum(m2(coefs, x)^2)
  
  numerator <- nu0 * sigma2
  denominator <- mu2^2 * sum_m_double_sq
  h_rot <- (numerator / denominator)^(1/5) * n^(-1/5)
  return(h_rot)
}

data(onions)
onions.data <- onions[, -3] # remove location
summary(onions.data)

h_rot_manual <- thumbBw_manual(onions.data$dens, onions.data$yield, kernel = "gauss")
h_rot_manual

xgrid <- seq(15, 185, by = 0.5)
res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                       kernel = "gauss", alpha = 0.05)
res_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                        kernel = "gauss", alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, type = "n", pch = 19,xlab = "Density", 
     ylab = "Yield", main = "LL Regression with Analytical and Simulated CIs")

# Simulated Approach
# lines(res_bootstrap$xgrid, res_bootstrap$est, col = "purple", lwd = 2)
lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)

# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("goldenrod1", "purple"), lwd = 2, lty = c(2, 3), bty = "n")