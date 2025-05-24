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
    sd_val <- sd(n, h, sqrt(sigma2_boot), nu0, f_boot)
    
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

sim.test.ll <- function(data, xgrid, h, kernel = "gauss", alpha = 0.05, B = 300) {
  nu0 <- switch(kernel, "gauss" = 1 / (2 * sqrt(pi)), "epa" = 3 / 5, "uni" = 1 / 2)
  mu2 <- switch(kernel, "gauss" = 1, "epa" = 1 / 5, "uni" = 1 / 3)
  
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  # Initial fit for residuals
  mhat <- my.ll.smoother(xdat = x, ydat = y, xgrid = x, h = h)
  residuals <- y - mhat
  
  # Initialize matrix for bootstrap estimates
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  
  # Residual bootstrap
  for (b in seq_len(B)) {
    # Resample residuals
    resid_boot <- sample(residuals, size = n, replace = TRUE)
    y_boot <- mhat + resid_boot
    
    # Fit smoother to bootstrap sample
    mhat_boot <- my.ll.smoother(xdat = x, ydat = y_boot, xgrid = xgrid, h = h)
    
    # Store bootstrap estimate
    bootMat[b, ] <- mhat_boot
  }
  
  # Compute mean and empirical SD from bootstrap
  est_mean <- apply(bootMat, 2, mean)
  sd_empirical <- apply(bootMat, 2, sd)
  
  # Compute CIs
  z <- qnorm(1 - alpha / 2)
  lower_ci <- est_mean - z * sd_empirical
  upper_ci <- est_mean + z * sd_empirical
  
  # Sort results
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est_mean[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(xgrid = xgrid_sorted, est = est_sorted,
              lower = lower_ci_sorted, upper = upper_ci_sorted))
}

res_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                        kernel = "gauss", alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "LL Regression with Analytical and Simulated CIs")

# Simulated Approach
lines(res_bootstrap$xgrid, res_bootstrap$est, col = "purple", lwd = 2)
lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("goldenrod1", "purple"), lwd = 2, lty = c(2, 3), bty = "n")

# Test Approach
res_test_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                        kernel = "gauss", alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "LL Regression with Analytical and Simulated CIs")

lines(res_test_bootstrap$xgrid, res_test_bootstrap$est, col = "purple", lwd = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$lower, col = "purple", lty = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$upper, col = "purple", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("goldenrod1", "purple"), lwd = 2, lty = c(2, 3), bty = "n")

res_bootstrap <- sim.test.ll(onions.data, xgrid = xgrid, h = h_rot, 
                        kernel = "gauss", alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "LL Regression with Analytical and Simulated CIs")

# Simulated Approach
lines(res_bootstrap$xgrid, res_bootstrap$est, col = "purple", lwd = 2)
lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)

# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "goldenrod1", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("goldenrod1", "purple"), lwd = 2, lty = c(2, 3), bty = "n")