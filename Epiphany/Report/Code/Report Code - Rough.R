# Simulation Approach - Lexie 
# The other helper Functions are from Report Code - Clean.R
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
lines(xgrid, NP_B$NP_m_est, lwd = 2)

sim.test.ll <- function(data, xgrid, h, kernel = "gauss", alpha = 0.05, B = 300) {
  nu0 <- switch(kernel, "gauss" = 1 / (2 * sqrt(pi)), "epa" = 3 / 5, "uni" = 1 / 2)
  mu2 <- switch(kernel, "gauss" = 1, "epa" = 1 / 5, "uni" = 1 / 3)
  
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  # Initialise matrices
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  biasMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  sdMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  
  # Bootstrap
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    
    model_boot <- lm(y_boot ~ poly(x_boot, 4, raw = TRUE))
    coefs_boot <- coef(model_boot)
    
    mhat_boot <- my.ll.smoother(xdat = x_boot, ydat = y_boot, xgrid = xgrid, h = h)
    mddash_boot <- m2(coefs_boot, xgrid)
    sigma_boot <- (sum((y_boot - m(coefs_boot, x_boot))^2) / (n - 5))^(1/2)
    f_boot <- est_density(x_boot, xgrid, h, kernel)
    f_boot[f_boot < 0.01] <- 0.01
    
    bias_val <- bias(h, mddash_boot, mu2)
    sd_val <- variance(n, h, sigma_boot, nu0, f_boot)^(1/2)
    
    bootMat[b, ] <- mhat_boot - bias_val
    sdMat[b, ] <- sd_val
  }
  
  est_mean <- apply(bootMat, 2, mean)
  sd_mean <- apply(sdMat, 2, mean)
  
  z <- qnorm(1 - alpha / 2)
  lower_ci <- est_mean - z * sd_mean
  upper_ci <- est_mean + z * sd_mean
  
  ord <- order(xgrid)
  return(list(
    xgrid = xgrid[ord],
    est = est_mean[ord],
    lower = lower_ci[ord],
    upper = upper_ci[ord]
  ))
}

res_test_bootstrap <- sim.test.v3.ll(as.matrix(onions.data), xgrid = xgrid, h = h_rot_manual, 
                        kernel = "gauss", alpha = 0.05, B = 500)
res_test_bootstrap$lower

which(is.na(est_mean))         # Grid indices with NA in the estimate
which(is.na(sd_mean))          # Grid indices with NA in the SD
which(rowSums(is.na(bootMat)) == ncol(bootMat))  # Bootstrap reps that failed

plot(onions.data$dens, onions.data$yield, type = "n", pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Test Simulation Approach")

# Simulated Approach
# lines(res_bootstrap$xgrid, res_bootstrap$est, col = "purple", lwd = 2)
lines(res_test_analytic$xgrid, res_test_analytic$est, col = "black", lwd = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$lower, col = "purple", lty = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$upper, col = "purple", lty = 2)

sim.test.v2.ll <- function(data, xgrid, h, kernel = "gauss", alpha = 0.05, B = 300) {
  # Kernel constants (for bias correction only)
  mu2 <- switch(kernel, "gauss" = 1, "epa" = 1 / 5, "uni" = 1 / 3)
  
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    
    mhat_boot <- my.ll.smoother(x_boot, y_boot, xgrid, h)
    bootMat[b, ] <- mhat_boot  # no bias correction inside loop
    if (any(is.na(mhat_boot))) {
      cat("Bootstrap", b, "has NA at grid indices:", which(is.na(mhat_boot)), "\n")
    }
  }
  
  # Centre estimate (bias-corrected only once)
  est_mean <- apply(bootMat, 2, mean)
  
  model <- lm(y ~ poly(x, 4, raw = TRUE))
  coefs <- coef(model)
  mddash <- m2(coefs, xgrid)
  bias_val <- bias(h, mddash, mu2)
  
  est_mean_bc <- est_mean - bias_val  # bias-corrected centre
  sd_mean     <- apply(bootMat, 2, sd)
  
  z <- qnorm(1 - alpha / 2)
  
  na_cols <- which(colSums(is.na(bootMat)) > 0)
  cat("NA grid points:", na_cols, "\n")
  
  lower_ci <- est_mean_bc - z * sd_mean
  upper_ci <- est_mean_bc + z * sd_mean
  
  ord <- order(xgrid)
  
  # Handle NA values safely
  keep <- complete.cases(est_mean_bc, sd_mean)
  
  return(list(
    xgrid = xgrid[ord][keep[ord]],
    est   = est_mean_bc[ord][keep[ord]],
    lower = lower_ci[ord][keep[ord]],
    upper = upper_ci[ord][keep[ord]]
  ))
}

ana.test.ll <- function(data, xgrid, h = 1, kernel = "gauss", alpha = 0.05) {
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
  sigma_hat <- (sum((y - m(coefs, x))^2) / (n - 5))^(1/2)
  f <- est_density(x, xgrid, h, kernel)
  f[f < 0.01] <- 0.01 # floor density for stability
  
  bias.ll <- bias(h, mddash, mu2)
  var.ll <- variance(n, h, sigma_hat, nu0, f)
  
  est <- mhat - bias.ll # Bias-corrected estimate
  z <- qnorm(1-alpha/2)
  lower_ci <- est - z * (var.ll)^(1/2) # Compute CI: (estimate - bias) - 2 * sd
  upper_ci <- est + z * (var.ll)^(1/2) # Compute CI: (estimate - bias) + 2 * sd
  
  # Sort by xgrid to ensure lines() function plots properly
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(xgrid = xgrid_sorted, est = est_sorted,
              lower = lower_ci_sorted, upper = upper_ci_sorted))
}

res_test_analytic <- ana.test.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                                  kernel = "gauss", alpha = 0.05)

plot(onions.data$dens, onions.data$yield, type = "n", pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Test Analytical Approach")

# Analytic Approach
lines(res_test_analytic$xgrid, res_test_analytic$est, col = "black", lwd = 2)
lines(res_test_analytic$xgrid, res_test_analytic$lower, col = "green", lty = 2)
lines(res_test_analytic$xgrid, res_test_analytic$upper, col = "green", lty = 2)

# Overall
plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Overall Approach")
lines(res_test_analytic$xgrid, res_test_analytic$est, col = "black", lwd = 2)
lines(res_test_analytic$xgrid, res_test_analytic$lower, col = "purple", lty = 2)
lines(res_test_analytic$xgrid, res_test_analytic$upper, col = "purple", lty = 2)

lines(res_test_bootstrap$xgrid, res_test_bootstrap$lower, col = "gold", lty = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$upper, col = "gold", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("goldenrod1", "purple"), lwd = 2, lty = c(2, 3), bty = "n")

sim.test.v3.ll <- function(data, xgrid, h, kernel = "gauss", alpha = 0.05, B = 300, min_n = 30) {
  # Validate inputs
  if (!is.matrix(data) || ncol(data) != 2) stop("Data must be a matrix with 2 columns (x, y).")
  if (nrow(data) < min_n) stop("Sample size too small for reliable bootstrap.")
  if (h <= 0) stop("Bandwidth h must be positive.")
  
  # Kernel constants
  mu2 <- switch(kernel, "gauss" = 1, "epa" = 1/5, "uni" = 1/3, stop("Unsupported kernel."))
  
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  
  # Initialize bootstrap matrix
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  
  # Bootstrap loop
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    
    mhat_boot <- tryCatch(
      my.ll.smoother(x_boot, y_boot, xgrid, h),
      error = function(e) {
        warning("Error in bootstrap ", b, ": ", e$message)
        return(rep(NA, length(xgrid)))
      }
    )
    
    bootMat[b, ] <- mhat_boot
    if (sum(is.na(mhat_boot)) > 0.1 * length(xgrid)) {
      warning("Bootstrap ", b, " has >10% NA values at grid indices: ", which(is.na(mhat_boot)))
    }
  }
  
  # Check for excessive NA values
  na_cols <- which(colSums(is.na(bootMat)) > 0.5 * B)
  if (length(na_cols) > 0) {
    warning("Grid points with >50% NA values: ", na_cols)
    xgrid <- xgrid[-na_cols]
    bootMat <- bootMat[, -na_cols, drop = FALSE]
  }
  
  # Centre estimate (bias-corrected)
  est_mean <- apply(bootMat, 2, mean, na.rm = TRUE)
  
  # Bias correction (using polynomial for second derivative)
  model <- lm(y ~ splines::bs(x, df = 5)) # More flexible than poly(4)
  coefs <- coef(model)
  mddash <- m2(coefs, xgrid) # Assumes m2 is defined
  bias_val <- bias(h, mddash, mu2) # Assumes bias is defined
  
  est_mean_bc <- est_mean - bias_val
  
  # Percentile-based confidence intervals
  ci <- apply(bootMat, 2, quantile, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  lower_ci <- ci[1, ]
  upper_ci <- ci[2, ]
  
  # Order results and handle remaining NA values
  ord <- order(xgrid)
  keep <- complete.cases(est_mean_bc, lower_ci, upper_ci)
  
  return(list(
    xgrid = xgrid[ord][keep[ord]],
    est   = est_mean_bc[ord][keep[ord]],
    lower = lower_ci[ord][keep[ord]],
    upper = upper_ci[ord][keep[ord]]
  ))
}
