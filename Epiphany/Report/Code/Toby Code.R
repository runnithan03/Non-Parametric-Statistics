install.packages(c("KernSmooth", "locpol", "SemiPar"))

library(KernSmooth)
library(locpol)
library(SemiPar)

my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

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

bias.ci <- function(h, m2, mu2) {
  (1/2) * (h^2) * m2 * mu2
}

sd.ci <- function(n, h, sigma2, nu0, f) {
  variance <- (1 / sqrt(n * h)) * (sigma2) * (nu0 / f)
  return(variance^(1/2))
}

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
  
  bias_val <- bias.ci(h, mddash, mu2)
  sd_val <- sd.ci(n, h, sigma2, nu0, f)
  
  est <- mhat - bias_val # Bias-corrected estimate
  
  z <- qnorm(1-alpha/2)
  lower_ci <- est - z * sd_val # Compute CI: (estimate - bias) - z * sd
  upper_ci <- est + z * sd_val # Compute CI: (estimate - bias) + z * sd
  
  # Sort by xgrid to ensure lines() function plots properly
  ord <- order(xgrid)
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  
  return(list(xgrid = xgrid_sorted, est = est_sorted,
              lower = lower_ci_sorted, upper = upper_ci_sorted))
}

sim.ll <- function(data, xgrid, h, kernel = "gauss", alpha = 0.05, B = 300)

thumbBw_manual <- function(x, y, kernel = "gauss") {
  n <- length(x)
  nu0 <- 1 / (2 * sqrt(pi)) # switch(kernel, "gauss" = 1 / (2 * sqrt(pi)), "epa" = 3 / 5, "uni" = 1 / 2)
  mu2 <- 1 # switch(kernel, "gauss" = 1, "epa"   = 1 / 5, "uni" = 1 / 3)
  
  pilot_fit <- lm(y ~ poly(x, degree = 4, raw = TRUE))
  coefs <- coef(pilot_fit)
  
  sigma2 <- sum((y - m(coefs, x))^2) / (n - 5)
  sum_m_double_sq <- sum(m2(coefs, x)^2)
  
  numerator <- nu0 * sigma2
  denominator <- n * mu2^2 * sum_m_double_sq
  h_rot <- (numerator / denominator)^(1/5)
  return(h_rot)
}

data(onions)
onions.data <- onions[, -3] # remove location
summary(onions.data)

h_rot_manual <- thumbBw_manual(onions.data$dens, onions.data$yield, kernel = "gauss")
h_rot_manual

h_rot <- thumbBw(onions.data$dens, onions.data$yield, deg = 1, kernel = gaussK)
h_rot

h_dpill <- dpill(onions.data$dens, onions.data$yield)

# Simulated Approach - don't work...
res_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                        kernel = "gauss", alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, type = "n", pch = 19, xlab = "Density", 
     ylab = "Yield", main = "Gaussian Simulated Approach")

# lines(res_bootstrap$xgrid, res_bootstrap$est, col = "purple", lwd = 2)
lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("goldenrod1", "purple"), lwd = 2, lty = c(2, 3), bty = "n")

my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

xgrid <- seq(15, 185, by = 0.5)
res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                       kernel = "gauss", alpha = 0.05)

plot(onions.data$dens, onions.data$yield, type = "n", pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Gaussian Analytic Approach Manual")

# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)

plot(onions.data$dens, onions.data$yield, type = "n", pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Gaussian Analytic Approach Automatic")

res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot, 
                       kernel = "gauss", alpha = 0.05)

# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)

sim.test.ll <- function(data, xgrid, h = 1, B = 1000, alpha = 0.05) {
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]

  mhat <- my.ll.smoother(xdat = x, ydat = y, xgrid = xgrid, h = h)
  
  mhat_at_x <- my.ll.smoother(xdat = x, ydat = y, xgrid = x, h = h)
  residuals <- y - mhat_at_x # to be resampled
  
  # Initialise bootstrap matrix
  boots <- matrix(NA_real_, nrow = B, ncol = length(xgrid))
  
  for (b in 1:B) {
    resampled_resid <- sample(residuals, replace = TRUE)
    y_star <- mhat_at_x + resampled_resid # new prediction
    boots[b, ] <- my.ll.smoother(xdat = x, ydat = y_star, xgrid = xgrid, h = h)
  }
  
  lower_ci <- apply(boots, 2, quantile, probs = alpha / 2)
  upper_ci <- apply(boots, 2, quantile, probs = 1 - alpha / 2)
  ord <- order(xgrid)
  
  return(list(
    xgrid = xgrid[ord],
    est = mhat[ord],
    lower = lower_ci[ord],
    upper = upper_ci[ord]
  ))
}


res_test_bootstrap <- sim.test.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                        alpha = 0.05, B = 1000)

plot(onions.data$dens, onions.data$yield, type = "n", pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Gaussian Simulated Approach Manual")

# Simulated Approach
# lines(res_bootstrap$xgrid, res_bootstrap$est, col = "purple", lwd = 2)
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$lower, col = "purple", lty = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$upper, col = "purple", lty = 2)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Simulated Approach")
lines(res_test_bootstrap$xgrid, res_test_bootstrap$lower, col = "purple", lty = 2)
lines(res_test_bootstrap$xgrid, res_test_bootstrap$upper, col = "purple", lty = 2)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Analytic Approach")
res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot, 
                       kernel = "gauss", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)

# Uniform
my.kernel<-function(u){
  u<- 0.5*ifelse(abs(u)<=1, 1,0)
  #u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

Sn <- function(xdat, x, h, j){sum(my.kernel((xdat - x)/h)*(xdat - x)^j)}

vix <- function(xdat, x, h){
  my.kernel((xdat - x)/h)*(Sn(xdat, x, h, 2) - (xdat - x)*Sn(xdat, x, h, 1))}

my.ll.smoother <- function(xdat, ydat, xgrid = xdat, h){
  G <- length(xgrid)
  est <- rep(0, G)
  for (j in 1:G){
    est[j] <- sum(ydat*vix(xdat, xgrid[j], h))/sum(vix(xdat, xgrid[j], h))}
  return(est)}


plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Uniform Analytic Approach")
res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot, 
                       kernel = "uni", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Uniform Analytic Approach")
res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_dpill, 
                       kernel = "uni", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)



x<-my.kernel(0.8)
x

Sn <- function(xdat, x, h, j){sum(my.kernel((xdat - x)/h)*(xdat - x)^j)}

vix <- function(xdat, x, h){
  my.kernel((xdat - x)/h)*(Sn(xdat, x, h, 2) - (xdat - x)*Sn(xdat, x, h, 1))}

my.ll.smoother <- function(xdat, ydat, xgrid = xdat, h){
  G <- length(xgrid)
  est <- rep(0, G)
  for (j in 1:G){
    est[j] <- sum(ydat*vix(xdat, xgrid[j], h))/sum(vix(xdat, xgrid[j], h))}
  return(est)}

# Epanechikov
my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  #u<- dnorm(u)
  u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Epanechikov Analytic Approach")
res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot, 
                       kernel = "epa", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Epanechikov Analytic Approach")
res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_dpill, 
                       kernel = "epa", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "black", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "goldenrod1", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "goldenrod1", lty = 2)