# EDA
library(locpol)
library(SemiPar)

data(onions)
str(onions)
onions.data <- onions[, -3] # remove location

any(is.na(onions.data))
plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield",main = "Scatterplot of Onion Yield vs. Density")

# Helper Functions
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

m <- function(coefs, x) {coefs[1] + coefs[2]*x + 
    coefs[3]*x^2 + coefs[4]*x^3 + coefs[5]*x^4}
m2 <- function(coefs, x) {2 * coefs[3] + 6 * coefs[4]*x + 12 * coefs[5]*x^2}

bias.ci <- function(h, m2, mu2) {(1/2) * (h^2) * m2 * mu2}
sd.ci <- function(n, h, sigma2, nu0, f) {
  variance <- (1 / (n * h)) * (sigma2) * (nu0 / f)
  return(variance^(1/2))}

est_density <- function(data, xgrid, h, kernel = "gauss") {
  my.kernel <- switch(kernel, "gauss" = function(x) dnorm(x),
                      "epa"   = function(x) { 3/4 * (1 - x^2) *
                          ifelse(abs(x) <= 1, 1, 0) },
                      "uni"   = function(x) { 0.5 * ifelse(abs(x) <= 1, 1, 0) })
  n <- length(data)
  est <- numeric(length(xgrid))
  denom <- n * h
  for (j in seq_along(xgrid)) {
    est[j] <- sum(my.kernel((data - xgrid[j]) / h)) / denom}
  return(est)}

## Analytical Approach
ana.ll <- function(data, xgrid, h = 1, kernel = "gauss", alpha = 0.05) {
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
  # f <- pmax(f, 0.005)  # Prevents small values from inflating variance
  
  bias_val <- bias.ci(h, mddash, mu2)
  sd_val <- sd.ci(n, h, sigma2, nu0, f)
  
  est <- mhat - bias_val # Bias-corrected estimate
  z <- qnorm(1-alpha/2)
  lower_ci <- est - z * sd_val # Compute CI: (estimate - bias) - z * sd
  upper_ci <- est + z * sd_val # Compute CI: (estimate - bias) + z * sd
  
  ord <- order(xgrid) # Sort by xgrid to ensure lines() function plots properly
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  return(list(xgrid = xgrid_sorted, est = est_sorted,
              lower = lower_ci_sorted, upper = upper_ci_sorted))
}

## Simulated Approach
sim.ll <- function(data, xgrid, h = 1, B = 1000, alpha = 0.05) {
  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]
  mhat <- my.ll.smoother(xdat = x, ydat = y, xgrid = xgrid, h = h)
  mhat_at_x <- my.ll.smoother(xdat = x, ydat = y, xgrid = x, h = h)
  residuals <- y - mhat_at_x # to be resampled
  
  # Initialise bootstrap matrix
  boots <- matrix(NA, nrow = B, ncol = length(xgrid))
  for (b in 1:B) {
    resampled_resid <- sample(residuals, replace = TRUE)
    y_star <- mhat_at_x + resampled_resid # new prediction
    boots[b, ] <- my.ll.smoother(xdat = x, ydat = y_star, xgrid = xgrid, h = h)}
  
  lower_ci <- apply(boots, 2, quantile, probs = alpha / 2)
  upper_ci <- apply(boots, 2, quantile, probs = 1 - alpha / 2)
  ord <- order(xgrid)
  return(list(xgrid = xgrid[ord], est = mhat[ord],
    lower = lower_ci[ord], upper = upper_ci[ord]))
}

# Analysis
h_rot <- thumbBw(onions.data$dens, onions.data$yield, deg = 1, kernel = gaussK)
h_rot # 3.293061

thumbBw_manual <- function(x, y, kernel = "gauss") {
  n <- length(x)
  nu0 <- switch(kernel, "gauss" = 1 / (2 * sqrt(pi)), "epa" = 3 / 5, "uni" = 1 / 2)
  mu2 <- switch(kernel, "gauss" = 1, "epa"   = 1 / 5, "uni" = 1 / 3)
  
  poly_fit <- lm(y ~ poly(x, degree = 4, raw = TRUE))
  coefs <- coef(poly_fit)
  sigma2 <- sum((y - m(coefs, x))^2) / (n - 5)
  sum_m_double_sq <- sum(m2(coefs, x)^2)
  
  numerator <- nu0 * sigma2
  denominator <- mu2^2 * sum_m_double_sq
  h_rot <- (numerator / denominator)^(1/5) * n^(-1/5)
  return(h_rot)
}

h_rot_manual <- thumbBw_manual(onions.data$dens, onions.data$yield, kernel = "gauss")
h_rot_manual # 1.374277

summary(onions.data)
xgrid <- seq(15, 185, by = 0.5)

res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                              kernel = "gauss", alpha = 0.05)
res_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_rot_manual, 
                        alpha = 0.05, B = 300)
plot(onions.data$dens, onions.data$yield, # type = "n", 
     pch = 19,xlab = "Density", ylab = "Yield", main = "LL Regression with the 2 Approaches")

# Simulated Approach
lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)
# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "skyblue", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "gold", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "gold", lty = 2)
legend("topright", legend = c("Estimate", "Analytical CI", "Simulated CI"),
       col = c("skyblue", "gold", "purple"), lwd = 2, lty = c(1, 3, 3), bty = "n")

width_analytic <- mean(res_analytic$upper - res_analytic$lower)
width_bootstrap <- mean(res_bootstrap$upper - res_bootstrap$lower)
interval_widths <- data.frame(Method = c("Analytical Approach", "Simulated Approach"),
                              Average_Width = c(width_analytic, width_bootstrap))
interval_widths

# Diagnosing the analytical approach:
f <- est_density(onions.data$dens, xgrid, h_rot_manual, kernel="gauss")

plot(xgrid, f, type = "l", col = "darkgreen", lwd = 2,
     ylab = "Estimated Density", main = "Kernel Density Estimate at xgrid")


model <- lm(onions.data$yield ~ poly(onions.data$dens, 4, raw = TRUE))
summary(model)

hist(residuals(model), breaks = 30, col = "grey", main = "Residuals Histogram")

nu0 <- 1 / (2 * sqrt(pi)) # switch(kernel,"gauss" = 1 / (2 * sqrt(pi)),"epa" = 3 / 5,"uni" = 1 / 2)
mu2 <- 1 # switch(kernel,"gauss" = 1,"epa" = 1 / 5,"uni" = 1 / 3)
n <- nrow(onions.data)
x <- onions.data[, 1]
y <- onions.data[, 2]

model <- lm(y ~ poly(x, 4, raw = TRUE))   # Fit degree 4 polynomial for derivatives
coefs <- coef(model)
mhat <- my.ll.smoother(xdat = x, ydat = y, xgrid = xgrid, h = h_rot_manual)
mddash <- m2(coefs, xgrid)
sigma2 <- sum((y - m(coefs, x))^2) / (n - 5)
f <- est_density(x, xgrid, h_rot_manual, kernel="gauss")

bias_val <- bias.ci(h_rot_manual, mddash, mu2)
sd_val <- sd.ci(n, h_rot_manual, sigma2, nu0, f)

est <- mhat - bias_val # Bias-corrected estimate

plot(xgrid, mhat, type = "l", col = "blue", lwd = 2,
     ylab = "Estimate", main = "Bias-Corrected vs Uncorrected Estimate")
lines(xgrid, est, col = "black", lty = 2, lwd = 2)
legend("topright", legend = c("Uncorrected (mhat)", "Bias-corrected (est)"),
       col = c("blue", "black"), lty = c(1, 2), lwd = 2, bty = "n")
# Bias is not an issue!

# Appendix

## Build of the `thumbBw` Function
thumbBw_automatic <- function (x, y, deg, kernel, weig = rep(1, length(y))) 
{
  k <- 3
  rd <- compDerEst(x, y, deg, weig)
  denom <- sum(rd$der^2)
  numer <- mean(rd$res^2)
  cte <- cteNuK(0, deg, kernel, lower = dom(kernel)[[1]], upper = dom(kernel)[[2]], 
      subdivisions = 100)
  res <- cte * (numer/denom)^(1/(2 * deg + k))
  return(res)
}

## Adjusted Analytical Function
ana.adj.ll <- function(data, xgrid, h = 1, kernel = "gauss", alpha = 0.05) {
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
  f <- pmax(f, 0.005)  # Prevents small values from inflating LL variance
  
  bias_val <- bias.ci(h, mddash, mu2)
  sd_val <- sd.ci(n, h, sigma2, nu0, f)
  
  est <- mhat - bias_val # Bias-corrected estimate
  z <- qnorm(1-alpha/2)
  lower_ci <- est - z * sd_val # Compute CI: (estimate - bias) - z * sd
  upper_ci <- est + z * sd_val # Compute CI: (estimate - bias) + z * sd
  
  ord <- order(xgrid) # Sort by xgrid to ensure lines() function plots properly
  xgrid_sorted <- xgrid[ord]
  est_sorted <- est[ord]
  lower_ci_sorted <- lower_ci[ord]
  upper_ci_sorted <- upper_ci[ord]
  return(list(xgrid = xgrid_sorted, est = est_sorted,
              lower = lower_ci_sorted, upper = upper_ci_sorted))
}

## Epanechikov Kernel
my.kernel<-function(u){
  u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

h_rot_epa_manual <- thumbBw_manual(onions.data$dens, onions.data$yield, kernel = "epa")
h_rot_epa_manual

h_rot_epa <- thumbBw(onions.data$dens, onions.data$yield, deg = 1, kernel = EpaK)
h_rot_epa

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Epanechikov Analytic Approach")
res_analytic <- ana.adj.ll(onions.data, xgrid = xgrid, h = h_rot_epa_manual, 
                       kernel = "epa", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "skyblue", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "gold", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "gold", lty = 2)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Epanechikov Analytic Approach")
res_analytic <- ana.adj.ll(onions.data, xgrid = xgrid, h = h_rot_epa, 
                           kernel = "epa", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "skyblue", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "gold", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "gold", lty = 2)



## Uniform Kernel
my.kernel<-function(u){u<- 0.5*ifelse(abs(u)<=1, 1,0)}

h_rot_uni_manual <- thumbBw_manual(onions.data$dens, onions.data$yield, kernel = "uni")
h_rot_uni_manual
# there does not exist an option to input the uniform kernel into thumbBw

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "Uniform Analytic Approach")
res_analytic <- ana.adj.ll(onions.data, xgrid = xgrid, h = h_rot_uni_manual, 
                       kernel = "uni", alpha = 0.05)
lines(res_analytic$xgrid, res_analytic$est, col = "skyblue", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "gold", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "gold", lty = 2)

## Automatic Implementation of the "Rule-of-Thumb" Bandwidth
my.kernel<-function(u){return(dnorm(u))}

res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_rot, 
                              kernel = "gauss", alpha = 0.05)
res_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_rot, 
                        alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield", main = "LL Regression with Analytical and Simulated CIs")

# Simulated Approach
lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)

# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "skyblue", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "gold", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "gold", lty = 2)

legend("topright", legend = c("Estimate", "Analytical CI", "Simulated CI"),
       col = c("skyblue", "gold", "purple"), lwd = 2, lty = c(1, 3, 3), bty = "n")

width_analytic <- mean(res_analytic$upper - res_analytic$lower)
width_bootstrap <- mean(res_bootstrap$upper - res_bootstrap$lower)

interval_widths <- data.frame(Method = c("Analytical Approach", "Simulated Approach"),
                              Average_Width = c(width_analytic, width_bootstrap))
interval_widths

## Dpill Bandwidth Implementation
library(KernSmooth)
h_dpill <- dpill(onions.data$dens, onions.data$yield)

res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_dpill, 
                              kernel = "gauss", alpha = 0.05)
res_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_dpill, 
                        alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield",main = "LL Regression with Analytical and Simulated CIs")

# Simulated Approach
# lines(res_bootstrap$xgrid, res_bootstrap$est, col = "purple", lwd = 2)
lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)

# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "skyblue", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "gold", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "gold", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("gold", "purple"), lwd = 2, lty = c(2, 3), bty = "n")

width_analytic <- mean(res_analytic$upper - res_analytic$lower)
width_bootstrap <- mean(res_bootstrap$upper - res_bootstrap$lower)

interval_widths <- data.frame(Method = c("Analytical Approach", "Simulated Approach"),
                              Average_Width = c(width_analytic, width_bootstrap))
interval_widths

# Rmd version of the above:

## Dpill Bandwidth Implementation
\label{dpill}
This uses the package `KernSmooth`.
```{r dpill}
library(KernSmooth)
h_dpill <- dpill(onions.data$dens, onions.data$yield)
h_dpill
```

```{r dpill plot}
my.kernel<-function(u){return(dnorm(u))}

res_analytic <- ana.ll(onions.data, xgrid = xgrid, h = h_dpill, 
                       kernel = "gauss", alpha = 0.05)
res_bootstrap <- sim.ll(onions.data, xgrid = xgrid, h = h_dpill, 
                        alpha = 0.05, B = 300)

plot(onions.data$dens, onions.data$yield, pch = 19,xlab = "Density", 
     ylab = "Yield",main = "LL Regression with Analytical and Simulated CIs")

# Simulated Approach

lines(res_bootstrap$xgrid, res_bootstrap$lower, col = "purple", lty = 2)
lines(res_bootstrap$xgrid, res_bootstrap$upper, col = "purple", lty = 2)

# Analytical Approach
lines(res_analytic$xgrid, res_analytic$est, col = "skyblue", lwd = 2)
lines(res_analytic$xgrid, res_analytic$lower, col = "gold", lty = 2)
lines(res_analytic$xgrid, res_analytic$upper, col = "gold", lty = 2)

legend("topright", legend = c("Analytical CI", "Simulated CI"),
       col = c("gold", "purple"), lwd = 2, lty = c(2, 3), bty = "n")
```
