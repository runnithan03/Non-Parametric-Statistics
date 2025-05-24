install.packages("SemiPar")
install.packages("KernSmooth")
install.packages("locpol")

require(KernSmooth)
require(locpol)
require(SemiPar)

my.kernel<-function(u){
  ker <- switch(kernel,
                gauss = function(u) dnorm(u),
                epa   = function(u) {0.75 * (1 - u^2) * ifelse(abs(u) <= 1, 1, 0)},
                uni   = function(u) {0.5 * ifelse(abs(u) <= 1, 1, 0)},
                stop("Kernel not recognized. Use 'gauss', 'epa', or 'uni'.")
  )
  return(ker)
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

# Now the biases and variances
# bias <- (1/2)*(h^2)*m2*mu2
# variance <- (1/(n*h))*(sigma^2)*(nu0/f)
# m2(x), f(x) and sigma(x) to be decided

# Define the estimated function m(x)
m <- function(coefs, x) {
  coefs[1] +
    coefs[2]*x +
    coefs[3]*x^2 +
    coefs[4]*x^3 +
    coefs[5]*x^4
}

# Define the first derivative m'(x)
m_dash <- function(coefs, x) {
  coefs[2] +
    2 * coefs[3]*x +
    3 * coefs[4]*x^2 +
    4 * coefs[5]*x^3
}

# Define the second derivative m''(x)
m_double_dash <- function(coefs, x) {
  2 * coefs[3] +
    6 * coefs[4]*x +
    12 * coefs[5]*x^2
}

bias <- function(h, m2, mu2){
  return((1/2)*(h^2)*m2*mu2)
}

sd <- function(n, h, sigma, nu0, f){
  return(((1/(n*h))*(sigma^2)*(nu0/f))^(1/2))
}

est_density <- function(data, xgrid, h = 1, kernel = "gauss") {
  # Define the kernel function based on the chosen kernel
  my.kernel <- switch(kernel,
                      "gauss" = function(x) dnorm(x),
                      "epa"   = function(x) { 3/4 * (1 - x^2) * ifelse(abs(x) <= 1, 1, 0) },
                      "uni"   = function(x) { 0.5 * ifelse(abs(x) <= 1, 1, 0) }
  )
  
  # Remove any missing values
  data <- na.omit(data)
  
  n <- length(data)
  g <- length(xgrid)
  est <- numeric(g)
  denom <- n * h
  
  # Loop over each point in the evaluation grid and compute the density estimate
  for (j in seq_along(xgrid)) {
    num <- sum(my.kernel((data - xgrid[j]) / h))
    est[j] <- num / denom
  }
  
  # Return only the density estimates
  return(est)
}

data(onions)

onions.data <- onions[, -3] # remove location variable (because it is binary)
head(onions.data)
colnames(onions.data)
dens.grid <- seq(15, 185, by=0.5)

# Simulation Approach
sim.ll <- function(data, xgrid = dat, h = 1, kernel = "gauss", alpha = 0.05, B = 300) {
  my.kernel <- switch(kernel,
                      "gauss" = function(x) dnorm(x),
                      "epa"   = function(x) 3/4 * (1 - x^2) * ifelse(abs(x) <= 1, 1, 0),
                      "uni"   = function(x) 1/2 * ifelse(abs(x) <= 1, 1, 0))
  
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
  
  # Storage matrices (B rows, n cols)
  bootMat <- matrix(NA, nrow = B, ncol = n)
  bootUpperMat <- matrix(NA, nrow = B, ncol = n)
  bootLowerMat <- matrix(NA, nrow = B, ncol = n)
  xMat <- matrix(NA, nrow = B, ncol = n)
  
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    
    # Sort x_boot and y_boot before fitting the model
    sorted_idx <- order(x_boot)
    x_sorted <- x_boot[sorted_idx]
    y_sorted <- y_boot[sorted_idx]
    
    # Fit model on sorted x and y
    model <- lm(y_sorted ~ poly(x_sorted, 4, raw = TRUE))
    coefs <- coef(model)
    
    # Evaluate m(x), m''(x), density and variance all on sorted x
    mhat <- m(coefs, x_sorted)
    mddash <- m_double_dash(coefs, x_sorted)
    sigma2 <- sum((y_sorted - mhat)^2) / (n - 5)
    f <- est.density(x_sorted, x_sorted, h = h, kernel = kernel)
    
    bias_val <- bias(h, mddash, mu2)
    sd_val <- sd(n, h, sqrt(sigma2), nu0, f)
    zval <- qnorm(1 - alpha / 2)
    
    est <- mhat - bias_val
    upper <- est + zval * sd_val
    lower <- est - zval * sd_val
    
    # Store everything in sorted order
    bootMat[b, ]       <- est
    bootUpperMat[b, ]  <- upper
    bootLowerMat[b, ]  <- lower
    xMat[b, ] <- x_sorted
  }
  
  return(list(
    est   = bootMat,
    lower = bootLowerMat,
    upper = bootUpperMat,
    x_est = xMat
  ))
}



# Direct plug-in methodology to select the bandwidth...
h_dpill <- dpill(onions.data$dens, onions.data$yield) # Gaussian Kernel Regression Estimate

# Rule of Thumb Bandwidth - use this because its formula is in the notes!
h_rot <- thumbBw(onions.data$dens, onions.data$yield, deg=0, kernel=EpaK)

# sigma<-min(sd(data), IQR(data)/1.34) - from notes?

res <- sim.ll(onions.data, xgrid=seq(15, 185, by=0.5), h=h_rot, kernel="gauss", alpha=0.05, B=300)
plot(onions.data$dens, onions.data$yield, pch = 19,
     main = "Locpol with Bootstrap CIs")
lines(res$x_est, res$est, col = "purple", lwd = 2)
lines(res$x_est, res$lower, col = "purple", lty = 2)
lines(res$x_est, res$upper, col = "purple", lty = 2)

dim(res$est)
dim(res$x_est)

sim.test.ll <- function(data, xgrid = dat, h = 1, kernel = "gauss", alpha = 0.05, B = 300) {
  my.kernel <- switch(kernel,
                      "gauss" = function(x) dnorm(x),
                      "epa"   = function(x) 3/4 * (1 - x^2) * ifelse(abs(x) <= 1, 1, 0),
                      "uni"   = function(x) 1/2 * ifelse(abs(x) <= 1, 1, 0))
  
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
  
  # Storage matrices
  bootMat <- matrix(NA, nrow = B, ncol = n)
  xMat <- matrix(NA, nrow = B, ncol = n)
  
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    
    sorted_idx <- order(x_boot)
    x_sorted <- x_boot[sorted_idx]
    y_sorted <- y_boot[sorted_idx]
    
    model <- lm(y_sorted ~ poly(x_sorted, 4, raw = TRUE))
    coefs <- coef(model)
    
    mhat <- m(coefs, x_sorted)
    mddash <- m_double_dash(coefs, x_sorted)
    sigma2 <- sum((y_sorted - mhat)^2) / (n - 5)
    
    f <- est.density(x_sorted, h = h, kernel = kernel)
    
    bias_val <- bias(h, mddash, mu2)
    sd_val <- sapply(f, function(f_i) sd(n, h, sqrt(sigma2), nu0, f_i))
    
    zval <- qnorm(1 - alpha / 2)
    
    est <- mhat - bias_val
    bootMat[b, ] <- est
    xMat[b, ] <- x_sorted
  }
  
  # Sort x and y estimates in each row for consistent plotting
  for (b in 1:B) {
    sort_idx <- order(xMat[b, ])
    xMat[b, ] <- xMat[b, sort_idx]
    bootMat[b, ] <- bootMat[b, sort_idx]
  }
  
  # Final CI results
  xgrid_final <- apply(xMat, 2, mean)
  est_mean <- apply(bootMat, 2, mean)
  lower_ci <- apply(bootMat, 2, quantile, probs = alpha / 2)
  upper_ci <- apply(bootMat, 2, quantile, probs = 1 - alpha / 2)
  
  return(list(
    xgrid = xgrid_final,
    est   = est_mean,
    lower = lower_ci,
    upper = upper_ci
  ))
}

res <- sim.test.ll(onions.data, xgrid=seq(15, 185, by=0.5), h=h_rot, kernel="gauss", alpha=0.05, B=300)
plot(onions.data$dens, onions.data$yield, pch = 19,
     main = "Locpol with Bootstrap CIs",
     xlab = "Density", ylab = "Yield")

lines(res$xgrid, res$est, col = "purple", lwd = 2)
lines(res$xgrid, res$lower, col = "purple", lty = 2)
lines(res$xgrid, res$upper, col = "purple", lty = 2)

