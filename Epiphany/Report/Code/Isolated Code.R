install.packages("SemiPar")
install.packages("KernSmooth")
install.packages("locpol")

require(KernSmooth)
require(locpol)
require(SemiPar)

my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

# Nadaraya-Watson (NW) estimator
my.kernel.smoother<- function(xdat, ydat, xgrid=xdat, h){
  G<- length(xgrid)
  est<-rep(0,G)
  for (j in 1:G){
    est[j]<-
      sum(ydat*my.kernel((xgrid[j]-xdat)/h))/sum(my.kernel((xgrid[j]-xdat)/h))
  }
  return(as.data.frame(cbind(xgrid, est)))
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

data(onions)

onions.data <- onions[, -3] # remove location variable (because it is binary)
head(onions.data)
colnames(onions.data)
dens.grid <- seq(15, 185, by=0.5)

# Initial Method:
fit1<- locpol(yield~dens, bw=2, kernel=gaussK, xeval=dens.grid, data=onions.data)
fit1$bw
fit1$deg
plot(fit1)

# Manual Approach to the above:
fit2<- my.ll.smoother(xdat=onions.data$dens, ydat=onions.data$yield, xgrid=dens.grid, h=2)
plot(fit2, main="h=2")
lines(fit2, col=2, lwd=3)
lines(fit1$lpFit[1:2], lty=2, lwd=2)
legend(75,225, c("our function", "locpol" ), lwd=c(3,2), lty=c(1,2), col=c(1,2))

# Now let’s look into bandwidth selection. For the “rule-of-thumb”, we get:
thumbBw(onions.data$dens, onions.data$yield, kernel=gaussK, deg=1) # 3.293061
fit3<- locpol(yield~dens, bw=3.293061, kernel=gaussK, xeval=dens.grid, data=onions.data)
plot(fit3, main="h=thumbBw")

# Extended Method
dpill(onions$dens, onions$yield) # 7.921906
fit<- locpol(yield~dens, bw=7.921906, kernel=gaussK, xeval=dens.grid, data=onions.data)
names(fit)
plot(fit, main="h=dpill")
fit$residuals

# Yes this is now a bit better. Note that the difference between the optimal bandwidths of the two approaches is smaller for the
# onions data, which are homoscedastic, and hence may more likely deliver a good-performing ROT bandwidth.

# Built-In method
locpol.bootstrap <- function(x, y, B = 200,
                                kernel = gaussK,
                                alpha = 0.05) {
  h <- dpill(x, y) # specify bandwidth using a non-parametric estimate
  xgrid <- seq(min(x), max(x), length.out = 100)
  
  # Fit the original model once first
  original_fit <- locpol(y ~ x, bw = h, kernel = kernel,
                         xeval = xgrid, data = data.frame(x, y))
  
  # Extract the fitted values at the grid
  est0 <- original_fit$lpFit[, "y"]  
  
  n <- length(x)
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid)) # initialise bootstrap matrix
  
  # Loop over bootstrap resamples
  for(b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_b <- x[idx]
    y_b <- y[idx]
    
    # Fit locpol on the bootstrap sample and Store the fitted curve
    fit_b <- locpol(y_b ~ x_b, bw = h, kernel = kernel,
                    xeval = xgrid, data = data.frame(x_b, y_b))
    bootMat[b, ] <- fitted(fit_b, deg = 0)
  }
  
  # Finally, compute pointwise bootstrap CIs
  lower <- apply(bootMat, 2, quantile, probs = alpha/2)
  upper <- apply(bootMat, 2, quantile, probs = 1 - alpha/2)
  
  list(xgrid = xgrid, est = est0, lower = lower, upper = upper, h = h)
}

res <- locpol.bootstrap(onions.data$dens, onions.data$yield, B=300)
# cat("Chosen bandwidth was:", res$h, "\n")

plot(onions.data$dens, onions.data$yield, pch=19,
     main="Built-In")
lines(res$xgrid, res$est, col="purple", lwd=2)
lines(res$xgrid, res$lower, col="purple", lty=2)
lines(res$xgrid, res$upper, col="purple", lty=2)
legend(169, 272, 
       legend = c("Fitted curve", "Bootstrap CI"),
       col = c("purple", "purple"), 
       lty = c(1, 2), 
       title = "Key:", 
       cex = 1.5)

# Manual Method
my.ll.bootstrap <- function(x, y, B = 200, h = NULL, alpha = 0.05) {
  # Determine the bandwidth
  h <- dpill(x,y)
  xgrid <- seq(min(x), max(x), length.out = 100)
  
  # Fit local linear smoother to the full data.
  original_fit <- my.ll.smoother(xdat = x, ydat = y, xgrid = xgrid, h = h)
  est0 <- original_fit$est
  
  # Initialise a matrix to store the fitted values from each bootstrap process.
  n <- length(x)
  bootMat <- matrix(NA, nrow = B, ncol = length(xgrid))
  
  # Carry out bootstrap loop
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]
    
    # Fit the manual local linear smoother
    fit_boot <- my.ll.smoother(xdat = x_boot, ydat = y_boot, xgrid = xgrid, h = h)
    bootMat[b, ] <- fit_boot$est
  }
  
  lower <- apply(bootMat, 2, quantile, probs = alpha / 2)
  upper <- apply(bootMat, 2, quantile, probs = 1 - alpha / 2)
  
  list(xgrid = xgrid, est = est0, lower = lower, upper = upper, h = h)
}

res <- my.ll.bootstrap(onions.data$dens, onions.data$yield, B = 300)
plot(onions.data$dens, onions.data$yield, pch = 19,
     main = "Manual")
lines(res$xgrid, res$est, col = "purple", lwd = 2)
lines(res$xgrid, res$lower, col = "purple", lty = 2)
lines(res$xgrid, res$upper, col = "purple", lty = 2)
legend(169, 272, 
       legend = c("Fitted curve", "Bootstrap CI"),
       col = c("purple", "purple"), 
       lty = c(1, 2), 
       title = "Key:", 
       cex = 1.5)
