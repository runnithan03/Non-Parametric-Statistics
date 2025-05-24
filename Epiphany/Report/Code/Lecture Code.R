# Going through 4.2 using his lecture code...

install.packages("SemiPar")
require(SemiPar)

data(fossil)
plot(fossil)

colnames(fossil)
hist(fossil$age)
hist(fossil$strontium.ratio)

data(lidar)
plot(lidar)

colnames(lidar)     
hist(lidar$range)
hist(lidar$logratio)

# 4.1 Introduction
# 4.1.2 Motivating kernel regression via KDE

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

# We apply this function for a range of bandwidths on the fossil data….
age.grid<- seq(90, 130, by=0.5)
fit1<- my.kernel.smoother(fossil$age, fossil$strontium.ratio, age.grid, h=1)
fit2<- my.kernel.smoother(fossil$age, fossil$strontium.ratio, age.grid, h=2)
fit3<- my.kernel.smoother(fossil$age, fossil$strontium.ratio, age.grid, h=3)
plot(fossil)
lines(fit1, col=1, lwd=2)
lines(fit2, col=2, lwd=2)
lines(fit3, col=3, lwd=2)
legend(95,0.70731, c("h=1", "h=2", "h=3"), lwd=c(2,2,2), lty=c(1,1,1), col=c(1,2,3))

range.grid<-350:750
fitl1<- my.kernel.smoother(lidar$range, lidar$logratio, range.grid, h=10)
fitl2<- my.kernel.smoother(lidar$range, lidar$logratio, range.grid, h=20)
fitl3<- my.kernel.smoother(lidar$range, lidar$logratio, range.grid, h=30)
plot(lidar)
lines(fitl1, col=1, lwd=2)
lines(fitl2, col=2, lwd=2)
lines(fitl3, col=3, lwd=2)
legend(400,-0.4, c("h=10", "h=20", "h=30"), lwd=c(2,2,2), lty=c(1,1,1), col=c(1,2,3))

# 4.2 Localized Regression

# 4.2.2 Local Linear Regression Estimator 
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

# We apply this function now again on the fossil data.
age.grid<- seq(90, 130, by=0.5)
fit1<- my.ll.smoother(fossil$age, fossil$strontium.ratio, age.grid, h=1)
fit2<- my.ll.smoother(fossil$age, fossil$strontium.ratio, age.grid, h=2)
fit3<- my.ll.smoother(fossil$age, fossil$strontium.ratio, age.grid, h=3)
plot(fossil)
lines(fit1, col=1, lwd=2)
lines(fit2, col=2, lwd=2)
lines(fit3, col=3, lwd=2)
legend(95,0.70730, c("h=1", "h=2", "h=3"), lwd=c(2,2,2), lty=c(1,1,1), col=c(1,2,3))

# We now apply the local linear estimator on the lidar data.
range.grid<-350:750
fitl1<- my.ll.smoother(lidar$range, lidar$logratio, range.grid, h=10)
fitl2<- my.ll.smoother(lidar$range, lidar$logratio, range.grid, h=20)
fitl3<- my.ll.smoother(lidar$range, lidar$logratio, range.grid, h=30)
plot(lidar)
lines(fitl1, col=1, lwd=2)
lines(fitl2, col=2, lwd=2)
lines(fitl3, col=3, lwd=2)
legend(400,-0.4, c("h=10", "h=20", "h=30"), lwd=c(2,2,2), lty=c(1,1,1), col=c(1,2,3))

# 4.2.6 Bandwidth selection

# We initially use R function locpol to reproduce our earlier, manually implemented, estimates. 
# Note that the default setting of locpol is a polynomial degree deg=1 , that is local linear, 
# and that the default kernel is Epanechnikov, which we change now to Gaussian. 
# Let’s begin with the fossil data.

install.packages("locpol")
require(locpol)

fit2a<- locpol(strontium.ratio~age, bw=2, deg=1, kernel=gaussK, xeval=age.grid, data=fossil)
fit2a

plot(fit2a)

#Let’s inspect the fitted object more closely
names(fit2a)

fit2a$bw # h = bandwidth

# Does the result correspond to our manual implementation?
age.grid<- seq(90, 130, by=0.5)
fit2<- my.ll.smoother(fossil$age, fossil$strontium.ratio, age.grid, h=2)
plot(fossil, main="h=2")
lines(fit2, col=2, lwd=3)
lines(fit2a$lpFit[1:2], lty=2, lwd=2)
legend(95,0.70730, c("our function", "locpol" ), lwd=c(3,2), lty=c(1,2), col=c(1,2))

# Yes, there is perfect agreement between the two implementations.
# We can also do this for the lidar data:
fitl2a<- locpol(logratio~range, bw=20, kernel=gaussK, data=lidar)

# with associated plots
plot(fitl2a)

# Now let’s look into bandwidth selection. For the “rule-of-thumb”, we get:
thumbBw(fossil$age, fossil$strontium.ratio, kernel=gaussK, deg=1) # 0.6231473
thumbBw(lidar$range, lidar$logratio, kernel=gaussK, deg=1) # 4.178616

# with resulting estimates:
fit2b<- locpol(strontium.ratio~age, bw=0.6232473, kernel=gaussK, xeval=age.grid, data=fossil)
plot(fit2b)

fit2lb<- locpol(logratio~range, bw=4.1786, kernel=gaussK, xeval=range.grid, data=lidar)
plot(fit2lb)

# In both cases this seems to have settled on bandwidths which are rather too small. 
# I would suspect that an underestimation of has played a role in both. 
# Perhaps a full plug-in technique performs better? We find function dpill in R package
# KernSmooth ( “Least squares quartic fits over blocks of data are used to obtain an initial estimate.”):

install.packages("KernSmooth")
require(KernSmooth)

dpill(fossil$age, fossil$strontium.ratio) # 0.7875147
dpill(lidar$range, lidar$logratio) # 16.20422

fit2b<- locpol(strontium.ratio~age, bw=0.7875, kernel=gaussK, xeval=age.grid, data=fossil)
plot(fit2b)

fit2lb<- locpol(logratio~range, bw=16.204, kernel=gaussK, xeval=range.grid, data=lidar)
plot(fit2lb)

# Yes this is now a bit better. Note that the difference between the optimal bandwidths of the two approaches is larger for the
# lidar data, which are heteroscedastic, and hence may more likely deliver a poorly performing ROT bandwidth.


