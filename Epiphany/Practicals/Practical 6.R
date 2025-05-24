data(LakeHuron)

?hist
hist(LakeHuron, freq = FALSE) # default value of NULL means that no shading lines are drawn


my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

est.density<- function(data, xgrid=data, h=1, kernel="gauss"){
  my.kernel<- switch(kernel, 
                     "gauss"= function(x){dnorm(x)},
                     "epa"=function(x){x<- 3/4*(1-x^2)*ifelse(abs(x)<=1,1,0)},
                     "uni"=function(x){x<- 1/2*ifelse(abs(x)<=1,1,0)}   
  )
  data<-na.omit(data)
  n<-length(data)
  g<-length(xgrid)
  est<-rep(0, g)
  denom<- n*h
  for (j in 1:g){
    num<- sum(my.kernel((data-xgrid[j])/h))
    est[j]<- num/denom
  }  
  return(est)
}

x<- seq(575,582, by=0.1)
hist(LakeHuron, freq=FALSE, ylim=c(0,0.4), main="KDEs with LakeHuron",  col="grey95")
lines(x,est.density(LakeHuron, x, h=1, "gauss"), col=1, lwd=2)
# lines(x,est.density(LakeHuron, x, h=1, "gauss"), col=2, lwd=2)
# lines(x,est.density(LakeHuron, x, h=2, "gauss"), col=3, lwd=2)
# lines(x,est.density(LakeHuron, x, h=3, "gauss"), col=4, lwd=2)
# legend(-3, 0.1, c("h=0.5","h=1", "h=2", "h=3"), lty=c(1,1,1,1), lwd=c(2,2,2,2), col=c(1,2,3,4))

h.sil<- function(data){
  data<-na.omit(data)
  n<-length(data)
  sigma<-min(sd(data), IQR(data)/1.34)
  hopt<- 0.9*sigma*n^(-1/5)
  return(hopt)
}

h.sil(LakeHuron) # 0.4671343

dens <- density(LakeHuron)
dens$bw

my.kernel<-function(u){
  #u<- 0.5*ifelse(abs(u)<=1, 1,0)
  u<- dnorm(u)
  #u<- 3/4*(1-u^2)*ifelse(abs(u)<=1,1,0)
}

est.density.new<- function(data, xgrid=data, h=1, kernel="gauss"){
  my.kernel<- switch(kernel, 
                     "gauss"= function(x){dnorm(x)},
                     "epa"=function(x){x<- 3/4*(1-x^2)*ifelse(abs(x)<=1,1,0)},
                     "uni"=function(x){x<- 1/2*ifelse(abs(x)<=1,1,0)}   
  )
  data<-na.omit(data)
  n<-length(data)
  g<-length(xgrid)
  est<-rep(0, g)
  denom<- n*h
  for (j in 1:g){
    num<- sum(my.kernel((data-xgrid[j])/h))
    est[j]<- num/denom
  }  
  return(as.data.frame(cbind(xgrid, est)))
}

hist(LakeHuron, freq=FALSE, ylim=c(0,0.4), main="KDEs with LakeHuron",  col="grey95")
lines(est.density.new(LakeHuron, xgrid = x), col=1, lwd=2)


kde.reject<- function(kde, n){
  
  x.min <- min(kde$xgrid)
  x.max <- max(kde$xgrid)
  y.max <- max(kde$est)        # Highest KDE value
  
  samples<-rep(0, n)    # empty vector into which to place the drawn samples
  count   <- 0            # Counter for accepted samples
  
  while (count < n) {
    
    x.sim <- runif(1,x.min, x.max)
    
    y.sim <- runif(1,0, y.max)
    
    y.dens <- approx(kde$xgrid, kde$est, xout = x.sim)$y
    
    # Accept the sample if it falls below the KDE curve
    if (y.sim < y.dens) {
      count <- count + 1
      samples[count] <- x.sim
    }
  }
  
  return(samples)
}

?runif
lake.dens <- est.density.new(LakeHuron, xgrid = x, h = h.sil(LakeHuron), "gauss")
sim1 <- kde.reject(lake.dens, 100)
sim1

sort(sim1)

hist(sim1, freq=FALSE, ylim=c(0,0.4))
lines(est.density.new(LakeHuron, xgrid = sort(sim1), h = h.sil(LakeHuron)), col=1, lwd=2)



