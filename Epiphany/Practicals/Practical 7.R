## Gleanings from Lecture

# Example 3.14
require(plot3D)
require(rgl)
require(carData)
require(ks)

soils3<- Soils[,c("Mg", "Ca", "K")]
plot(soils3)

hd.rot<- function(data){
  data<-na.omit(data)
  n<-dim(data)[1]
  d<-dim(data)[2]
  hopt<-rep(0,d)
  for (j in 1:d){
    sigmaj<-sd(data[,j])
    hopt[j]<- sigmaj*n^{-1/(d+4)}
  }  
  return(hopt)    
}

h.soil3 <- hd.rot(soils3)
fit.soil3<-kde(soils3, H=diag(h.soil3^2))
plot(fit.soil3)

plot(fit.soil3, display="rgl")

soils4<- Soils[,c("Mg", "Ca", "K", "Na")]
plot(soils4)

h.soil4 <- hd.rot(soils4) 
fit.soil4<-kde(soils4, H=diag(h.soil4^2))
# plot(fit.soil4)     # commented out as it gives an error
names(fit.soil4)

# fit.soil4$estimate  # commented out as it needs lots of space

soils5<- Soils[,c("Mg", "Ca", "K", "Na", "P")]
plot(soils5)

h.soil5 <- hd.rot(soils5)
fit.soil5<-kde(soils5, H=diag(h.soil5^2))
soils6<- Soils[,c("Mg", "Ca", "K", "Na", "P", "N")]
plot(soils6)

h.soil6 <- hd.rot(soils6)
# fit.soil6<-kde(soils6, H=diag(h.soil6^2))
# this one took too long and in fact crashed my R session. Feel free to try :-)

data(Soils)
soils<- Soils[,c("Mg", "Ca")]
plot(soils)

h.soil2<-hd.rot(soils)
h.soil2
## [1] 0.717701 1.707931
h.soil3
## [1] 0.7869971 1.8728366 0.1288146
h.soil4
## [1] 0.8433256 2.0068830 0.1380344 2.0272460
h.soil5
## [1]  0.8899095  2.1177399  0.1456592  2.1392278 52.7886334
# These are gently increasing as expected. Let’s inspect a bit further.
# Which grid sizes (G) have been used for each coordinate axis? (Note that, the number of evaluation points of the KDE is Gd!)

length(fit.soil2$eval.points[[1]])
## [1] 151
length(fit.soil3$eval.points[[1]])
## [1] 51
length(fit.soil4$eval.points[[2]])
## [1] 21
length(fit.soil5$eval.points[[3]])
## [1] 21
# How many estimated density values does this result in?
  
length(fit.soil2$estimate)
## [1] 22801
length(fit.soil3$estimate)
## [1] 132651
length(fit.soil4$estimate)
## [1] 194481
length(fit.soil5$estimate)
## [1] 4084101
# Observe, for instance, that 214=194481.

 #How many (precisely) zero estimates are there?
  
sum(fit.soil2$estimate==0)
## [1] 3269
sum(fit.soil3$estimate==0)
## [1] 45473
sum(fit.soil4$estimate==0)
## [1] 0
sum(fit.soil5$estimate==0)
## [1] 0

# How many de-facto zero estimates are there? (Let’s call a KDE value f^H(x) “de facto zero” if f^H(x)<10−6×maxzf^H(z)).
fit.soil2 <- kde(soils, H=diag(h.soil2^2))
sqrt(fit.soil2$H)
sum(fit.soil2$estimate<1e-6*max(fit.soil2$estimate))
## [1] 3681
sum(fit.soil3$estimate<1e-6*max(fit.soil3$estimate))
## [1] 57419
sum(fit.soil4$estimate<1e-6*max(fit.soil4$estimate))
## [1] 124389
sum(fit.soil5$estimate<1e-6*max(fit.soil5$estimate))
## [1] 3346029

# Which proportion of estimates is precisely zero?
sum(fit.soil2$estimate==0)/length(fit.soil2$estimate)
## [1] 0.1433709
sum(fit.soil3$estimate==0)/length(fit.soil3$estimate)
## [1] 0.3428018
sum(fit.soil4$estimate==0)/length(fit.soil4$estimate)
## [1] 0
sum(fit.soil5$estimate==0)/length(fit.soil5$estimate)

