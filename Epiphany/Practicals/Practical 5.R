install.packages("carData") # you may need to uncomment this one, if using the package first time
require(carData)
data(Davis)

summary(Davis$weight)

attach(Davis)
weight

sum(is.na(weight))
hist(weight)

bopt<-function(vec){
  vec<-na.omit(vec)
  n<-length(vec)
  s<-sd(vec)
  b<- 2.42*s* n^(-1/3)
  return(b)
}
b_opt <- bopt(weight)

bscott<-function(vec){
  vec<-na.omit(vec)
  n<-length(vec)
  s<-sd(vec)
  b<- 3.5*s* n^(-1/3)
  return(b)
}
b_scott <- bscott(weight)

par(mfrow = c(1,1))
x0 <- 0
index <- 0:50
breaks= x0 + index*b_opt
manual <- hist(weight, breaks=breaks, freq=FALSE)

x0 <- 0
index <- 0:50
breaks= x0 + index*b_scott

scott_hist <- hist(weight, breaks="scott", freq=FALSE)
FD_hist <- hist(weight, breaks="FD", freq=FALSE)

manual$breaks
FD_hist$breaks
# scott_hist$breaks

