h.rot<- function(data){
  data<-na.omit(data)
  n<-length(data)
  sigma<-min(sd(data), IQR(data)/1.34)
  hopt<- 1.06*sigma*n^(-1/5)
  return(hopt)
}

length(onions.data$dens)
sd(onions.data$dens)
IQR(onions.data$dens)/1.34

h.rot(onions.data$dens)

summary(onions.data)

x <- onions.data$dens
y <- onions.data$yield
x_scaled <- scale(x)[, 1]  # mean 0, sd 1

h.rot(x_scaled)
