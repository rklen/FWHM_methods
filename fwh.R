#Methods for computing FWHM

#x is a vector whose ith element is a middlepoints of the ith interval
#y is a vector whose ith element is the number of observations on the ith interval
#some elements of y might be zero but they are dealt in the methods if necessary

#Methods F1-F5 from Markevich-Gertner

f1<-function(x,y){
  #let j be the index of the middle-most maximum
  j=which(y==max(y))
  if(length(j)>1){
    j=j[round(length(j)/2+0.5)]
  }
  #let l,r be the closest indexes with values of y less than half maximum on each side of j
  l<-j-1
  while(y[l]>=y[j]/2){
    l<-l-1
  }
  r<-j+1
  while(y[r]>=y[j]/2){
    r<-r+1
  }
  #use these indexes to estimate x at half maximum
  c_l=(x[l]+x[l+1])/2
  c_r=(x[r-1]+x[r])/2
  #FWHM is the difference of c_r-c_l
  return(c_r-c_l)
}

f2<-function(x,y){
  #let j be the index of the middle-most maximum
  j=which(y==max(y))
  if(length(j)>1){
    j=j[round(length(j)/2+0.5)]
  }
  #let d be the width of the jth interval
  d<-((x[j+1]+x[j])/2-(x[j-1]+x[j])/2)
  #let N be all the observations
  N<-sum(y)
  #count the y_j/d
  #return FWHM as sqrt(8ln(2)/(2pi))Nd/y_j
  sqrt(8*log(2)/(2*pi))*N*d/y[j]
}

f3<-function(x,y){
  #count the following sums
  T1<-T2<-0
  for(i in 1:length(y)){
    T1<-T1+x[i]*y[i]
    T2<-T2+x[i]^2*y[i]
  }
  #estimate standard deviation
  sigma3<-sqrt(T2/sum(y)-(T1/sum(y))^2)
  #return FWHM as sigma*sqrt(8ln(2))
  return(sigma3*sqrt(8*log(2)))
}

f4<-function(x,y){
  #create a new vector without counts <=3
  y1=y[y>3]
  #fit a parabola to the data a*log(y1)^2+b*log(y1)+c
  fit<-lm(log(y1)~poly(x[y>3],2,raw=T))
  #find the coefficient of the 2nd degree term
  a<-as.numeric(fit$coefficients[3])
  #estimate the standard deviation
  sigma4<-sqrt(-1/(2*a))
  #return FWHM as sigma*sqrt(8ln(2))
  return(sigma4*sqrt(8*log(2)))
}

f5<-function(x,y){
  #remove the counts <=3
  y1=y[y>3]
  x1=x[y>3]
  #create a new data ln(y_{i+1}/y_{i-1})/(x_{i+1}-x_{i-1}) for i=2,...,n-1
  l1<-c()
  for(i in 2:(length(y1)-1)){
    l1<-c(l1,log(y1[i+1]/y1[i-1])/(x1[i+1]-x1[i-1]))
  }
  #fit a line to this data
  fit<-lm(l1~poly(x1[2:(length(x1)-1)],1,raw=T))
  A<-as.numeric(fit$coefficients[2])
  #estimate the standard deviation
  sigma5<-sqrt(1/abs(-A))
  #return FWHM as sigma*sqrt(8ln(2))
  return(sigma5*sqrt(8*log(2)))
}

#method of using the intersection points found with linear interpolation
f6<-function(x,y){
  #let j be the index of the middle-most maximum
  j=which(y==max(y))
  if(length(j)>1){
    j=j[round(length(j)/2+0.5)]
  }
  #let l,r be the closest indexes with values of y less than half maximum on each side of j
  l<-j-1
  while(y[l]>=y[j]/2){
    l<-l-1
  }
  r<-j+1
  while(y[r]>=y[j]/2){
    r<-r+1
  }
  #compute the intersection points
  c_l=x[l]+(y[j]/2-y[l])/(y[l+1]-y[l])*(x[l+1]-x[l])
  c_r=x[r]-(y[j]/2-y[r])/(y[r-1]-y[r])*(x[r]-x[r-1])
  #FWHM is the difference of c_r-c_l
  return(c_r-c_l)
}

#method of using middle-most indexes for l,r
f7<-function(x,y){
  #let j be the index of the middle-most maximum
  j=which(y==max(y))
  if(length(j)>1){
    j=j[round(length(j)/2+0.5)]
  }
  #find all intervals containing the half maximum on each side of j
  L<-c()
  for(i in 1:(j-1)){
    if(y[i]<y[j]/2 & y[i+1]>=y[j]/2){
      L<-c(L,i)
    }
  }
  R<-c()
  for(i in j:length(y)){
    if(y[i]<y[j]/2 & y[i-1]>=y[j]/2){
      R<-c(R,i)
    }
  }
  #choose middle-most indexes for l,r
  if(length(L)==1){
    l=L
  }else{
    l<-L[round(length(L)/2+0.5)]
  }
  if(length(R)==1){
    r=R
  }else{
    r<-R[round(length(R)/2+0.5)]
  }
  #interpolate x at half maximum
  c_l=(x[l]+x[l+1])/2
  c_r=(x[r-1]+x[r])/2
  #FWHM is the difference of c_r-c_l
  return(c_r-c_l)
}

#method from NEMA
f8<-function(x,y){
  #let j be the index of the middle-most maximum
  j=which(y==max(y))
  if(length(j)>1){
    j=j[round(length(j)/2+0.5)]
  }
  #fit parabola to the points y_{j-1},y_j,y_{j+1}
  fit<-lm(y[(j-1):(j+1)]~poly(x[(j-1):(j+1)],2,raw=T))
  a<-as.numeric(fit$coefficients[3])
  b<-as.numeric(fit$coefficients[2])
  c<-as.numeric(fit$coefficients[1])
  #find the height of the parabola
  h<--b^2/(4*a)+c
  if(is.nan(h)){
    h<-max(y)
  }
  #let l,r be the closest indexes with values of y less than half of h
  l<-j-1
  while(y[l]>=h/2){
    l<-l-1
  }
  r<-j+1
  while(y[r]>=h/2){
    r<-r+1
  }
  #compute the intersection points
  c_l=x[l]+(h/2-y[l])/(y[l+1]-y[l])*(x[l+1]-x[l])
  c_r=x[r]-(h/2-y[r])/(y[r-1]-y[r])*(x[r]-x[r-1])
  #FWHM is the difference of c_r-c_l
  return(c_r-c_l)
}

# Function to be fitted (Gaussian distribution)
# A is numeric vector of parameters (A[1]=shift, A[2]=mean, A[3]=sd)
# x is numeric vector of function values
fit_function<-function(A,x){
  A[1]+1/(A[3]*sqrt(2*pi))*exp(-(1/2)*((x-A[2])/A[3])^2)
}
# Function for error (non-negative value)
# A is numeric vector of parameters
# x and y are numeric vectors with same length
error_function<-function(A,x,y){
  sum((fit_function(A,x)-y)^2)
}
#method of fitting Gaussian curve
f9<-function(x,y){
  #scale y so that the area under curve is 1
  y<-y/(sum(y)*mean(x[2:length(x)]-x[1:(length(x)-1)]))
  #count the following sums
  N1<-sum(y)
  T1<-T2<-0
  for(i in 1:length(y)){
    T1<-T1+x[i]*y[i]
    T2<-T2+x[i]^2*y[i]
  }
  #estimate mean and standard deviation
  mu9<-T1/N1
  sigma9<-sqrt(T2/N1-mu9^2)
  #let these be the initial values of parameters
  a_init=c(0,mu9,sigma9)
  #fit the curve
  tmp=optim(a_init,error_function,x=x,y=y)$par
  #use the estimate of sigma
  sigma9=as.numeric(tmp[3])
  #return FWHM as sigma*sqrt(8ln(2))
  return(sigma9*sqrt(8*log(2)))
}

n<-100
N<-100000
sigma<-1
s<-0
t<-0
h<-rep(0,9)
for(j in 1:1000){
  z<-rnorm(N,0,sigma)#+runif(N,-s/2,s/2)
  k<-seq(min(z),max(z),length.out=n+1)
  #k<-c(min(z),sort(runif(n-1,min(z),max(z))),max(z))
  y<-hist(z,breaks=k,plot=F)$counts#+rnorm(n,0,t)
  x<-c()
  for(i in 1:n){
    x<-c(x,(k[i]+k[i+1])/2)
  }
  h<-rbind(h,c(f1(x,y),f2(x,y),f3(x,y),f4(x,y),f5(x,y),f6(x,y),f7(x,y),f8(x,y),f9(x,y)))
  if(j%%100==0){
    print(j)
  }
}
h<-h[2:length(h[,1]),]
h<-h-sigma*sqrt(8*log(2))
write.csv(h,file='h001.csv',row.names=F)
h<-read.csv('h001.csv')
h1<-abs(h)/(sigma*sqrt(8*log(2)))
colMeans(h1)
y1<-colMeans(h1)
#y1<-rbind(y1,colMeans(h1))

write.csv(y1,file='y004.csv',row.names=F)

#h001 n=100,N=100000
#h002-h006 n=100,N=1000,3162,10000,31623,100000
#h007-h013 N=100000,n=10,32,100,316,1000,3162,10000
#h014-h024 n=100,N=1000,s=0,0.1,...,1
#h025-h035 n=100,N=1000,t=0,10,...,100

#Table 1
h<-read.csv('h001.csv')
i=9
mean(h[,i])
sd(h[,i])
mean(abs(h[,i]))
mean(abs(h[,i]))/sqrt(8*log(2))*100

#t-test: does the difference between estimated and real fwhm have zero mean?
h<-read.csv('h001.csv')
t.test(h[,9],mu=0)

#figure 1a
cols=c('black','lightblue','steelblue','blue','turquoise',
       'gray','purple','darkgreen','green')

y1<-read.csv('y001.csv')
plot(c(3,3.5,4,4.5,5),y1[,7],xlim=c(3,5),ylim=c(0,0.6),type='l',col=cols[7],ylab='',
     xlab='log_10(N)',cex.axis=1.5,cex.lab=2)
points(c(3,3.5,4,4.5,5),y1[,7],pch=8,col=cols[7])
for(i in c(1:6,8,9)){
  par(new=TRUE)
  plot(c(3,3.5,4,4.5,5),y1[,i],xlim=c(3,5),ylim=c(0,0.6),type='l',
       axes=F,xlab='',ylab='',col=cols[i])
  points(c(3,3.5,4,4.5,5),y1[,i],pch=8,col=cols[i])
}
legend('topright',legend=c('F1','F2','F3','F4','F5','F6','F7','F8','F9'),pch=8,col=cols,cex=1.6)

#figure 1b
y2<-read.csv('y002.csv')
plot(c(1,1.5,2,2.5,3,3.5,4),y2[,7],xlim=c(1,4),ylim=c(0,1),type='l',col=cols[7],ylab='',
     xlab='log_10(n)',cex.axis=1.5,cex.lab=2)
points(c(1,1.5,2,2.5,3,3.5,4),y2[,7],pch=8,col=cols[7])
for(i in c(1:6,8,9)){
  par(new=TRUE)
  plot(c(1,1.5,2,2.5,3,3.5,4),y2[,i],xlim=c(1,4),ylim=c(0,1),type='l',
       axes=F,xlab='',ylab='',col=cols[i])
  points(c(1,1.5,2,2.5,3,3.5,4),y2[,i],pch=8,col=cols[i])
}
#legend('topleft',legend=c('F1','F2','F3','F4','F5','F6','F7','F8','F9'),pch=8,col=cols,cex=1.5)

#figure 1c
y3<-read.csv('y003.csv')
plot(0.1*c(0:10),y3[,7],xlim=c(0,1),ylim=c(0,0.07),type='l',col=cols[7],ylab='',
     xlab='s',cex.axis=1.5,cex.lab=2)
points(0.1*c(0:10),y3[,7],pch=8,col=cols[7])
for(i in c(1:6,8,9)){
  par(new=TRUE)
  plot(0.1*c(0:10),y3[,i],xlim=c(0,1),ylim=c(0,0.07),type='l',
       axes=F,xlab='',ylab='',col=cols[i])
  points(0.1*c(0:10),y3[,i],pch=8,col=cols[i])
}
#legend('topleft',legend=c('F1','F2','F3','F4','F5','F6','F7','F8','F9'),pch=8,col=cols,cex=1.5)

#figure 1d
y4<-read.csv('y004.csv')
plot(10*c(0:10),y4[,7],xlim=c(0,100),ylim=c(0,0.55),type='l',col=cols[7],ylab='',
     xlab='t',cex.axis=1.5,cex.lab=2)
points(10*c(0:10),y4[,7],pch=8,col=cols[7])
for(i in c(1:6,8,9)){
  par(new=TRUE)
  plot(10*c(0:10),y4[,i],xlim=c(0,100),ylim=c(0,0.55),type='l',
       axes=F,xlab='',ylab='',col=cols[i])
  points(10*c(0:10),y4[,i],pch=8,col=cols[i])
}
#legend('topleft',legend=c('F1','F2','F3','F4','F5','F6','F7','F8','F9'),pch=8,col=cols,cex=1.5)

#t-test: is the relative error of f3 smaller than f9?
h<-read.csv('h001.csv')
h1<-abs(h)/(sigma*sqrt(8*log(2)))
mean(h1[,9]-h1[,3])
t.test(h1[,3],h1[,9],paired=T)

#t-test: does it affect to the error if breaks are uneven?
h<-read.csv('h001.csv')
h1<-abs(h)/(sigma*sqrt(8*log(2)))
colMeans(h1)
h<-read.csv('h023.csv')
h23<-abs(h)/(sigma*sqrt(8*log(2)))
colMeans(h23)
t.test(h1[,1],h23[,1],paired=F)
