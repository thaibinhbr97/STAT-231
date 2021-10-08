###################################################################################
# Run this code only once 
library(MASS)     # truehist is in the library MASS
###################################################################################

###################################################################################
# Problem 4: One sample Gaussian model
id<-20687353
try(suppressWarnings(RNGkind(sample.kind="Rounding"), silent = TRUE))
set.seed(id) 
model<-sample(c(1:3),1)
cat("Model = ", model)
# Data are randomly generated from Gaussian distribution (model=1), Gamma distribution (model=2), or 
# Poisson distribution (model=3)
if (model==1) {
  mu<-id-10*trunc(id/10)                                     # mu = last digit of ID
  sigma<-max(1,trunc(id/10)-10*trunc(id/100))    # sig = second last digit of ID unless last digit is zero
  cat("mu = ", mu, ", sigma = ", sigma)       # display values of mu and sigma
  y<-sort(round(rnorm(30,mu,sigma),digits=2))  # 30 observations from G(mu,sig)
} else if (model==2) {
  mu<-max(1,id-10*trunc(id/10))                     # mu = last digit of ID unless it is zero
  y<-sort(round(rgamma(30,3,1/mu),digits=2))   # 30 observations from Gamma(3,1/mu)
  sigma<-(3*mu^2)^0.5
  cat("mu = ", mu, ", sigma = ", sigma)       # display values of mu and sigma
} else if (model==3) {
  mu<-max(1,id-10*trunc(id/10))                     # mu = last digit of ID unless it is zero
  y<- sort(round(rpois(30,mu),digits=2))   # 30 observations from Poisson(mu)
  sigma<-mu^0.5
  cat("mu = ", mu, ", sigma = ", sigma)       # display values of mu and sigma
} 
# Check Gaussian assumption using qqplot
qqnorm(y,xlab="Standard Normal Quantiles",main="Qqplot of Data")
qqline(y,col="red",lwd=1.5)  # add line for comparison
mu0<-mu+1
cat("mu0 = ", mu0)       # display value of mu0
# test hypothesis mu=mu0 and obtain 95% confidence interval for mu
t.test(y,mu=mu0,conf.level=0.95) 
#
sigma0<-sigma+2  
cat("sigma0 = ", sigma0)       # display value of sigma0
# test hypothesis sigma=sigma0 and obtain 95% confidence interval for sigma
df<-length(y)-1      # degrees of freedom
s2<-var(y)               
cat("sample variance = ", s2)          # display sample variance
chitest<-s2*df/sigma0^2
q<-pchisq(chitest,df)
cat("p-value for testing sigma=sigma0: ", min(2*q,2*(1-q)))       
p<-0.95                   # p=0.95 for 95% confidence interval
a<-qchisq((1-p)/2,df) # lower value from Chi-squared dist'n
b<-qchisq((1+p)/2,df) # upper value from Chi-squared dist'n
cat("95% confidence interval for sigma squared:  ",c(s2*df/b,s2*df/a))          
cat("95% confidence interval for sigma: ",c(sqrt(s2*df/b),sqrt(s2*df/a)))
###################################################################################
#
###################################################################################
