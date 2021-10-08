###################################################################################
# Run this code only once 
library(MASS)     # truehist is in the library MASS
###################################################################################

###################################################################################
# Problem 2: Exponential confidence intervals and sampling distribution of likelihood ratio statistic
id<-20687351
set.seed(id)
theta<-max(1,id-10*trunc(id/10))                   # theta = last digit of ID unless it is zero
if (theta<1) {theta<-theta+1}   # avoid small values of theta
n<-5
cat("n = ",n," theta = ",theta)   # display values
ye<-rexp(5000*n,1/theta)    
# each of the 5000 rows of the matrix ye contains n independent observations from 
# Exponential(theta) distribution
ye<-matrix(ye,ncol=n,byrow=TRUE) 
that<-apply(ye,1,mean)                      # vector of 5000 means
pm<-1.96*that/sqrt(n)                     # used to get approximate 95% confidence interval 
# each approximate 95% confidence interval is stored in a row of matrix ciexp
ciexp<-matrix(c(that-pm,that+pm),nrow=5000,byrow=F)
ciexp[1:10,1:2]          # Look at first 10 approximate 95% confidence intervals 
# display proportion of approximate 95% confidence intervals which contain true value of theta
prop<- mean(abs(theta-that)<pm)
cat("proportion of approximate 95% confidence intervals which contain true value of theta = ",prop)
#
# create function to calculate Exponential relative likelihood function
ExpRLF<-function(x) {(thetahat/x)^n*exp(n*(1-thetahat/x))}
li<-rep(0,2*5000)
li<- matrix(li,ncol=2,byrow=TRUE)                   # initialize matrix to store likelihood intervals
# For the 5000 simulations determine 15% likelihood intervals which are also
# approximate 95% likelihood intervals 
for (i in 1:5000) {
  thetahat<-that[i]
  result<-uniroot(function(x) ExpRLF(x)-0.15,lower=max(0,thetahat-4*thetahat/(n^0.5)),upper=thetahat)
  li[i,1]<-result$root
  result<-uniroot(function(x) ExpRLF(x)-0.15,lower=thetahat,upper= thetahat+4*thetahat/(n^0.5))
  li[i,2]<-result$root
}
li[1:10,1:2]         # Look at first ten 15% likelihood intervals 
# display proportion of 15% likelihood intervals which contain the value of theta
prop<- mean(theta>=li[,1] & theta<=li[,2])
cat("proportion of 15% likelihood intervals which contain true value of theta = ",prop)
#
# calculate the likelihood ratio statistic for all 5000 simulations and plot a relative histogram of values
# the histogram approximates the sampling distribution of the likelihood ratio statistic
lambda<- -2*log((that/theta)^n*exp(n*(1-that/theta)))
truehist(lambda,h=0.5,xlab="Likelihood Ratio Statistic",main="Sampling Distribution of Likelihood Ratio Statistic")
curve(dchisq(x,1), from=0.001,to=12,add=TRUE,col="red",lwd=2) # superimpose Chi-squared (1) pdf
###################################################################################
#
###################################################################################