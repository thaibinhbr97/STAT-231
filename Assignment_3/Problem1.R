###################################################################################
# Run this code only once 
library(MASS)     # truehist is in the library MASS
###################################################################################

###################################################################################
# Problem 1: Binomial confidence intervals and sampling distribution of likelihood ratio statistic
id<-20687535
set.seed(id)
#   generate a random value of theta
theta<-rbeta(1, max(1,id-10*trunc(id/10)), max(1,trunc(id/10)-10*trunc(id/100))) 
if (theta<0.1) {theta<-theta+0.1}   # avoid small values of theta
if (theta>0.9) {theta<-theta-0.1}   # avoid large values of theta
n<-30
cat("n = ",n," theta = ",theta)   # display values
# vector of observations for 5000 simulations from a Binomial(n,theta) distribution
yobs<- rbinom(5000,n,theta)
# corresponding vector of thetahat values
that<-yobs/n
# values used to construct an approximate 95% confidence interval  based on Gaussian approximation                               
pm<-1.96*sqrt(that*(1-that)/n)       
# each approximate 95% confidence interval is stored in a row of matrix cibi
cibi<-matrix(c(that-pm,that+pm),nrow=5000,byrow=F)
cibi[1:10,1:2]         # Look at first 10 approximate 95% confidence intervals
# display proportion of approximate 95% confidence intervals which contain true value of theta
prop<- mean(abs(theta-that)<pm)
cat("proportion of approximate 95% confidence intervals which contain true value of theta = ",prop)
#
# create function to calculate Binomial relative likelihood function
BinRLF <- function(x) {dbinom(y,n,x)/dbinom(y,n,thetahat)}
li<-rep(0,2*5000)
li<- matrix(li,ncol=2,byrow=TRUE)                   # initialize matrix to store likelihood intervals
# For the 5000 simulations determine 15% likelihood intervals which are also 
# approximate 95% likelihood intervals 
for (i in 1:5000) {
  y<-yobs[i]
  thetahat<-that[i]
  if (thetahat==0) { li[i,1]<-0}        # if thetahat=0 then likelihood interval has left endpoint = 0
  else {result<-uniroot(function(x) BinRLF(x)-0.15,lower=0,upper=thetahat)
  li[i,1]<-result$root}
  if (thetahat==1) { li[i,2]<-1}       # if thetahat=1 then likelihood interval has right endpoint = 1
  else {result<-uniroot(function(x) BinRLF(x)-0.15,lower=thetahat,upper=1)
  li[i,2]<-result$root}
}
li[1:10,1:2]         # Look at first ten 15% likelihood intervals 
# display proportion of 15% likelihood intervals which contain the true value of theta
prop<- mean(theta>=li[,1] & theta<=li[,2])
cat("proportion of 15% likelihood intervals which contain true value of theta = ",prop)
#
# calculate the likelihood ratio statistic for all 5000 simulations and plot a relative histogram of values
# the histogram approximates the sampling distribution of the likelihood ratio statistic
lambda<-(-2*log(dbinom(yobs,n,theta)/dbinom(yobs,n,that)))
truehist(lambda,h=0.5,xlab="Likelihood Ratio Statistic",main="Sampling Distribution of Likelihood Ratio Statistic")
curve(dchisq(x,1), from=0.001,to=12,add=TRUE,col="red",lwd=2) # superimpose Chi-squared (1) pdf
#
# use number of trials = 100 for the Binomial experiment
n<-100
cat("n = ",n," theta = ",theta)   # display values
# vector of observations for 5000 simulations from a Binomial(n,theta) distribution
yobs<- rbinom(5000,n,theta)
# corresponding vector of thetahat values
that<-yobs/n
# values used to construct an approximate 95% confidence interval  based on Gaussian approximation                               
pm<-1.96*sqrt(that*(1-that)/n)       
# each approximate 95% confidence interval is stored in a row of matrix cibi
cibi<-matrix(c(that-pm,that+pm),nrow=5000,byrow=F)
cibi[1:10,1:2]         # Look at first 10 approximate 95% confidence intervals for theta
# display proportion of approximate 95% confidence intervals which contain true value of theta
prop<- mean(abs(theta-that)<pm)
cat("proportion of approximate 95% confidence intervals which contain true value of theta = ",prop)
#
# For the 5000 simulations determine 15% likelihood intervals which are also
# approximate 95% likelihood intervals 
for (i in 1:5000) {
  y<-yobs[i]
  thetahat<-that[i]
  if (thetahat==0) { li[i,1]<-0}       # if thetahat=0 then likelihood interval has left endpoint = 0
  else {result<-uniroot(function(x) BinRLF(x)-0.15,lower=0,upper=thetahat)
  li[i,1]<-result$root}
  if (thetahat==1) { li[i,2]<-1}     # if thetahat=1 then likelihood interval has right endpoint = 1
  else {result<-uniroot(function(x) BinRLF(x)-0.15,lower=thetahat,upper=1)
  li[i,2]<-result$root}
}
li[1:10,1:2]         # Look at first ten 15% likelihood intervals 
# display proportion of 15% likelihood intervals which contain true value of theta
prop<- mean(theta>=li[,1] & theta<=li[,2])
cat("proportion of 15% likelihood intervals which contain true value of theta = ",prop)
#
# calculate the likelihood ratio statistic for all 5000 simulations and plot a relative histogram of values
# the histogram approximates the sampling distribution of the likelihood ratio statistic
lambda<-(-2*log(dbinom(yobs,n,theta)/dbinom(yobs,n,that)))
truehist(lambda,h=0.5,xlab="Likelihood Ratio Statistic",main="Sampling Distribution of Likelihood Ratio Statistic")
curve(dchisq(x,1), from=0.001,to=12,add=TRUE,col="red",lwd=2) # superimpose Chi-squared (1) pdf
###################################################################################
#
###################################################################################