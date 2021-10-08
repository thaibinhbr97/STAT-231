###################################################################################

# Run this code only once. 

library(MASS)     # truehist is in the library MASS

#

###################################################################################

id<-20687353

###################################################################################

# Problem 2: Run this code which plots the Exponential relative likelihood function and the line y = 0.15

# on the same graph. The line can be used to determine a 15% likelihood interval.

set.seed(id)

theta<-rexp(1, 1/(max(1,id-10*trunc(id/10))))    #generate a random value of theta

if (theta<1) {theta<-theta+1}   # avoid small values of theta

n<-40

# generate random sample of n observation from Exponential(theta) distribution 

y<-sort(rexp(n,1/theta))    

y[1:3]                 # display first 3 numbers in the data set

thetahat<-mean(y)    #maximum likelihood estimate of theta

cat("n = ",n,", theta = ",theta,", thetahat = ",thetahat)   # display values

s<-thetahat/(n^0.5)

# interval of values for plotting relative likelihood function

th<-seq(max(0,thetahat-4*s), thetahat+4*s,0.001)

# create function to calculate Exponential relative likelihood function 

ExpRLF<-function(x) {(thetahat/x)^n*exp(n*(1-thetahat/x))}

plot(th,ExpRLF(th),xlab="theta",ylab="Rtheta",type="l",lwd=2)    # plot relative likelihood function

# draw a horizontal line at 0.15

abline(a=0.15,b=0,col="red",lwd=2)

title(main="Exponential Relative Likelihood Function")

###################################################################################



###################################################################################

# R code for finding the endpoints of a 15% likelihood interval for Exponential example

uniroot(function(x) ExpRLF(x)-0.15,lower=1,upper=1.5)

uniroot(function(x) ExpRLF(x)-0.15,lower=2,upper=2.5)

###################################################################################
