###################################################################################
# Run this code only once 
library(MASS)     # truehist is in the library MASS
###################################################################################

###################################################################################
# Problem 3: Gaussian confidence intervals
# The following R code runs a simulation in which 95% confidence intervals for the mean mu and the      
# standard deviation sigma are  calculated for 5000 randomly generated Gaussian data sets
id<-20687353
set.seed(id)
mu<-id-10*trunc(id/10)                                    # mu = last digit of ID
sig<-max(1,trunc(id/10)-10*trunc(id/100))    # sig = second last digit of ID unless last digit is zero
cat("mu = ", mu, ", sigma = ", sig)       #display values of mu and sigma
yn<-rnorm(5000*25,mu,sig)               # generate G(mu,sig) observations
# each of the 5000 rows of the matrix yn contains 25 independent observations from a 
# G(mu,sig) distribution
yn<-matrix(yn,ncol=25,byrow=TRUE) 
ybar<-apply(yn,1,mean)                      # vector of 5000 means                       
s<-apply(yn,1,sd)                                  # vector of 5000 sample standard deviations
a<-qt(0.975,24)                # value from t tables for 95% confidence interval for mu 
pm<-a*s/sqrt(25)            # used to get 95% confidence interval for mu   
# each confidence interval for mu is stored in a row of matrix cimu
cimu<-matrix(c(ybar-pm,ybar+pm),nrow=5000,byrow=F)
cimu[1:10,1:2]       #Look at first ten 95% confidence intervals for mu
# proportion of 95% confidence intervals which contain the true value of mu
prop<- mean(abs(mu-ybar)<pm)
cat("proportion of 95% confidence intervals which contain true value of mu = ",prop)
# values from Chi-square distribution for 95% confidence interval for sigma
a<-qchisq(0.025,24)         
b<-qchisq(0.975,24)
# each confidence interval for sigma is stored in a row of matrix cisig
cisig<-matrix(c(sqrt(24*s^2/b), sqrt(24*s^2/a)),nrow=5000,byrow=F)
cisig[1:10,1:2]        #Look at first ten 95% confidence intervals for sigma
# proportion of 95% confidence intervals which contain the true value of sigma
prop<-mean(sig>=cisig[,1] & sig<=cisig[,2])
cat("proportion of 95% confidence intervals which contain true value of sigma = ",prop)
###################################################################################
#
###################################################################################