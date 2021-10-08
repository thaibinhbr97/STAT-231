###################################################################################
# Run this code only once 
library(MASS)     # truehist is in the library MASS
###################################################################################

###################################################################################
# Problem 5: Linear regression model
id<-20687353
try(suppressWarnings(RNGkind(sample.kind="Rounding"), silent = TRUE))
set.seed(id)
x<-round(runif(100,1,20),digits=1)
alpha<-rnorm(1,0,5)
beta<- rnorm(1,0,5)
# display values of alpha and beta
cat("alpha = ", alpha, ", beta = ", beta)       
model<-sample(c(1:4),1)
cat("Model = ", model)
# Data are randomly generated depending on the value of the variable model
if (model==1) {
  y<-round(alpha+beta*x+rnorm(100,0,10),digits=1)
} else if (model==2) {
  y<-round(alpha+beta*x+rnorm(100,0,x),digits=1)
} else if (model==3) {
  y<-round(alpha+beta*((x-10)/5)^2+rnorm(100,0,3),digits=1)
} else if (model==4)  {
  y<-round(alpha+beta*x+3*rt(100,2),digits=1) }
# display sample correlation
cat("sample correlation = ", cor(x,y)) 
# run regression y = alpha+beta*x
RegModel<-lm(y~x)
summary(RegModel)     # parameter estimates and p-value for test of no relationship
alphahat<-RegModel$coefficients[1] # estimate of intercept
betahat<-RegModel$coefficients[2]  # estimate of slope
muhat<-RegModel$fitted.values      # fitted responses
r<- RegModel$residuals             # residuals
se<-summary(RegModel)$sigma        # estimate of sigma
rstar <- r/se                      # standardized residuals
# Scatterplot of data with fitted line
plot(x,y,col="blue")
title(main="Scatterplot with Fitted Line")
abline(a=alphahat,b=betahat,col="red",lwd=2)
# Residual plots
plot(x,rstar,xlab="x",ylab="Standardized Residual")
title(main="Residual vs x")
abline(0,0,col="red",lwd=1.5)
plot(muhat,rstar,xlab="Muhat",ylab="Standardized Residual")
title(main="Residual vs Muhat")
abline(0,0,col="red",lwd=1.5)
qqnorm(rstar,main="")
qqline(rstar,col="red",lwd=1.5)  # add line for comparison
title(main="Qqplot of Residuals")
# 95% Confidence interval for slope
confint(RegModel,level=0.95)
# 90% confidence interval for mean response at x=5
predict(RegModel,data.frame("x"=5),interval="confidence",lev=0.90)
# 99% Prediction interval for response at x=2
predict(RegModel,data.frame("x"=2),interval="prediction",lev=0.99)
# 95% confidence interval for sigma
df<-length(y)-2
a<-qchisq(0.025,df)
b<-qchisq(0.975,df)
cat("95% confidence interval for sigma: ",c(se*sqrt(df/b),se*sqrt(df/a)))
###################################################################################
print(se)