 library(PLCE)
#devtools::install_github("ratkovic/PLCE")
library(profvis)
library(microbenchmark)
library(ranger)
# devtools::load_all('~/Dropbox/Github/PLCE')


options(device="quartz")
rm(list=ls())
ran.num<-round(runif(1)*1e8)

## Trying random forest

theta.run<-se.run<-cover.run<-NULL
n<-2000
k<-5

meanhet<-errhet<-inter<-T
inter <-F
tol<-1e-3
verbose <- EM <- T

##Format data
var.mat   <- diag(k)
var.mat[var.mat==0] <- 0.5
res.true <- 0
# covariates


X <- MASS::mvrnorm(n, rep(0, k), Sig = var.mat)
X<-apply(X,2,scale)
X[,1]<-X[,1]/(mean(X[,1]^2)^.5)
X<-apply(X,2,FUN=function(x) x/(mean(x^2))^.5)
if(errhet)	sd.use<-((1+X[,1]+(length(unique(X[,1]))==2) )^2/2)^.5 else sd.use<-1
errst<-rnorm(n,sd=sd.use)#rt(n,8)*sd.use#
treat<-(X[,1])+errst +res.true
errsy<-rnorm(n,sd=1)#rt(n,8)*sd.use#
if(meanhet) Y<-1+treat*(X[,1]^2)+errsy+res.true else Y<-1+treat+(X[,1]^2)+errsy+res.true
##Shift X
if(inter){
  m1<-exp(-10*(outer(X[,1],X[,1],"-")^2))
  diag(m1)<-0
  oY<-m1%*%X[,1]/rowSums(m1)
  oT<-m1%*%X[,1]/rowSums(m1)
  
  Y<-Y+oY
  treat<-treat-oT
}

X2.model<-X
X2.model[,1]<-X[,1]+.5*X[,2]
X2.model[,2]<-X[,2]+.5*X[,1]
X <- X2.model
## Rotates X's

## GRF
g1 <-grf::causal_forest(Y=Y, treat,X=X)
grf::average_treatment_effect(g1)

replaceme <- rep(1:3, n)
replaceme <- sample(replaceme[1:n])

X<-X2.model
colnames(X) <- paste("X", 1:ncol(X),sep="_")

lm(Y~treat+X)


# 
# generate.bases(res2.b,
#                res2.b,
#                X,
#                X,
#                id = NULL,
#                replaceme,
#                alwaysinter = res2.b)

replaceme2<- replaceme+sample(c(0,.5),n,T)
replaceme2 <- 2*replaceme2-1
## Construct treatment errors and model
est.run <- est.run2 <- NULL
for(i.outer in 1:3){

d1 <- data.frame(treat,X)
r.treat.errs <- ranger(treat~.,data=d1[replaceme==1,])
treat.errs <- treat - predict(object=r.treat.errs, data=d1)$predictions
 


## Make bases?
# eff.het.bases <- PLCE:::generate.bases(treat.errs, treat.errs, X, X, id=NULL,
                                       # replaceme, alwaysinter=res2.b)

d2 <- data.frame(treat.errs,X)
r.treat <- ranger(abs(treat.errs)~.,data=d2[replaceme==2,])
treat3 <- treat[replaceme==3] - predict(object=r.treat, data=d2[replaceme==3,])$predictions

control.treat <- predict(object=r.treat, data=d2[replaceme==3,])$predictions


d2.2 <- data.frame(treat.errs,X,X*treat.errs)
r.treat2 <- ranger(treat.errs~.,data=d2.2[replaceme==2,])

control.treat2 <- predict(object=r.treat2, data=d2.2[replaceme==3,])$predictions
treat3 <- treat3 - predict(object=r.treat2, data=d2.2[replaceme==3,])$predictions


Y2 <- lm(Y~treat, weights=1*(replaceme==2))$res
dy <- data.frame(Y2,X)
ry <-ranger(Y2~., data=dy[replaceme==2,])
control.y <- predict(object=ry, data=dy[replaceme==3,])$predictions
y3 <- Y[replaceme==3]-predict(object=ry, data=dy[replaceme==3,])$predictions

cov(treat3,y3)/cov(treat3, treat[replaceme==3])
lm(y3[replaceme==3] ~treat3[replaceme==3] + control.y+control.treat+control.treat2)$coef[2]

est.run[i.outer]<-cov(treat3,Y[replaceme==3])/cov(treat3, treat[replaceme==3])
est.run2[i.outer]<-
  lm(Y[replaceme==3] ~treat[replaceme==3] + control.y+control.treat+control.treat2)$coef[2]


lm(Y[replaceme==3] ~treat[replaceme==3]+control.y + control.treat2+control.treat)$coef[2]

if(i.outer < 3) replaceme <- (replaceme+1)%%3+1
}

grf::average_treatment_effect(g1)
mean(est.run)
mean(est.run2)

mean(c(est.run, est.run2))
