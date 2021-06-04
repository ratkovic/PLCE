#library(PLCE)
library(profvis)
library(microbenchmark)
devtools::load_all('~/Dropbox/Github/PLCE')


options(device="quartz")
rm(list=ls())
ran.num<-round(runif(1)*1e8)

## Starting time 430480

## First point estimate: 0.8677028 
theta.run<-se.run<-cover.run<-NULL
n<-200
k<-5

wts<-1
alpha.prior<-"parametric"
id<-sample(letters[1:5],n,T)
beta.init<-NULL
type <- "linear"
meanhet<-errhet<-inter<-F
inter <-F
tol<-1e-3
verbose <- EM <- T

#set.seed(1)
##Format data
var.mat   <- diag(k)
var.mat[var.mat==0] <- 0.5
# covariates

ids.map<-as.factor(sample(letters[1:20],n,T))
res.map <- rnorm(length(unique(ids.map)))*0
names(res.map)<-(unique(ids.map))
res.true <-  res.map[ids.map]

y<-res.true+rnorm(n)

X <- MASS::mvrnorm(n, rep(0, k), Sig = var.mat)
X<-apply(X,2,scale)

s1<-sparsereg:::sparsereg(y,X,EM=T,verbose=F,sparseregweights=F,alpha="parametric")
s2<-PLCE:::sparsereg(y,X,alpha.prior="custom",alpha.use=s1$parameters$alpha)
s3<-PLCE:::bayesLasso(y,cbind(1,X),alpha=s1$parameters$alpha)

l1<-lmer(y~(1|ids.map))
t(ranef(l1)$ids)
s1$REs

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
## Rotates X's
X<-X2.model

y<-treat

alpha.max <- max(n*log(ncol(X)), ncol(X)+1)
alpha.min <- min(ncol(X), alpha.max/100)
alpha.list <- seq((alpha.max), alpha.min,length=10)
GCV.out<-sapply(alpha.list, FUN=function(z) sparsereg(y,X, edf=T, alpha.prior = "custom", alpha.use = z )$GCV)
a.min <-alpha.list[GCV.out==min(GCV.out)]
model.out <- sparsereg(y,X, edf=T, alpha.prior = "custom", alpha.use = a.min )


d1<-diag(5)
d1[1,1]<-0

eigen(ginv(d1))
