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
res.map <- rnorm(length(unique(ids.map)))
names(res.map)<-(unique(ids.map))
res.true <-  res.map[ids.map]

y<-res.true+rnorm(n)

X <- MASS::mvrnorm(n, rep(0, k), Sig = var.mat)
X<-apply(X,2,scale)

s1<-sparsereg(y,X,id=ids.map)

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


### Figure

obj<-h1


gg_color <- function(n,alpha0=1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100,alpha=alpha0)[1:n]
}

cols<-gg_color(2)

par(mar = c(3, 3,1,0.2), # Dist' from plot to side of page
    mgp = c(2, 0.4, 0), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.01, # Reduce tick length
    xaxs = "i", yaxs = "i", # Remove plot padding
    oma =c(.2,.2,0,0))

obj$diff<-pos_measure(obj,trim=0.025)


plot(make.rank(obj$diff),obj$diff,type="n",col="red",lwd=2,ylim=range(c(obj$diff,0)),
     xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(-5,105))
#abline(c(0,1),col="gray50")
abline(h=0,lwd=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
abline(h=(-100:100)/2,col="white")
abline(v=(0:100)*10,col="white")

mtext(side=1,line=1.7,font=2,text="Percentile",las=0)
mtext(side=2,line=2,font=2,text="Excess Kurtosis",las=0)
axis(1)
axis(2)
lines(make.rank(obj$diff),obj$diff,col=cols[1],lwd=2)

#lines(make.rank(p2$diff),p2$diff,col=cols[2],lwd=2)


legend.txt<-c("Hawks Experiment","Doves Experiment")
legend("topright",legend=legend.txt,col=gg_color(2),lty=1,xpd=T,lwd=2,bg="gray90",bty="n",box.col="gray90")

