#library(PLCE)
library(profvis)
library(microbenchmark)
library(lme4)
library(grf)
# devtools::load_all('~/Dropbox/InfluenceFunctions/APSRsubmission/02_Resubmission/04a_Replication/Code/PLCE')
 

options(device="quartz")
# rm(list=ls())
ran.num<-round(runif(1)*1e8)

## Starting time 430480

## First point estimate: 0.8677028 
theta.run<-se.run<-cover.run<-NULL
n<-500
k<-5

wts<-1
alpha.prior<-"parametric"
id<-sample(letters[1:5],n,T)
beta.init<-NULL
type <- "linear"
addREs<-meanhet<-errhet<-inter<-F

errhet <- meanhet <- addREs<-T
inter <- errhet <- meanhet <- addREs<-T
inter <- addREs <- F

tol<-1e-3
verbose <- EM <- T
        
          #set.seed(1)
          ##Format data
          var.mat   <- diag(k)
          var.mat[var.mat==0] <- 0.5
          # covariates
          
          ids.map<-sample(as.factor(letters[1:20]),n,T)
          ids.map <-sample(as.factor(1:50),n,T)
          res.map <- rnorm(length(unique(ids.map)))
          names(res.map)<-sort(unique(ids.map))
          res.true <-  (res.map[ids.map])*addREs

          X <- MASS::mvrnorm(n, rep(0, k), Sig = var.mat)
          X<-apply(X,2,scale)
          X[,1]<-X[,1]/(mean(X[,1]^2)^.5)
          X<-apply(X,2,FUN=function(x) x/(mean(x^2))^.5)
          if(errhet)	sd.use<-((1+X[,1]+(length(unique(X[,1]))==2) )^2/2)^.5 else sd.use<-1
          errst<-rnorm(n,sd=sd.use)#rt(n,8)*sd.use#
          treat<-(X[,1])+errst +res.true
          errsy<-rnorm(n,sd=1)#rt(n,8)*sd.use#
          if(meanhet) Y<-1+treat*(X[,1]^2)+errsy-res.true else Y<-1+treat+(X[,1]^2)+errsy-res.true
          ##Shift X
          if(inter){
            m1<-exp(-10*(outer(X[,1],X[,1],"-")^2))
            diag(m1)<-0
            oY<-m1%*%X[,1]/rowSums(m1)
            oT<-m1%*%X[,1]/rowSums(m1)
            
            Y<-Y+oY*.5
            treat<-treat-oT*.5
          }
          
          X0.linear<-X2.model<-X
          X2.model[,1]<-X[,1]+.5*X[,2]
          X2.model[,2]<-X[,2]+.5*X[,1]
          ## Rotates X's
          X<-X2.model
          
          y<-treat

          

          
          #ids.map<-sample(ids.map)
          
          X.REs<-model.matrix(~ids.map)
          
          #g1<-causal_forest(cbind(X,X.REs),Y,treat,clusters=as.factor(ids.map))
          g1<-causal_forest(cbind(X),Y,treat)
          average_treatment_effect(g1)
          s1<-PLCE:::sparsereg_GCV(treat,X,id0=NULL)
          s2<-sparsereg::sparsereg(treat,X, EM=T)
          
          mean((X[,1]-s1$fitted.values)^2)^.5
          mean((X[,1]-s2$fitted.values)^2)^.5
          
          #alpha.min<-PLCE:::sparsereg_findalpha(treat,X,id=ids.map)
          #s2<-PLCE:::sparsereg(y,X,id,alpha.prior="custom",alpha.use=alpha.min,edf=T)
          
           #h2<-hoe(Y,treat,X)
         h1<-plce(Y,treat,X,num.fit = 5,printevery=1,var.type = "HC3",fit.interference = F);h1$point;h1$se
          #h2<-plce(Y,treat,X,num.fit = 20,printevery=10);h2$point;h2$se
          
          
        l1<-lme4::lmer(treat~X+(1|ids.map))
        s1$REs
        
        #s1$coef
        t((lme4::ranef(l1)$ids.map))
        #set.seed(1);profvis( {h1<-plce(Y,treat,X,id=ids.map,num.fit = 5);h1$point;h1$se})
        
        # h1<-plce(Y,treat,X,num.fit = 5,printevery=1,var.type = "HC3");h1$point;h1$se
        # h2<-plce(Y,treat,X,num.fit = 5,printevery=1,var.type = "HC3",
                  # fit.interference = F, fit.treatment.heteroskedasticity = T);h2$point;h2$se
       