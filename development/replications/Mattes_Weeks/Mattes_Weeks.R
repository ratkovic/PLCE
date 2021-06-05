
library(readstata13)
library(dplyr)
library(grf)
library(PLCE)
library(xtable)
library(sandwich)
#devtools::load_all('~/Dropbox/InfluenceFunctions/APSRsubmission/02_Resubmission/04a_Replication/Code/PLCE')


runall<-FALSE
Data1<-read.dta13("mattesweeksAJPS.dta")
options(device="quartz")

#reg disapprove1 hawkrap hawkneut doverap doveneut, robust noconst // estimate underlying model so we can compare using lincom command

Data2<- Data1 %>% subset(.,6-hddv1>0)
y0<-((6-Data2$hddv1)-1)/4*100
vars.control<-c("hawk","intl","trust","voted16","birthyr","gender","educ","pid7","ideo5","newsint","pew_religimp")
X<-Data2[,vars.control]
set.seed(1);X.shuff<-apply(X,2,sample)
treat1<-Data2$hawk_t==1 & Data2$rapproche_t==1
treat2<-Data2$hawk_t==1 & Data2$rapproche_t==2
treat3<-Data2$hawk_t==2 & Data2$rapproche_t==1
treat4<-Data2$hawk_t==2 & Data2$rapproche_t==2

sub1<-(Data2$hawk_t)==1
sub2<-(Data2$hawk_t)==2

  if(runall){
    set.seed(1);h1<-plce((y0>50)[sub1],treat=treat1[sub1],X=(X)[sub1,],var.type="HC0",num.fit = 250,printevery = 10)
    set.seed(1);h2<-plce((y0>50)[sub2],treat=treat3[sub2],X=X[sub2,],var.type="HC0",num.fit = 250,printevery = 10)
    
    output<-list("h1"=h1,"h2"=h2)
  save(output,file="mattesHOE.Rda")
  }


set.seed(1)
a1<-causal_forest(Y=(y0>50)[sub2],W=treat3[sub2],X=X[sub2,],num.trees = 4000)
a2<-predict(a1,estimate.variance = TRUE)
d1<-average_treatment_effect(a1)
d1

a1<-causal_forest(Y=(y0>50)[sub1],W=treat1[sub1],X=(X)[sub1,])


plotdata<-data.frame(X[sub2,],a2)



load("mattesHOE.Rda")


lm1.0<-(lm((y0>50)~treat1,sub=Data2$hawk_t==1))

lm1.1<-(lm((y0>50)~treat1+as.matrix(X),sub=Data2$hawk_t==1))


lm2.0<-(lm((y0>50)~treat3,sub=Data2$hawk_t==2))

lm2.1<-(lm((y0>50)~treat3+as.matrix(X),sub=Data2$hawk_t==2))

points1<-c(output$h1$point,lm1.0$coef[2],lm1.1$coef[2])
points2<-c(output$h2$point,lm2.0$coef[2],lm2.1$coef[2])

se.hc<-function(x) vcovHC(x,"HC0")[2,2]^.5
ses1<-c(output$h1$se,se.hc(lm1.0),se.hc(lm1.1))
ses2<-c(output$h2$se,se.hc(lm2.0),se.hc(lm2.1))

uppers1<-points1+1.96*ses1
lowers1<-points1-1.96*ses1

tab1<-rbind(points1,ses1,points2,ses2)*100
rownames(tab1)<-c("Hawks","s.e.","Doves","s.e.")
colnames(tab1)<-c("PLCE","Diff-in-Mean","OLS w Covariates ")
xtable(tab1,digits=2)

#plot(1:3,points1,ylim=range(c(uppers1,lowers1)))
#segments(x0=1:3,x1=1:3,y0=lowers1,y1=uppers1)
