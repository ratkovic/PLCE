#library(PLCE)
devtools::load_all('~/Dropbox/Github/PLCE')

## This example takes you through an implementation
# of the PLCE model. We first create the sample size and
# number of covariates.
set.seed(1234)
n <- 1000
p <- 5

## Generate covariate matrix
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
X<- apply(X,2,scale)

## Generate random effects
ids.map <- sample(letters[1:10], n, TRUE)
res.map <- rnorm(10)
names(res.map) <- letters[1:10]
res.true <-  (res.map[ids.map])

## Generate the treatment and outcome
treat <- (X[, 1]) +  res.true + rnorm(n)
Y <- 1 + treat * (X[, 1] ^ 2) + res.true + rnorm(n)

## Fit the PLCE model
plce1 <-
  plce(
    y=Y,
    treat=treat,
    X=X, id = ids.map,
    printevery = 1,
    num.fit = 5
  )

## Results: Point estimate, standard error,
##   sensitivity analysis
plce1$point
plce1$se
plce1$sens


obj<-plce1

pos_plot <- function(obj, trim=0.025){

p1<-pos_measure(obj,trim)



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


r1<-diff(range(c(p1$diff,0)))
plot(make.rank(p1$diff),p1$diff,type="n",col="red",lwd=2,ylim=range(c(p1$diff,0)),
     xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(-5,105))
#abline(c(0,1),col="gray50")
abline(h=0,lwd=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
abline(h=(-100:100)/10,col="white")
abline(v=(0:100)*10,col="white")

mtext(side=1,line=1.7,font=2,text="Percentile",las=0)
mtext(side=2,line=2,font=2,text="Excess Kurtosis",las=0)
axis(1)
axis(2)
lines(make.rank(p1$diff),p1$diff,col=cols[1],lwd=2)

#legend.txt<-c("Hawks Experiment","Doves Experiment")
# legend("topright",legend=legend.txt,col=gg_color(2),lty=1,xpd=T,lwd=2,bg="gray90",bty="n",box.col="gray90")

output <-list("ExcessKurtosis"=p1$diff)
invisible(output)
}

pos_plot(plce1)
