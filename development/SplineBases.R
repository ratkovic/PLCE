library(microbenchmark)
library(MASS)
devtools::load_all('~/Dropbox/InfluenceFunctions/APSRsubmission/02_Resubmission/04a_Replication/Code/PLCE')

n<-1000
p<-50

X<-Xt<-mvrnorm(n, rep(0,p),diag(p))
        
treat<-rnorm(n)

## Implementation rn in PLCE

nt <- ncol(Xt)
inter.schedule <- cbind(rep(1:nt, each = nt),
                        rep(1:nt, by = nt))
inter.schedule <-
  inter.schedule[inter.schedule[, 1] > inter.schedule[, 2],]
sample.cor <- 1:length(treat)
sample.cor <- 1:n

cors <- function(inter.schedule){
  apply(
    inter.schedule,
    1,
    FUN = function(x, Xt0 = Xt[sample.cor,], treat0 = treat[sample.cor]) {
      inter.try <- Xt0[, x[1]] * Xt0[, x[2]]
      cor.out <- 0
      if (sd_cpp2(inter.try) > 0)
        # cor.out <- abs(cor(treat0, inter.try, method = "spearman"))
        cor.out <- abs(corSpearman(treat0, inter.try))
       #cor.out <- abs(corPearson(treat0, inter.try))
      
      cor.out
    }
  )
}



