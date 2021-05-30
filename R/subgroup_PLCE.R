#' Sparse regression for experimental and observational data.
#'
#' Function for fitting a Bayesian LASSOplus model for sparse models with
#' uncertainty.
#' Function takes a dependent variable and an optional matrix of (pre-treatment)
#' covariates, or a formula and data frame.
#' The function \env{sparsereg} allows for estimation of a broad range of
#' sparse regressions.  The method allows for continuous and binary variables.
#' In experimental data, it can be used for subgroup analysis.
#' In observational data, it can be used
#' in place of a standard regression, especially in the presence of a large
#' number of variables.  The method also adjusts uncertainty estimates when
#' there are repeated observations through using random effects.
#'
#' Additional arguments are used to extend the method to the nonparametric
#' regression setting; see \env{sparseregNP} for details.
#'
#'
#' @importFrom coda mcmc as.mcmc
#' @importFrom tibble as_tibble
#' @importFrom Matrix Matrix t
#'
#' @export




#' @rdname sparsereg
#' #' @param data An optional dataframe with column names matching the terms in \code{formula},
#' when provided.
#'
#'
#'

sparsereg <- function(y, X, treat = NULL,
                      id = NULL,
                      data = NULL,
                      weights = 1,
                      iter = 2000, burnin = 2000, thin = 10,
                      type = "linear",
                      EM = TRUE,
                      prior_covariance = "half-t",
                      beta.init = NULL,
                      unpen = NULL,
                      trim.X = TRUE,
                      verbose = FALSE,
                      scale.type = "none",
                      tol = 1e-3,## 1e-5 in package
                      blockupdate = TRUE,
                      alpha.prior = "parametric",
                      alpha.use = NULL, 
                      sparseregweights=FALSE,
                      iter.initialize = 2, ## 20 in package
                      qfunction.calc = FALSE,
                      edf = FALSE
) {
  
  
  # Preproces ------
  # Remove NA rows and constant columns to get complete (cmp) dataset
  if (length(weights) == 1) weights <- rep(1, length(y))
  
  # Ensure X is named
  name_X<-function(X){
    if(length(colnames(X)) == ncol(X)) return(X)
    if(length(colnames(X)) != ncol(X)) {colnames(X)<-paste("X",1:ncol(X),sep="_"); return(X)}
  }
  
  X<-name_X(X)
  y<-as.vector(y)
  p.re<-length(unique(id))
  
  
  # stop if n still does not match ---
  if (sum(!is.finite(y)) > 0) 
    stop("Please remove missing or infinite values from y")
  if (sum(!is.finite(X)) > 0) 
    stop("Please remove missing or infinite values from X")
  if (NROW(id) > 0 & NROW(id) != length(y)) 
    stop("id and y need to have the same number of observations")
  
  
  
  
  # Declare variables ----
  n <- length(y)
  p <- max(2, ncol(X))
  alpha <- n * log(p)
  alpha <- p + 1
  alpha0 <- p + 1
  alpha1 <- (n * log(p))
  alpha2 <- ((p / n) * alpha1 + (n / p) * alpha0) / (p / n + n / p)
  if (alpha.prior == "oracle") {
    alpha <- alpha1+p#ifelse(EM,alpha1,alpha1 +p)
  }
  if (alpha.prior == "parametric") alpha <- alpha0
  if (alpha.prior == "balanced") alpha <- alpha2
  #if(!EM) alpha<-n*log(p)
  
  if(alpha.prior == "custom") alpha <- alpha.use
  
  conv <- F
  link <- function(x) x
  
  # Formatting type TX, TTX, etc. goes back in here ----
  X0 <- X 
  
  # Eliminate correlated vars from X goes back in here ----
  X <- as.matrix(X)
  scale.back <- apply(X, 2, FUN = function(x) sd_cpp(x))
  X <- as.matrix(apply(X, 2, FUN = function(x) (x-mean(x)) / sd_cpp(x)))
  X <- X[,colMeans(is.infinite(X))<1]
  
  X <- as.matrix(X)
  X <- cbind(1, X)
  X0<-as.matrix(X0[,apply(X0,2,sd)>0])
  X0 <- as.matrix(cbind(1, X0))
  
  X<-as.matrix(X)
  colnames(X0) <- colnames(X)
  colnames(X0)[1] <- colnames(X)[1] <- "(Intercept)"
  y0 <- y
  
  colmap.vec<-rep("main",ncol(X))
  colmap.vec[1]<-"intercept"
  ## Add in Random Effects! ----
  if(p.re>0){
    counts.res <- tapply(id,id, length)    
  }
  
  # Apply weights ----
  y <- y * weights^.5
  X <- apply(X, 2, FUN = function(x) x * weights^.5)
  
  ## Declare crossprods
  XprimeX <- crossprod(X) 
  Xprimey <- crossprod(X,y)
  
  scale.y <- 1#sd_cpp(y)
  y <- y / scale.y
  
  sd.y <-sd_cpp(y)
  tauthresh <- min(max(1e4, sd_cpp(y)^2 * n^3), 1e10)
  
  
  # Format containers; initialize ----
  iter <- 500
  burnin <- 1
  thin <- 1
  # mcmc.obj.init<-matrix(NA, nrow = p, ncol = iter)
  mcmc.obj.init <- matrix(NA, nrow = floor(iter / thin), ncol = ncol(X))
  colnames(mcmc.obj.init) <- colnames(as.matrix(X))
  mcmc.obj.init <- as.mcmc(mcmc.obj.init, start = burnin, end = burnin + iter, thin = thin)
  
  wtsaccept.run <- tausq.run <- wts.run <- beta.run <- mcmc.obj.init
  qfunction.run <- NULL
  
  seq.obj<-list()
  seq.obj[[1]]<- seq.m <- seq(-6, 10, by = 0.025)
  seq.obj[[2]]<- seq.gamma <- seq(-6, 12, by = 0.025)
  seq.gamma<-seq(-5,3,1)
  
  fits.run <- matrix(NA, nrow = n, ncol = floor(iter / thin))
  gammaaccept.run <- gamma.run <- lambda.run <- sigma.sq.run <- rep(NA, floor(iter / thin))
  eps.gamma <- eps.wts <- 1e-6
  Q1 <- -1e10
  
  beta <- beta.last  <- beta.init 
  up.wts<-matrix(1,nrow=3,ncol=length(beta))
  
  if (length(beta.init) != ncol(X)) {
    priormat <- sd_cpp(y)^2 * diag(ncol(X))
    priormat[1, 1] <- 0
    
    beta <- beta.last  <- beta.init <- as.vector(solve_cpp(crossprod(X) + priormat, crossprod(X, y))$beta)
    up.wts<-matrix(1,nrow=3,ncol=length(beta))
    fits.curr<-fits.X <- X%*%beta
    
    
    if (sum(is.na(beta)) > 0) beta[is.na(beta)] <- beta.last[is.na(beta)] <- rnorm(sum(is.na(beta.last))) / 5
    if (sum(beta == 0) > 0) beta[beta == 0] <- beta.last[beta == 0] <- rnorm(sum(beta == 0)) / 5
    
  } else {
    beta <- beta.last <- beta.init
  }
  
  ## Initializing lambda, weights ----
  sigma.sq <- gamma <- lambda <- lambda.sq <- 1
  prec <- ifelse(type == "linear", 1 / sd_cpp(as.vector(y - X %*% beta))^2, 1)
  sigma.sq <- 1 / prec
  tau.sq <- wts <- 1 / abs(beta)
  wts.sq <- wts^2
  #tau.sq[tau.sq > 1e5] <- 1e5
  #tau.sq[tau.sq < 1e-5] <- 1e-5
  tau.sq<-rcppClamp(tau.sq,1e-5,1e5)
  Dtau <- tau.sq#as.matrix(diag(tau.sq))
  Etausqinv <-1/Dtau#1/diag(Dtau)
  if(length(Dtau)==0) Dtau<-as.matrix(1e5)
  
  Dtau[1] <- 1e5
  
  ## Initialize Random Effects ----
  prec.b <- sigma.sq.b <- 1
  means.res<-fits.res <- 0
  
  if(p.re>0) {
    sums.res<-tapply(y-fits.X,id,sum)
    means.res<-sums.res/(counts.res+1)
    means.res<-means.res-mean(means.res)
    fits.res<-means.res[id]
  }
  
  
  
  usesparseregweights <-  TRUE #sparseregweights
  sparseregweights <- TRUE
  ## Main fitting routine for sparsereg ----
  for (i.iter in 1:floor((burnin + iter) / thin)) {
    for (i.thin in 1:thin) {
      sparseregweights <- FALSE
      
      fits.last<-X %*% beta 
      # E step: Update weights and related quantities ----
      gammastar<-gamma
      if(sparseregweights){
        # x.wts<-rep(NA,p)
        # x.wts<-abs(lambda * beta[2:(p+1)] * prec^.5)
        # x.wts <- rcppClamp(x.wts,exp(-5.99),exp(9.95))
        # wts.mat <- cbind(x.wts, gammastar)
        # wts[2:(p+1)] <- apply(wts.mat, 1, updatewts.EM,seq.obj0=seq.obj)
        # if(p.re>0) wts[-c(1:(p+1))] <- 1
        wts<-rep(1,p+1)
        
      } else{
        wts<-rep(1,p+1)
      }
      wts[1]<-1e-5
      
      if (i.iter < 10) wts <- rcppClamp(wts,1e-5,1e3)#pmax(pmin(wts, 1e3), 1e-5)
      beta.trim<- rcppClamp(abs(beta), 1e-5,1e10)
      Ewtsqtausq <- (beta.trim[1:(p+1)] / lambda * sqrt(prec) * wts[1:(p+1)]  + 1 / lambda.sq)
      Etausqinv <- (lambda / beta.trim[1:(p+1)]  * sqrt(sigma.sq) * wts[1:(p+1)] )
      if(length(beta)>1){
        Dtau <- (1 / Etausqinv)
        Dtau[Dtau > 1e5] <- 1e5
        Dtau[1] <- 1e5
        if(i.iter < 4) Dtau <- rcppClamp(Dtau,1e-3,1e3)      
      } else {
        Dtau <-as.matrix(1e5)
      }
      
      conv.res <- TRUE
      if(p.re>0){
        means.res.last<-means.res
        #Dtau<-c(Dtau,rep(sigma.sq.b/sigma.sq,p.re))
        sums.res <- tapply(y-fits.curr,id,sum)
        means.res <- sums.res/(counts.res+sigma.sq/sigma.sq.b)
        means.res <- means.res-mean(means.res)
        fits.res <- means.res[id]
        conv.res <- max(abs(means.res-means.res.last))
      }
      conv.res <- 0
      
      # Zeroes out prior!!!  ----
      # update beta -----
      conv <- F
      beta.last <- beta
      beta2.last <- c(beta.last)
      
      
      up.beta <- updatebeta(y-fits.res, 
                            Dtau, 
                            sigma.sq,
                            TRUE,
                            XprimeX,
                            Xprimey,
                            type, 
                            blockupdate, 
                            beta, 
                            trim = TRUE,
                            i.iter,
                            burnin,
                            tauthresh,
                            sd.y)
      beta<-as.vector(up.beta$beta)
      if(length(beta)!=length(beta2.last)) beta <- beta2.last
      if(sum(!is.finite(beta))>0 ) beta <- beta2.last
    
      fits.curr <- as.vector(X %*% beta )
      beta2 <- c(beta)
      if (i.iter < 10) beta <- sign(beta) * rcppClamp(abs(beta), 1e-5,1e10)#pmax(abs(beta), 1e-5)
      
      # print(max(abs(beta2.last - beta2)) )
      # print(usesparseregweights)
      # print(conv.res)
      
      tryCatch(
      if (max(abs(beta2.last - beta2)) < tol & 
          #usesparseregweights == TRUE & 
          conv.res < tol) {
        beta[(abs(beta) < 1e-3 * sd_cpp(y))]<-0
        conv <- TRUE
      }
      )
      
      # 
      # if (max(abs(beta2.last - beta2)) < tol & usesparseregweights == FALSE & conv.res < tol ) {
      #   usesparseregweights <- TRUE
      #   beta[(abs(beta) < 1e-3 * sd_cpp(y))]<-0
      #   conv <- TRUE
      # }
      # 
      if (i.iter == (burnin + iter) / thin) {
        #cat("Maximum number of EM iter reached\n")
        conv <- TRUE
      }
      # Update precision for main model and random effects----
      
      if(p.re > 0) {
        b.re <- sum(means.res^2)
        b.re <- pmin(pmax(b.re,1e-6),1e6)
      }
      prec <- (
        ( n/2+ p/2 +1)/
          ( sum((y-fits.curr-fits.res)^2) +sum(beta[2:(p+1)]^2 / (Dtau[2:(p+1)])) / 2)
      )
      
      if(p.re>0){
        prec.b <-(
          (p.re/2+1)/
            (b.re/2)   
        )
        # }
        
        sigma.sq<-1/prec
        sigma.sq.b <- 1/prec.b
      } # Close out precision loop
      # Update lambda, lambda.sq ----
      lambda.sq <- (alpha - 1) / (sum(Ewtsqtausq) / 2 + 1)
      lambda <- lambda.sq^.5
      ## Stops lambda from collapsing on zero
      if (i.iter < 10) {
        lambda <- max(1, lambda)
        lambda.sq <- max(1, lambda)
      }
      # Update gamma ----
      if(sparseregweights){
        gamma <- 1
        #gamma <- updategamma.EM(wts, lambda, beta, gamma, prec, approx.EM, n0 = n, seq.gamma)
      } else{
        gamma<-1
      }
      
      
      
      
    } # End thin loop
    
    # Gather output ----
    if (i.iter > ceiling(burnin / thin)) {
      if (EM & conv) {
        beta.run <- beta.run[1, ] <- as.vector(beta)
        names(beta.run) <- colnames(X0)
        beta.run <- as.mcmc(t(beta.run), start = 1, end = i.iter * thin, thin = 1)
        fits.run <- matrix(link(as.vector(X0 %*% beta) +fits.res), ncol = 1)
        sigma.sq.run <- sigma.sq
        wts.run <- wts
        lambda.run <- lambda
        gamma.run <- gamma
        tausq.run <- (Dtau)
        chol.adj <- chol.inv <- GCV <- dof <- NULL
        if(edf)  {
          if(sum(beta!=0)>1){
            XprimeXdf<- as.matrix(XprimeX[beta!=0,])
            XprimeXdf<- as.matrix(XprimeXdf[,beta!=0])
            Dtau[!is.finite(Dtau)] <- 1e-5
            Dtau[is.na(Dtau)]<- 1e-5
            Dtauinv.df <- diag(1/Dtau[beta!=0])
            var.inv <- MASS::ginv(XprimeXdf+Dtauinv.df)
            svd.inv <- svd(var.inv)
            svd.inv$d <-sqrt(svd.inv$d)
            keeps <- (cumsum(svd.inv$d)/sum(svd.inv$d))<0.9#svd.inv$d/svd.inv$d[1]>1e-4
            if(sum(keeps)<=1) {
              chol.inv <- NULL
              chol.adj <-NULL
            } else{
              diag.svd <- (svd.inv$d[keeps])
              chol.inv <-svd.inv$u[,keeps]%*%(diag(diag.svd)%*%t(svd.inv$v[,keeps]))
              #svd2 <- svd(X[,beta.temp!=0]%*%chol.inv)
              #chol.inv<- svd2$v[,svd2$d>0]
            }
            dof <- sum(diag(
              XprimeXdf %*% var.inv
            )
            )
          }
          else {
            chol.inv<-diag(1)
            dof <- 1
          }
          GCV <- (mean((y0-fits.curr-fits.res)^2)/(1-dof/n*(log(n)/2))^2)^.5
          #GCV <- (mean((y0-fits.run)^2)/(1-dof/n)^2)^.5
          
          
          
          
        } 
      }
      if (EM & conv) break
    }
    
  } # End iterations loop -------
  
  ## Use lmer!!!
  
  #if(length(id)>0) means.res <- suppressWarnings(ranef(lmer(y-fits.curr ~(1|id)))$id)
  # print(means.res)
  # names(means.res) <- sort(unique(id))
  
  # Output -----
  output <- list(
    "intercept" = beta.run[,1],
    "coefficients" = beta.run[,colmap.vec=="main"],
    "beta.all" = as.vector(beta.run),
    "REs" = means.res,
    "fitted.values" = fits.curr,
    "sigma.sq" = sigma.sq.run,
    "weights" = weights,
    "dof" = dof,
    "GCV"= GCV,
    "chol.inv" = chol.inv,
    "parameters" = list("wts" = wts.run, 
                        "tausq" = tausq.run,
                        "lambda" = lambda.run, 
                        "gamma" = gamma.run 
    ),
    "colmap"=colmap.vec
  )
  
  return(output)
}
