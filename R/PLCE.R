# devtools::load_all('~/Dropbox/InfluenceFunctions/APSRsubmission/02_Resubmission/04a_Replication/Code/PLCE')
#' Sparse regression for experimental and observational data.
#'
#'
#' @param y Vector of the dependent, or outcome variable.
#' @param treat Treatment variable, for which the average partial effect is desired.
#' May be continous, binary, or count.
#' @param X Matrix of covariates. Typical vocabulary would refer to these as
#'   "pre-treatment" covariates.
#' @param inner.boot The number of bootstrapped estimates to do in each subsample.  Default is 30.
#' @param num.fit The number of cross-fitting iterations to run.  Default is 25.
#'            inner.boot = 30,
#' @param var.type The type of variance estimate, passed to \code{vcovHC} in sandwich. Options are
#' \code{HC0}, \code{HC1}, \code{HC2}, \code{HC3}, with \code{HC3} the default.
#'
#'
#' @export
#'
#' @return \describe{
#' \item{point}{The estimated average partial effect.}
#' \item{se}{Standard error of the estimate of the average partial effect.}
#'
#' \item{sens}{Results from the sensitivity analysis.}
#'}
#'
#' @references Ratkovic, Marc.  2020. "Relaxing Assumptions, Improving Inference:
#' Utilizing Machine Learning for Valid Causal Inference."  Working Paper.
#'
#' @rdname plce
#'
#' @importFrom splines bs
#' @importFrom MASS ginv
#' @importFrom sandwich vcovHC
#' @importFrom sensemakr sensemakr
#' @importFrom splines2 bSpline
#' @importFrom splines2 dbs
#' @importFrom RcppEigen fastLm
#' @importFrom ccaPP corSpearman
#' @importFrom ccaPP corPearson
#' @importFrom lme4 lmer
#' @examples
#' \dontrun{
#'
#' }
#'

plce <-
  function(y,
           treat,
           X,
           id=NULL,
           inner.boot = 30,
           num.fit = 25,
           var.type = "HC3",
           sens = TRUE,
           fit.pc = TRUE,
           printevery = 10,
           fit.interference = TRUE, 
           fit.treatment.heteroskedasticity = TRUE
           ) {
    n <- length(y)
    k <- ncol(X)
    y0 <- y
    treat0 <- treat
    X0 <- X
    
    REs<-F
    if(length(id)>0) {
      #X.REs <- model.matrix(~id-1)
      REs<-T
    }
    ## Format inputs  -----
    treat <- as.vector(scale2(treat0))
    y <- as.vector(scale2(lm(y0 ~ treat)$res))
    if(length(id)>0){
      treat <- suppressMessages(scale2(treat-predict(lmer(treat~(1|id)))))
      y <- suppressMessages(scale2(y-predict(lmer(y~(1|id)))))
      
    }
    
    X <- cleanX(X)
    k <- ncol(X)
    
    if(REs) {
      y<-suppressMessages(scale2(y-predict(lmer(y~(1|id)))))
      treat<-suppressMessages(scale2(treat-predict(lmer(treat~(1|id)))))
      
    }
    
    y2.b <- y
    treat.b <- treat
    
    ## Make X.t, X.y ----
    X.y <- X.t <- X
    X <- apply(
      X,
      2,
      FUN = function(z)
        as.vector(scale2(rank(z)))
    )
    
    
    ## Containers ----
    treatpoints.run <- points.run <- NULL
    se.run <- NULL
    sens.run <- NULL
    loo.run <- NULL
    treat.res.run <- matrix(NA, nr = n, nc = (num.fit * 2))
    
    beta.t.init<-beta.y.init <- NULL
    # if(length(id)==0){
    #   beta.t.init <- solve_cpp(crossprod(cbind(1,basest0))+sd_cpp2(treat)^2*diag(ncol(basest0)+1),crossprod(cbind(1,basest0),treat) )
    #   beta.y.init <- solve_cpp(crossprod(cbind(1,basesy0))+sd_cpp2(y)^2*diag(ncol(basesy0)+1),crossprod(cbind(1,basesy0),y0) )
    # }
    cat("#####  Beginning cross-fitting for the average partial effect\n")
    
    for (i.hoe in 1:num.fit) {
      ## Create matrices for fitting ----
      replaceme <- make.replaceme(y,id)
      allbases.obj<-allbases(y,y2.b,treat,treat.b,treat.y,X,id, replaceme,
                             fit.interference, 
                             fit.treatment.heteroskedasticity)

      basest <- basest0 <- allbases.obj$basest0
      basesy <- basesy0 <- allbases.obj$basesy0
      basest2 <- basest2.0 <- allbases.obj$basest2.0
      X.interfy <- allbases.obj$X.interfy
      X.interft <- allbases.obj$X.interft
      replaceme <- allbases.obj$replaceme
      
      
      
      id.temp<-NULL
      #if(length(id)>0) id.temp <- id[replaceme>4]
      ## Propensity model for interference, heteroskedasticity estimated off splits 5-6 ----
      colnames(basest) <- paste("X", 1:ncol(basest), sep = "_")
      st <- sparsereg_GCV(treat[replaceme > 4], X = basest[replaceme > 4, ], id = id.temp
      )
      fits.all <-
        cbind(basest[, intersect(colnames(basest), names(st$coef))]) %*% st$coef[intersect(colnames(basest), names(st$coef))]
      if(length(id.temp)>0){
        fits.all <- fits.all+as.vector(st$REs[id])#X.REs[,names(st$REs)]%*%st$REs
      }
      fits <- fits.all[replaceme < 5 ]
      # m1<-mget(ls())
      # save(m1,file="diagnose.Rda")
      res <- as.vector(scale2(treat[replaceme < 5] - fits))
      
      ## Outcome model for interference estimated off splits 5-6 ----
      colnames(basesy) <- paste("X", 1:ncol(basesy), sep = "_")
      sy <- sparsereg_GCV(y[replaceme > 4], X = basesy[replaceme > 4, ], id = id.temp)
      names.inter.y <- intersect(colnames(basesy), names(sy$coef))
      
      fits.y <- scale2(basesy[, names.inter.y] %*% sy$coef[names.inter.y])
      if(length(id.temp)>0){
        fits.y <- fits.y+as.vector(sy$REs[id])#X.REs%*%sy$REs
      }
      
      resy <- scale2(y - fits.y)
      
      
      
      ## Construct interference bases off this split----
      ## First, the residuals we need
      res1 <- scale2(treat - fits.all)
      colnames(basest2.0) <- paste("X", 1:ncol(basest2.0), sep = "_")
      
      ## No need for REs--already taken out.
      st.2 <-
        sparsereg_GCV(res1[replaceme > 4], apply(
          basest2.0[replaceme > 4, ],
          2,
          FUN = function(z)
            scale2(z * scale2(res1[replaceme > 4]))
        ))
      fits.st.2 <- basest2.0 %*% st.2$coef
      res1.2 <- scale2(res1 - fits.st.2)
      
      ## Make interference bases ----
      if(fit.interference){
      X.resy <- cleanNAs(makebases.interf(X.interfy, resy, replaceme))
      X.rest <- cleanNAs(makebases.interf(X.interft, res1.2, replaceme))
      } else{
        X.resy <- as.matrix(rnorm(n))
        X.rest <- as.matrix(rnorm(n))
      }
      ## Generate objects for cross-fitting ----
      y.temp <- y[replaceme < 5]
      y0.temp <- y0[replaceme < 5]
      treat.temp <- treat[replaceme < 5]
      treat0.temp <- treat0[replaceme < 5]
      replaceme.temp <- replaceme[replaceme < 5]
      id.temp<-NULL
      if(length(id)>0) id.temp <- id[replaceme<5]
      
      # m1<-mget(ls())
      # save(m1,file="diagnose.Rda")
      
      ## Bases to be passed to selection and estimation set ----
      if(fit.treatment.heteroskedasticity){
      bases.t.temp <-
        cbind(basest[replaceme < 5, ], res * basest2[replaceme < 5, ])
      }else{
        bases.t.temp <- basest[replaceme < 5, ]
      }
      bases.y.temp <- basesy0[replaceme < 5, ]
      bases.t.temp <- cleanNAs(bases.t.temp)
      bases.y.temp <- cleanNAs(bases.y.temp)
      

      if(length(X.rest)>0) X.rest.temp <- as.matrix(X.rest)[replaceme < 5, ] else X.rest.temp <- NULL
      if(length(X.resy)>0) X.resy.temp <- as.matrix(X.resy)[replaceme < 5, ] else X.resy.temp <- NULL
      

      h1 <-
        hoe.inner(
          y = y.temp,
          treat = treat.temp,
          y0 = y0.temp,
          treat0 = treat0.temp,
          basest = bases.t.temp,
          basesy = bases.y.temp,
          inner.boot = inner.boot,
          replaceme = replaceme.temp,
          var.type,
          lin.adj = NULL,
          lin.adj.t = NULL,
          sens = sens,
          fit.pc = fit.pc,
          X.rest.temp,
          X.resy.temp,
          beta.t.init,
          beta.y.init,
          id=id.temp,
          fit.interference, 
          fit.treatment.heteroskedasticity
        )
      points.run <- c(h1$point, points.run)
      se.run <- c(h1$se, se.run)
      sens.curr <- h1$sens
      sens.run <- cbind(sens.run, sens.curr)
      #loo.run<-cbind(loo.run,h1$loo)
      rownames(sens.run) <- rownames(h1$sens)
      treatpoints.run <- c(h1$treat, treatpoints.run)
      treat.res.run[replaceme < 5, (2 * i.hoe - 1):(2 * i.hoe)] <-
        h1$treat.res
      
      if (i.hoe %% printevery == 0 | i.hoe %in% c(1, num.fit)) {
        cat("##  Iteration ", i.hoe, "\n")
        cat("#  Point estimates: ", h1$point, "\n")
        cat("#  Standard error estimates: ", h1$se, "\n")
        treatpoints.run <- c(h1$treat, treatpoints.run)
        cat("#  Average point estimate to this point: ",
            hl.mean(points.run),
            "\n\n")
      }
    }
    
    point.est <- hl.mean(points.run)
    var.est <-
      hl.mean(se.run ^ 2) + (hl.mean(points.run ^ 2) - point.est ^ 2) / length(points.run)
    out1 <-
      list(
        "point" = (points.run),
        "se" = (se.run),
        "treatpoint" = treatpoints.run
      )
    out1 <-
      list(
        "point" = point.est,
        "se" = (var.est) ^ .5,
        "jk" = (var(points.run) / 2) ^ .5,
        "treatpoint" = points.run,
        "sens" = sens.run,
        "treat.res" = treat.res.run,
        "loo.hoe" = NULL,
        "loo.lm" = NULL
      )
    
    return(out1)
  }



hoe.inner <-
  function(y,
           treat,
           y0,
           treat0,
           basest,
           basesy,
           inner.boot,
           replaceme,
           var.type,
           lin.adj = NULL,
           lin.adj.t = NULL,
           sens = NULL,
           fit.pc = NULL,
           X.rest,
           X.resy,
           beta.t.init,
           beta.y.init,
           id,
           fit.interference, 
           fit.treatment.heteroskedasticity
  ) {
    
    
    ## Save original inputs ----
    n <- length(y)
    
    ##Declare holders for loop
    sens.run <- point.out <- se.out <- NULL
    treat.res.inner <- matrix(NA, nr = n, nc = 2)
    
    if(fit.interference){
      basest.all <- cbind(basest, X.rest)
      basesy.all <- cbind(basesy, X.resy)
    }else{
      basest.all <- cbind(basest)
      basesy.all <- cbind(basesy)
    }
    colnames(basest.all) <- paste("X", 1:ncol(basest.all), sep = "_")
    colnames(basesy.all) <- paste("X", 1:ncol(basesy.all), sep = "_")
    
    basest.all<-basest.all[,check.cor(basest.all,thresh=0.001)$k]
    basesy.all<-basesy.all[,check.cor(basesy.all,thresh=0.001)$k]
    
    basest.all <- apply(basest.all, 2, scale2)
    basesy.all <- apply(basesy.all, 2, scale2)
    
    y.orthog<-y
    treat.orthog <-treat
    if(length(id)>0) {
      X.dum<-cbind(sample(y),sample(y))
      colnames(X.dum)<-c("X1","X2")
      X.dum<-apply(X.dum,2,FUN=function(z) fastLm(z~cbind(1,y))$res)
      y.orthog <- y-sparsereg_GCV(y,X.dum,id=NULL)$fit
      
      colnames(X.dum)<-c("X1","X2")
      X.dum<-apply(X.dum,2,FUN=function(z) fastLm(z~cbind(1,treat))$res)
      treat.orthog <- treat-sparsereg_GCV(treat,X.dum,id=NULL)$fit
      
    }
    
    
    for (i.outer in 1:2) {
      
      
      ## Select bases off splits 1,2 ----
      basest.all.inner <- orthog.me(treat.orthog, basest.all, weights.lm = (replaceme<3))
      basesy.all.inner <- orthog.me(y.orthog, basesy.all, weights.lm = (replaceme<3))
      
      #basest.all.inner <-svd(apply(basest.all,2,scale2))$u
      basest.all.inner <- basest.all[,check.cor(basest.all[replaceme<3,])$k] #orthog.me(treat.orthog, basest.all, weights.lm = (replaceme<3))
      basesy.all.inner <- basesy.all[,check.cor(basesy.all[replaceme<3,])$k] #orthog.me(y.orthog, basesy.all, weights.lm = (replaceme<3))
      
      id.temp <-NULL
      if(length(id)>0) id.temp<-id[replaceme < 3]
      
      st <- sparsereg(treat[replaceme < 3], basest.all.inner[replaceme < 3, ], id.temp,edf=T, alpha.prior = "parametric")
      select.lin.t <- (abs(st$coef) > 1e-3)
      
      sy <- sparsereg(y[replaceme < 3], basesy.all.inner[replaceme < 3, ], id.temp,edf=T)
      select.lin.y <- (abs(sy$coef) > 1e-3)
      
      soe.t<-soe.y<-NULL
      
      threshind <-
        floor(min(n / 2 * .8, length(st$coef) + length(sy$coef) - 2))
      
      thresh1 <-
        sort(c(abs(st$coef), abs(sy$coef)), dec = T)[threshind]
      select.lin.y[abs(sy$coef) < thresh1] <- F
      select.lin.t[abs(st$coef) < thresh1] <- F
      
      fits.t <- cbind(basest.all.inner) %*% st$coef
      fits.y <- cbind(basesy.all.inner) %*% sy$coef
      
      fits.REs.t <- fits.REs.y <- NULL
      
      
      soe.t<-make.soe(apply(basest.all.inner,2,scale2),id,st,treat,fits.t,fits.REs.t,replaceme)
      #y.temp <-lm(y~soe.t,weights=1*(replaceme < 3) )$res
      soe.y<-make.soe(apply(basesy.all.inner,2,scale2),id,sy,y,fits.y,fits.REs.y,replaceme)

            if (fit.pc | TRUE) {
        Xmat.adjust <- as.matrix(cbind(
            # basest.all.inner[replaceme > 2,],
            # basesy.all.inner[replaceme > 2,],
            fits.t[replaceme > 2],
            fits.y[replaceme > 2],
            soe.t,
            soe.y,
          NULL
        ))
        
        
        if (length(Xmat.adjust) == 0) {
          Xmat.adjust <- matrix(rnorm(sum(replaceme > 2)))
          Xmat.adjust <-
            matrix(lm(Xmat.adjust ~ y[replaceme > 2] + treat[replaceme > 2])$res)
        }
        
      } else {
        Xmat.adjust <- as.matrix(cbind(lin.adj[replaceme > 2, ]))
        
      }
      
      if(F){
        cat("#########################\n")
        cat("Dimension Check: Candidate Bases\n")
        cat(ncol(basesy.all.inner)+ncol(basest.all.inner),"\n")
        cat("Dimension Check:: Selected Bases\n")
        cat(ncol(Xmat.adjust),"\n")
      }
      
      if(length(id)>0){
        id.2 <- id[replaceme>2]
        y2 <- lm(y0[replaceme > 2]~Xmat.adjust)$res
        treat2 <- lm(treat0[replaceme > 2]~Xmat.adjust)$res
        res.y.temp <- suppressMessages(predict(lmer(y2~(1|id.2))))
        res.t.temp <- suppressMessages(predict(lmer(treat2~(1|id.2))))
        
        
        Xmat.adjust<-cbind(Xmat.adjust,res.y.temp,res.t.temp)
        lm1 <-
          lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      } else {
        lm1 <-
          lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      }
      lm1$coef[2]
      keeps <- which(!is.na(lm1$coef[-c(1:2)]))
      Xmat.adjust <- cleanNAs(Xmat.adjust[, keeps])
      lm1 <- lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      lm.inf <- lm(treat0[replaceme > 2] ~ Xmat.adjust)
      name.use <- names(lm1$coef)[2]
      
      treat.res.inner[replaceme > 2, i.outer] <-
        lm(treat0[replaceme > 2] ~ Xmat.adjust)$res
      sens1 <- NULL
      sens1$sensitivity_stats <- rnorm(3)
      if (sens) {
        y2 <- c(y0[replaceme > 2], y0[replaceme > 2], y0[replaceme > 2])
        treat2 <-
          c(treat0[replaceme > 2], treat0[replaceme > 2], treat0[replaceme > 2])
        Xmat.2 <- rbind(Xmat.adjust, Xmat.adjust, Xmat.adjust)
        lm.s <- lm(y2 ~ treat2 + Xmat.2)
        name.use <- names(lm.s$coef)[2]
        sens1 <- (sensemakr(lm.s, name.use)$sensitivity_stats)
      }
      sens1 <- as.numeric(unlist(sens1)[-1])
      sens.run <- cbind(sens.run, sens1)
      #rownames(sens.run)<-names(sensemakr(lm1,name.use)$sensitivity_stats)[-1]
      point.out[i.outer] <- lm1$coef[2]
      se.out[i.outer] <- vcovHC(lm1, var.type)[2, 2] ^ .5
      n.se <- length(lm1$res)
      k.se <- sum(!is.na(lm1$coef))
      se.adj <- 1 / (3 - k.se / n.se)
      ## Switch indices for crossfit ----
      replaceme <- 5 - replaceme
    }# Close outer loop
    output <-
      list(
        "point" = (point.out),
        "se" = (se.out ^ 2 / 3) ^ .5,
        "Ut" = NULL,
        "Uy" = NULL,
        "U.adjust" = NULL,
        "subset" = 5 - replaceme,
        "treatpoint" = NULL,
        "sens" = sens.run,
        "treat.res" = treat.res.inner,
        "loo" = NULL
      )
    return(output)
  }


hoe.inner_cholesky <-
  function(y,
           treat,
           y0,
           treat0,
           basest,
           basesy,
           inner.boot,
           replaceme,
           var.type,
           lin.adj = NULL,
           lin.adj.t = NULL,
           sens = NULL,
           fit.pc = NULL,
           X.rest,
           X.resy,
           beta.t.init,
           beta.y.init,
           id,
           fit.interference, 
           fit.treatment.heteroskedasticity
  ) {
    
    
    ## Save original inputs ----
    n <- length(y)
    
    ##Declare holders for loop
    sens.run <- point.out <- se.out <- NULL
    treat.res.inner <- matrix(NA, nr = n, nc = 2)

    if(fit.interference){
    basest.all <- cbind(basest, X.rest)
    basesy.all <- cbind(basesy, X.resy)
    }else{
      basest.all <- cbind(basest)
      basesy.all <- cbind(basesy)
    }
    colnames(basest.all) <- paste("X", 1:ncol(basest.all), sep = "_")
    colnames(basesy.all) <- paste("X", 1:ncol(basesy.all), sep = "_")
    
    basest.all<-basest.all[,check.cor(basest.all,thresh=0.001)$k]
    basesy.all<-basesy.all[,check.cor(basesy.all,thresh=0.001)$k]
    
    basest.all <- apply(basest.all, 2, scale2)
    basesy.all <- apply(basesy.all, 2, scale2)
    
    y.orthog<-y
    treat.orthog <-treat
    if(length(id)>0) {
      X.dum<-cbind(sample(y),sample(y))
      colnames(X.dum)<-c("X1","X2")
      X.dum<-apply(X.dum,2,FUN=function(z) fastLm(z~cbind(1,y))$res)
      y.orthog <- y-sparsereg_GCV(y,X.dum,id=NULL)$fit
      
      colnames(X.dum)<-c("X1","X2")
      X.dum<-apply(X.dum,2,FUN=function(z) fastLm(z~cbind(1,treat))$res)
      treat.orthog <- treat-sparsereg_GCV(treat,X.dum,id=NULL)$fit
      
    }
    
    
    for (i.outer in 1:2) {
      
      
      ## Select bases off splits 1,2 ----
      basest.all.inner <- orthog.me(treat.orthog, basest.all, weights.lm = (replaceme<3))
      basesy.all.inner <- orthog.me(y.orthog, basesy.all, weights.lm = (replaceme<3))
      
      basest.all.inner <- basest.all[,check.cor(basest.all[replaceme<3,])$k] #orthog.me(treat.orthog, basest.all, weights.lm = (replaceme<3))
      basesy.all.inner <- basesy.all[,check.cor(basesy.all[replaceme<3,])$k] #orthog.me(y.orthog, basesy.all, weights.lm = (replaceme<3))
      
      id.temp <-NULL
      if(length(id)>0) id.temp<-id[replaceme < 3]

      st <- sparsereg_GCV(treat[replaceme < 3], basest.all.inner[replaceme < 3, ], id.temp)
      select.lin.t <- (abs(st$coef) > 1e-3)

      sy <- sparsereg_GCV(y[replaceme < 3], basesy.all.inner[replaceme < 3, ], id.temp)
      select.lin.y <- (abs(sy$coef) > 1e-3)
      
      soe.t<-soe.y<-NULL
      
      threshind <-
        floor(min(n / 2 * .8, length(st$coef) + length(sy$coef) - 2))
      
      thresh1 <-
        sort(c(abs(st$coef), abs(sy$coef)), dec = T)[threshind]
      select.lin.y[abs(sy$coef) < thresh1] <- F
      select.lin.t[abs(st$coef) < thresh1] <- F
      
      fits.t <- cbind(basest.all.inner) %*% st$coef
      fits.y <- cbind(basesy.all.inner) %*% sy$coef
      
      fits.REs.t <- fits.REs.y <- NULL
      
      soe.y<-make.soe(basesy.all.inner,id,sy,y,fits.y,fits.REs.y,replaceme)
      soe.t<-make.soe(basest.all.inner,id,st,treat,fits.t,fits.REs.t,replaceme)
      
      if (fit.pc | TRUE) {
        Xmat.adjust <- as.matrix(cbind(#basest.all[replaceme > 2, select.lin.t],
          fits.t[replaceme > 2],
          fits.y[replaceme > 2],
          soe.t,
          soe.y
        ))
        
        
        if (length(Xmat.adjust) == 0) {
          Xmat.adjust <- matrix(rnorm(sum(replaceme > 2)))
          Xmat.adjust <-
            matrix(lm(Xmat.adjust ~ y[replaceme > 2] + treat[replaceme > 2])$res)
        }
        
      } else {
        Xmat.adjust <- as.matrix(cbind(lin.adj[replaceme > 2, ]))
        
      }
      
      if(F){
        cat("#########################\n")
        cat("Dimension Check: Candidate Bases\n")
        cat(ncol(basesy.all.inner)+ncol(basest.all.inner),"\n")
        cat("Dimension Check:: Selected Bases\n")
        cat(ncol(Xmat.adjust),"\n")
      }
      
      if(length(id)>0){
        id.2 <- id[replaceme>2]
        y2 <- lm(y0[replaceme > 2]~Xmat.adjust)$res
        treat2 <- lm(treat0[replaceme > 2]~Xmat.adjust)$res
        res.y.temp <- suppressMessages(predict(lmer(y2~(1|id.2))))
        res.t.temp <- suppressMessages(predict(lmer(treat2~(1|id.2))))
        

        Xmat.adjust<-cbind(Xmat.adjust,res.y.temp,res.t.temp)
        lm1 <-
          lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      } else {
        lm1 <-
          lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      }
      lm1$coef[2]
      keeps <- which(!is.na(lm1$coef[-c(1:2)]))
      Xmat.adjust <- cleanNAs(Xmat.adjust[, keeps])
      lm1 <- lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      lm.inf <- lm(treat0[replaceme > 2] ~ Xmat.adjust)
      name.use <- names(lm1$coef)[2]
      
      treat.res.inner[replaceme > 2, i.outer] <-
        lm(treat0[replaceme > 2] ~ Xmat.adjust)$res
      sens1 <- NULL
      sens1$sensitivity_stats <- rnorm(3)
      if (sens) {
        y2 <- c(y0[replaceme > 2], y0[replaceme > 2], y0[replaceme > 2])
        treat2 <-
          c(treat0[replaceme > 2], treat0[replaceme > 2], treat0[replaceme > 2])
        Xmat.2 <- rbind(Xmat.adjust, Xmat.adjust, Xmat.adjust)
        lm.s <- lm(y2 ~ treat2 + Xmat.2)
        name.use <- names(lm.s$coef)[2]
        sens1 <- (sensemakr(lm.s, name.use)$sensitivity_stats)
      }
      sens1 <- as.numeric(unlist(sens1)[-1])
      sens.run <- cbind(sens.run, sens1)
      #rownames(sens.run)<-names(sensemakr(lm1,name.use)$sensitivity_stats)[-1]
      point.out[i.outer] <- lm1$coef[2]
      se.out[i.outer] <- vcovHC(lm1, var.type)[2, 2] ^ .5
      n.se <- length(lm1$res)
      k.se <- sum(!is.na(lm1$coef))
      se.adj <- 1 / (3 - k.se / n.se)
      ## Switch indices for crossfit ----
      replaceme <- 5 - replaceme
    }# Close outer loop
    output <-
      list(
        "point" = (point.out),
        "se" = (se.out ^ 2 / 3) ^ .5,
        "Ut" = NULL,
        "Uy" = NULL,
        "U.adjust" = NULL,
        "subset" = 5 - replaceme,
        "treatpoint" = NULL,
        "sens" = sens.run,
        "treat.res" = treat.res.inner,
        "loo" = NULL
      )
    return(output)
  }




##  Function for making basis functions
makebases <- function(treat, X, het, SIS.use = NULL, replaceme,alwaysinter=NULL) {
  X <- X[, check.cor(X, thresh = 0.000001)$k]
  colnames(X) <- paste("X", 1:ncol(X), sep = "_")
  ## Create bases ----
  Xr <- apply(X, 2, rank)
  expand <- apply(
    Xr,
    2,
    FUN = function(x)
      length(unique(x)) > 2
  )
  #Xt<-lapply(data.frame(Xr[,expand]),make.spline,name.var="z",degree=c(3,5))
  #Xt<-matrix(unlist(Xt),nr=nrow(X))
  #Xt<-cbind(X,Xt)
  Xt <- NULL
  for (i in 1:ncol(X)) {
    Xt <- cbind(Xt, bs.me(X[, i]))
  }
  Xt <- Xt[, apply(Xt, 2, sd_cpp2) > 0]
  Xt <- Xt[, check.cor(Xt, thresh = 1e-5)$k]
  nt <- ncol(Xt)
  inter.schedule <- cbind(rep(1:nt, each = nt),
                          rep(1:nt, by = nt))
  inter.schedule <-
    inter.schedule[inter.schedule[, 1] > inter.schedule[, 2],]
  sample.cor <- 1:length(treat)
  sample.cor <- which(replaceme > 4)
  if(length(sample.cor)>1000){
    sample.cor <- sample(sample.cor, 1000, FALSE)
  }
  ## Old code using apply
  # if(length(alwaysinter)==0) alwaysinter <- rep(1,n)
  # cors <-
  #   apply(
  #     inter.schedule,
  #     1,
  #     FUN = function(x, Xt0 = Xt[sample.cor,], treat0 = treat[sample.cor],alwaysinter0=alwaysinter[sample.cor]) {
  #       inter.try <- Xt0[, x[1]] * Xt0[, x[2]]*alwaysinter0
  #       cor.out <- 0
  #       if (sd_cpp2(inter.try) > 0)
  #         # cor.out <- abs(cor(treat0, inter.try, method = "spearman"))
  #         cor.out <- abs(corSpearman(treat0, inter.try))
  #         #cor.out <- abs(corPearson(treat0, inter.try))
  # 
  #       cor.out
  #     }
  #   )

  
  ## Rcpp implementation
  cors <- corbases(treat=treat[sample.cor], X = Xt[sample.cor,],inter.schedule,1:length(sample.cor))$cors
  ## Sure screen ----
  n <- length(treat)
  if (length(SIS.use) == 0) {
    SIS.use <- floor(25 * (1 + n ^ (1 / 5)))
    SIS.use <- min(SIS.use, 400)
    SIS.use <- max(SIS.use, 50)
    #SIS.use<-min(SIS.use,n/4)
    #SIS.use<-floor(SIS.use/2)
  }
  
  keeps <- sort(cors, dec = T, ind = T)$ix[1:SIS.use]
  t1 <-
    apply(
      inter.schedule[keeps,],
      1,
      FUN = function(x, Xt0 = Xt, treat0 = treat) {
        #inter.try<-Xt0[,x[1]]*Xt0[,x[2]]
        #lm(Xt0[,x[1]]*Xt0[,x[2]]~Xt0[,x[1]]+Xt0[,x[2]])$res
        Xt0[, x[1]] * Xt0[, x[2]]
      }
    )
  t1 <- cbind(X, t1)
  t1 <- t1[, check.cor(t1, thresh = 1e-5)$k]
  colnames(t1) <- paste("X", 1:ncol(t1), sep = "_")
  
  output <- list("bases" = cbind(t1))
  return(output)
}


#' Double Machine Learning
#' Run Double Machine Learning Method with continuous treatment
#'
#'
#' @param Y Outcome vector.
#' @param D Treatment variable (scalar).
#' @param X Covariates.
#' @param method_out Method used to predict the outcome. Random forest (\code{"RF"}) and LASSO (\code{"LASSO"}) are allowed.
#' @param method_treat See \code{method_out} options.
#' @examples
#'      rfrf_ce <- DML(Y = obs, D = treat, X = X, method_out = 'RF', method_treat = 'RF')
#'      rfls_ce <- DML(Y = obs, D = treat, X = X, method_out = 'RF', method_treat = 'LASSO')
#'      lsrf_ce <- DML(Y = obs, D = treat, X = X, method_out = 'LASSO', method_treat = 'RF')
#'      lsls_ce <- DML(Y = obs, D = treat, X = X, method_out = 'LASSO', method_treat = 'LASSO')
DML <-
  function(Y,
           D,
           X,
           n_split = 2,
           method_out = "RF",
           method_treat = "RF",
           fit.lm = "IV",
           ...) {
    # sample splitting stage
    idx <- sample(1:n_split, size = length(Y), replace = TRUE)
    res_point <- rep(NA, n_split - 1)
    save_diff <- list()
    
    Y.ss <- treat.ss <- NULL
    for (split in 1:n_split) {
      # get training samples' index
      train_idx <- idx == split
      
      # predict outcome
      Y.ss[idx == (split)] <-
        Y_diffed <- DML_predict(Y, X, method_out, train_idx, ...)
      
      # predict treatment
      treat.ss[idx == (split)] <-
        D_diffed <- DML_predict(D, X, method_treat, train_idx, ...)
      save_diff[[split]] <- cbind(Y_diffed, D_diffed, D[train_idx])
      
      if (fit.lm == "IV") {
        out <-
          summary(ivreg(Y_diffed ~ D[idx == split],  ~ D_diffed))$coef[2, 1:2]
      } else {
        # OLS regress (Y-Y_hat) on (D - D_hat)
        out <- summary(lm(Y_diffed ~ D_diffed))$coef[2, 1:2]
      }
      
      res_point[split] <- out[1]
      # res_se[split]    <- out[2]
    }
    
    res_point <- mean(res_point)
    
    ## compute se here
    diff_mat <- do.call("rbind", save_diff)
    Vhat     <- diff_mat[, 2]
    
    if (fit.lm == "IV") {
      Zhat <- diff_mat[, 1] - diff_mat[, 3] * res_point
    } else {
      Zhat <- diff_mat[, 1] - diff_mat[, 2] * res_point
    }
    res_se   <-
      sqrt(mean(Vhat ^ 2 * Zhat ^ 2)) / (mean(Vhat ^ 2) * sqrt(length(Y)))
    
    return(list('point' = res_point, 'se' = res_se))
    
  }


DML_predict <-
  function(outcome, predictor, method, train_idx, ...) {
    if (method == "RF") {
      # double machine learning
      args <- list(...)
      # train and predict the outcome with covariates X
      rfout <- randomForest(x = predictor[!train_idx,],
                            y = outcome[!train_idx],
                            ntrees = args[['n_trees']])
      predict_diffed <-
        outcome[train_idx] - predict(rfout, newdata = predictor[train_idx,])
      
    } else if (method == "LASSO") {
      #lsout <- cv.glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1, nfolds = 5)
      lsout <- rlasso(outcome[!train_idx] ~ predictor[!train_idx, ])
      # lsout <- glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1)
      predict_diffed <-
        outcome[train_idx] - predict(lsout, newdata = predictor[train_idx,])
    } else if (method == "SL") {
      #lsout <- cv.glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1, nfolds = 5
      p1 <- predictor[!train_idx, ]
      p2 <- predictor[train_idx, ]
      lsout <-
        SuperLearner(
          outcome[!train_idx] ,
          data.frame(p1),
          newX = data.frame(p2),
          SL.library = c(
            "SL.glm",
            "SL.glmnet",
            "SL.randomForest",
            "SL.loess",
            "SL.polymars",
            "SL.mean"
          )
        )
      # lsout <- glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1)
      predict_diffed <- outcome[train_idx] - lsout$SL.predict
    } else if (method == "KRLS") {
      p1 <- predictor[!train_idx, ]
      p2 <- predictor[train_idx, ]
      lsout <- krls(y = outcome[!train_idx] , X = (p1))
      predict_diffed <-
        outcome[train_idx] - predict(lsout, newdata = p2)$fit
    } else {
      stop("Proper method required!\n")
    }
    
    return(predict_diffed)
  }



cleanX <- function(X) {
  X <- as.matrix(X[, apply(X, 2, sd_cpp2) > 0])
  if (ncol(X) > 1)
    X <- X[, check.cor(X)$k]
  X <- apply(X, 2, rank)
  X <- as.matrix(apply(X, 2, scale2))
  colnames(X) <- paste("X", 1:ncol(X), sep = "_")
  #return(X)
  #X<-cbind(X,svd(X)$u)
  #X <- cbind(X, make.twosvd(X))#,svd(X)$u)
  X <- X[, check.cor(X)$k]
  X <- as.matrix(apply(X, 2, rank))
  
  #X <- cbind(X,svd(X)$u)
  X <- as.matrix(apply(X, 2, scale2))
  colnames(X) <- paste("X", 1:ncol(X), sep = "_")
  return(X)
  
  #X<-make.interXXsame(X)$X
  
  X <- as.matrix(apply(X, 2, scale))
  X <- as.matrix(X[, apply(X, 2, sd_cpp2) > 0])
  if (ncol(X) > 1)
    X <- X[, check.cor(X, thresh = 0.001)$k]
  X <- as.matrix(apply(X, 2, scale))
  colnames(X) <- paste("X", 1:ncol(X), sep = "_")
  X
}

make.replaceme_old <- function(x) {
  # n <- length(x)
  # out <- rep(sample(1:6), length = n)
  # return(sample(out))
  out <-NULL
  for(i in 1:(ceiling(n/6) +2)){
    out[6*(i-1)+1:6]<-sample(1:6)
  }
  return(out)
  # out<-out[1:n]
  # ind.sort <- sort(x, ind = T)$ix
  # return(out[ind.sort])
}

make.replaceme <- function(x, id) {
  n <- length(x)
  out <-NULL
  if(length(id)==0){
    for(i in 1:(ceiling(n/6) +2)){
      out[6*(i-1)+1:6]<-sample(1:6)
    }
    out<-out[1:n]
  } else{
    unique.id <- unique(id)
    for(i in unique.id){
      out.temp <- which(id==i)
      length.temp <- length(out.temp)
      out[out.temp]<-sample(rep(sample(1:6), ceiling(length.temp/6)) )[1:length.temp]
    }

  }
  return(out)
  
  # ind.sort <- sort(x, ind = T)$ix
  # return(out[ind.sort])
}




orthog.me <- function(y,
                      X,
                      mult = NULL,
                      weights.lm = NULL) {
  if (length(X) > 0) {
    X <- as.matrix(cleanNAs(X))
  } else{
    return(NULL)
  }
  
  
  if (length(weights.lm) == 0)
    weights.lm <- rep(1, length(y))
  X0 <- X
  X2 <- X * 0
  X <- apply(X, 2, scale2)
  y <- scale2(y)
  if (length(mult) > 0)
    X <- apply(
      X,
      2,
      FUN = function(z)
        z * scale2(mult)
    )
  for (i in 1:ncol(X)) {
    covs.curr <- abs(cov(X, y))
    covs.curr[covs.curr == 0] <- -1
    which.curr <- sort(covs.curr, dec = T, ind = T)$ix[1]
    x.temp <- X[, which.curr]
    X[, which.curr] <- 0
    X2[, i] <- x.temp
    X <- fastres_cpp(X = cbind(1, x.temp),
                     y = X,
                     w = weights.lm)$res
  }
  if (length(mult) > 0)
    X2 <- apply(
      X2,
      2,
      FUN = function(z)
        z / scale2(mult)
    )
  X2[, apply(X2, 2, sd_cpp2) > 0.001]
}

hl.mean <- function(x) {
  x <- x[!is.na(x)]
  diff1 <- as.vector(outer(x, x, "+") / 2)
  median_cpp(c(diff1,x))
}





trimbases.boot_quick <-
  function(y,
           basesy0,
           trimboot.num = 30,
           wts.use = NULL,
           id = NULL) {
    n <- length(y)
    if (length(wts.use) == 0)
      wts.use <- rep(1, n)
    
    #alpha.min<-sparsereg_findalpha(y,basesy0,id=id)
    #s1 <- sparsereg(y,basesy0,id,alpha.prior="custom",alpha.use=alpha.min,edf=T)
    
    keep.cols<- rep(0,ncol(basesy0))
    for(i in 1:1){
      samp.curr <- 1:n
      #s1 <- sparsereg(y[samp.curr],basesy0[samp.curr,],id,alpha.prior="balanced",edf=T)
      #alpha.min<-sparsereg_findalpha(y,basesy0,id=id)
      s1 <- sparsereg(y,basesy0,id,alpha.prior="parametric",edf=T)
      
      
      select.vals<-apply(n*abs(crossprod(s1$chol.inv)),1,max)[-1]
      keep.cols <- keep.cols+((select.vals > 0.01)|(abs(s1$coefficients)>0.01))
    }
    keep.cols <-1*(keep.cols>0)
    if(sum(keep.cols)==0) keep.cols[sort(select.vals,decreasing = TRUE,ind=T)$ix[1:2]]<-TRUE
    return(as.matrix(basesy0[,which(keep.cols==1)]))
  }



trimbases.boot <-
  function(y,
           basesy0,
           trimboot.num = 10,
           wts.use = NULL,
           id = NULL,
           replaceme) {
    n <- length(y)
    if (length(wts.use) == 0)
      wts.use <- rep(1, n)
    coefsy.samp <- NULL
    basesy0 <- apply(basesy0,2,scale2)
    y<- scale2(y)
    basesy0int <- cbind(1, basesy0)
    beta.init <- NULL
    
    
    XpX <- crossprod(basesy0int) + sd_cpp2(y) ^ 2 * diag(ncol(basesy0int))
    Xpy <- crossprod(basesy0int, y)
    beta.init <- solve_cpp(XpX, Xpy)$beta
    
    indsamp <- which(replaceme > 4)
    numsamp <- length(indsamp)
    numsamp <- min(numsamp,1000)
    sy0 <-
      sparsereg_GCV(y,
                    basesy0)
    trimboot.num <- rcppClamp(max(ceiling(sy0$dof),2)*5,10,50)
    
    for (i.samp in 1:trimboot.num) {

      #samp.curr <- which(sample(1:3, n, T) == 1)

      samp.curr <- sample(indsamp, numsamp, TRUE)
      basesy0.boot <- basesy0[samp.curr,]
      sds.boot <- apply(basesy0.boot,2,sd_cpp2)
      basesy0.boot[,sds.boot<0.01] <- rnorm(length(samp.curr)*sum(sds.boot<0.01))
      sy0 <-
        sparsereg_GCV(y[samp.curr],
                  basesy0.boot)
      # ,
      #             weights = wts.use[samp.curr],
      #             beta.init = beta.init,
      #             sparseregweights = TRUE
      #   )
      coefsy.samp <- rbind(abs(sy0$coef), coefsy.samp)
    }
    
    cutoff.sis <-
      mean(sapply(
        1:500,
        FUN = function(x)
          max(abs(rnorm(trimboot.num))) / n ^ .5
      )) * 2
    SIS.keep.y <- apply(coefsy.samp, 2, max) > cutoff.sis
    basesy0 <- basesy0[, SIS.keep.y]
    as.matrix(basesy0)
  }



## Generate pairwise SVD ----

make.twosvd <- function(X) {
  X0 <- apply(X, 2, rank)
  X0 <- apply(X0, 2, scale2)
  col.sched <- cbind(rep(1:ncol(X0), each = ncol(X0)),
                     rep(1:ncol(X0), times = ncol(X0)))
  col.sched <- t(apply(
    col.sched,
    1,
    FUN = function(z)
      sort(z)
  ))
  col.sched <- unique(col.sched, MAR = 1)
  col.sched <- col.sched[col.sched[, 1] != col.sched[, 2],]
  X.svd <- NULL
  for(i in 1:nrow(col.sched)) {
    inter.temp <- as.vector(apply(X0[,col.sched],1,prod))
     X.svd<-cbind(X.svd,svd(X0[,col.sched])$u[,1:2])
  }
  #svd.X <- svd(X)
  #spec <- cumsum(svd(X)$d) / sum((svd(X)$d))
  #keeps <- rep(T, ncol(X))
  #keeps[keeps > 0.9] <- F
  #X.svd <- cbind(X, apply(svd.X$u, 2, scale)[, keeps])
  X.svd
}

## Construct bases ----
generate.bases <- function(y2.b, y, basesy0, X, id=NULL, replaceme,alwaysinter=NULL) {
  n <- length(y2.b)
  X <- basesy0
  if(length(id)>0) {
    X.dum<-cbind(sample(y),sample(y))
    colnames(X.dum)<-c("X1","X2")
    X.dum<-apply(X.dum,2,FUN=function(z) fastLm(z~cbind(1,y))$res)
    y<-y2.b <- y2.b-sparsereg(y2.b,X.dum,id=NULL)$fit
  }
  #makebases <- function(treat, X, het, SIS.use = NULL, replaceme)
  ## Stage 1 SIS
  basesy0 <-
    suppressMessages(makebases(y2.b, basesy0, het = F, SIS.use=NULL,replaceme,alwaysinter)$bases)
  keeps.try <- check.cor(basesy0, thresh = 0.0000001)$k
  keeps.try[1:ncol(X)] <- T
  #basesy0<-basesy0[,keeps.try]
  
  basesy <- cbind(X, basesy0[, keeps.try])
  keeps.try <- check.cor(basesy, thresh = 0.000001)$k
  keeps.try[1:ncol(X)] <- T
  basesy <- as.matrix(basesy[, keeps.try])
  basesy <- orthog.me(y, basesy)
  basesy <- apply(basesy, 2, scale2)
  
  #wts1<-(t(basesy)%*%scale2(y))^2
  #basesy<-basesy[,sort(wts1,dec=T,ind=T)$ix]
  
  # wts.ls <- as.vector(t(basesy)%*%scale2(y))
  # bases.ls <-t(apply(basesy,1,FUN=function(z) cumsum(z*wts.ls)))
  # return(basesy[,as.vector(cor(y,bases.ls))< 0.9])
  # 
  #  
  # cs1 <- cumsum(wts1)/sum(wts1)
  # 
  # return(basesy[,cs1<0.9])
  
  basesy <- trimbases.boot(y, basesy, replaceme=replaceme)
  basesy <- cbind(X, basesy)
  basesy <- basesy[, check.cor(basesy, thresh = 1e-4)$k]
  basesy
}


## Generate interference bases ----
## Structure:
## 1) generate X.interf for t, y off full sample (use rot bandwidth)
## 2) Use X.interf for each
## Generates bases, calls estimate bases
##
generate.Xinterf <- function(resy, X, treat, replaceme) {
  n <- length(resy)
  
  ## Make bases
  X.interf <- NULL
  for (i.interf in 1:ncol(X)) {
    x.curr <- X[, i.interf]
    X.interf <- cbind(X.interf, bs.me(x.curr))
  }
  
  X.interf <- X.interf[, apply(X.interf, 2, sd_cpp2) > 0.0001]
  X.interf <-
    cbind(1, X.interf[, check.cor(X.interf, thresh = 0.001)$k])
  ## Typo from earlier version!
  #rot.est <- -log(1.06 * length(n) ^ .2 * (1 / 3) ^ .2)
  #rot.est <- -log(1.06 * (n) ^ .2 * (1 / 3) ^ .2)
  rot.est <- log(1.06 * (n) ^ -.2 * (1 / 3) ^ -.2)
  
  
  X.interf.outer <- X.interf.mat <- X.interf
  if (length(treat) > 0)
    X.interf.outer <- cbind(X.interf)
  
  
  ## add intercept
  ncX.mat <- ncol(X.interf.mat)
  ncX.outer <- ncol(X.interf.outer)
  
  inter.schedule.Xint <-
    cbind(rep(1:ncX.mat, by = ncX.outer),
          rep(1:ncX.outer, each = ncX.mat))
  X.interf.outer <- apply(X.interf.outer, 2, scale2)
  X.interf.mat <- apply(X.interf.mat, 2, scale2)
  resy.2 <- scale2(resy)
  sample.use <- 1:n
  if(n>200){
    sample.use <- sample(1:n,200,FALSE)
  }
  cors.screen <-
    apply(
      inter.schedule.Xint,
      1,
      FUN = function(x,
                     X.interf.mat.0 = X.interf.mat,
                     X.interf.outer.0 = X.interf.outer,
                     resy.0 = resy,
                     replaceme.0 = replaceme,
                     sample.use.0=sample.use) {
        z1 <- scale2(X.interf.mat.0[sample.use.0, x[1]])
        z2 <- scale2(X.interf.outer.0[sample.use.0, x[2]])
        create.inter.basis(theta = rot.est, resy.0[sample.use.0], z1, z2, replaceme.0[sample.use.0])$logcor
      }
    )
  ind2 <- sort(cors.screen, dec = T, index = T)$ix[1:10]
  X.out <- Xmat.out <- Xbasis <- NULL
  for (i.ind in 1:10) {
    x.temp <-
      scale2(X.interf.mat[, inter.schedule.Xint[ind2[i.ind], 1]])
    x.temp2 <-
      scale2(X.interf.outer[, inter.schedule.Xint[ind2[i.ind], 2]])
    x.temp3 <-
      create.inter.basis(theta = rot.est, resy, x.temp, x.temp2, replaceme)$basis
    X.out <- cbind(X.out, x.temp)
    Xmat.out <- cbind(Xmat.out, x.temp2)
    Xbasis <- cbind(Xbasis, x.temp3)
  }
  
  keeps.basis <- check.cor(X.out, thresh = 0.00001)$k
  Xbasis <- Xbasis[, keeps.basis]
  X.out <- X.out[, keeps.basis]
  Xmat.out <- Xmat.out[, keeps.basis]
  create.temp <-
    function(x)
      create.inter.basis(rot.est, resy, res1.2, x, replaceme)$basis
  X.basis.temp <- Xbasis
  ##If no covariates, return this
  if (ncol(X.basis.temp) == 0)
    return(NULL)
  which.keep <- NULL
  for (i.keep in 1:ncol(X.basis.temp))
    which.keep[i.keep] <-
    max(abs(cor(Xbasis[, i.keep], X.basis.temp))) == 1
  Xvars <- X.out[, which.keep]
  Xmatvars <- Xmat.out[, which.keep]
  output <-
    list(
      "Xvars" = as.matrix(Xvars),
      "Xmatvars" = as.matrix(Xmatvars),
      "Xbasis" = X.basis.temp
    )
  output <- lapply(output, cleanNAs)
  return(output)
}



create.inter.basis <- function(theta, resy, res1.2, x, replaceme) {
  ## Create interf kernel off basis using replaceme==1
  matvar <- x
  resvar <- res1.2#as.vector(scale2((res1.2)))
  yvar <- resy#as.vector(scale2(resy))
  basis <-
    basisvec_rcpp(
      matvec = matvar,
      resvec = resvar * (replaceme > 4),
      onesvec = 1 * (replaceme > 4),
      theta
    )$basis
  basis <- as.vector(basis)
  cor1 <- ifelse(sd_cpp2(basis) > 0,
                 abs(corSpearman(basis[replaceme > 4], resy[replaceme > 4])),
                 #abs(cor(basis[replaceme > 4], resy[replaceme > 4], method = "spearman")),
                 1e-5)
  cor.out <- log(max(1e-5, cor1))
  output <- list("basis" = basis, "logcor" = cor.out)
  return(output)
}

makebases.interf <- function(X.interfy, resy, replaceme) {
  if (length(X.interfy$Xvars) == 0)
    return(NULL)
  if (ncol(X.interfy$Xvars) > 0) {
    cor1 <- function(theta, ...)
      create.inter.basis(theta, ...)$logcor
    X.resy <- 0 * X.interfy$Xmatvars
    
    replaceme.temp<-replaceme
    if(sum(replaceme>4)>1000){
      which.gt4<-which(replaceme.temp>4)
      replaceme.temp[replaceme.temp>4]<-0
      replaceme.temp[sample(which.gt4,1000,FALSE)]<-5
    }
    
    if(length(replaceme)>3000){
      which.zero<-sample(1:length(replaceme),3000,FALSE)
      replaceme.temp[which.zero]<-0
    }
    
    
    rot.est <- log(1.06 * (sum(replaceme.temp>4)) ^ -.2 * (1 / 3) ^ -.2)
    
    for (i.interf in 1:ncol(X.interfy$Xvars)) {
      opt.y <-
        optimize(
          cor1,
          interval = c(-2, 2)+rot.est,
          resy = resy[replaceme.temp!=0],
          res1.2 = X.interfy$Xvars[replaceme.temp!=0, i.interf],
          x = X.interfy$Xmatvars[replaceme.temp!=0 , i.interf],
          replaceme.temp[replaceme.temp!=0],
          maximum = T,
          tol = 0.001
        )
      # opt.y <- NULL
      # opt.y$maximum <- rot.est
      # create.inter.basis<-function(theta,resy,rest,x,replaceme)
      resy.s <- scale2(resy)
      X.interfy$Xvars <- apply(X.interfy$Xvars, 2, scale2)
      X.interfy$Xmatvars <- apply(X.interfy$Xmatvars, 2, scale2)
      
      X.resy[, i.interf] <-
        create.inter.basis(
          theta = opt.y$maximum,
          resy = resy.s,
          res1.2 = X.interfy$Xvars[, i.interf],
          x = X.interfy$Xmatvars[, i.interf],
          replaceme.temp
        )$basis
    }
  }
  X.resy <-
    cbind(lm(resy ~ X.resy, weights = 1 * (replaceme %in% c(3, 4)))$fit, X.resy)
  orthog.me(resy, X.resy, weights.lm = 1 * (replaceme %in% c(3, 4)))
}


scale2 <- function(z) {
  sdz <- sd_cpp2(z)
  if (sdz == 0) {
    out <- rep(1, length(z))
  } else{
    out <- (z-mean(z))/sdz
    #out <- as.vector(scale(z))
  }
  out
}

bSpline2 <- function(x, ...) {
  b1 <- bSpline(x, ...)
  b2 <- bSpline(-x, ...)
  cbind(b2[, ncol(b2)], b1)
}


bs.me<-function(x,degree=5){
  n<-length(x)
  x<-(rank(x)-1)/(n)
  x<-2*x-1
  basis.out<-NULL
  for(i.basis in 1:degree){
   # basis.out<-cbind(basis.out,cos(i.basis*acos(x)),cos(-i.basis*acos(x)))
  }
  basis.out<-cbind(x,bs(x,deg=3,knots=median(x)),bs(-x,deg=3,knots=median(x)))
  basis.out[,check.cor(basis.out,0.0001)$k]
  
}

bs.me_old <- function(x, degree = 3) {
  if (length(unique(x)) <= 2)
    return(cbind(x))
  if (length(unique(x)) <= 3)
    return(cbind(x, abs(x - mean(x))))

  basis.out <-
    cbind(x, bSpline2(x, deg = 3), bSpline2(x, deg = 4)[, -c(1, 5)])#, bSpline2(x, deg = 5), bSpline2(x, deg = 6))
  basis.out <- basis.out[, apply(basis.out, 2, sd_cpp2) > 0.00001]
  basis.out
}



## Calculate loo error
loo.est <- function(y0, X) {
  X <- as.matrix(X)
  y0 <- as.vector(y0)
  n <- length(y0)
  k <- ncol(X)
  lm.y <- lm(y0 ~ X)
  res2 <- lm.y$res / (1 - influence(lm.y)$hat)
  rss <- mean(res2 ^ 2)
  #rss<-mean(lm.y$residuals^2)/(1-k/n)^2
  rss
}


bSpline2 <- function(x, ...) {
  b1 <- bSpline(x, ...)
  b2 <- bSpline(-x, ...)
  cbind(b2[, ncol(b2)], b1)
}


bs.me_new <- function(x) {
  #sd.x<-sd_cpp2(x)
  #x<-x/sd.x
  x <- (x - min(x)) / abs(diff(range(x)))
  m2 <-
    cbind(x,
          bSpline2(x, df = 3),
          bSpline2(x, df = 4),
          bSpline2(x, df = 5),
          bSpline2(x, df = 6))
  return(m2)
  
}


cleanNAs <- function(X){
  return(X)
 if(length(X)==0) return(X)
  X<-as.matrix(X)
 if(ncol(X) < 2) return(as.matrix(X))
   return(as.matrix(X[, colSums(is.na(X)) == 0]))
}


makebin <- function(k) {
  m1 <- matrix(NA, nr = 2 ^ k, nc = k)
  runvec <- 1:(2 ^ k)
  for (i in 1:ncol(m1)) {
    m1[, i] <- runvec %% (2)
    runvec <- ceiling(runvec / 2)
  }
  (m1[rowSums(m1) > 0,] == 1)
}

#Check for positivity!


pos_measure <- function(obj, trim = .0125) {
  if (class(obj) == "list")
    obj <- obj$treat.res
  #obj2<-t(apply(obj,1,FUN=function(x) (x-mean(x,na.rm=T))))
  
  #d1<-density(obj[!is.na(obj)])
  #mode.dens<-d1$x[d1$y==max(d1$y)][1]
  #obj<-obj-mode.dens
  obj <- obj - mean(obj, na.rm = T)
  obj <- obj / sd_cpp2(as.vector(obj))
  
  vec1 <- rowMeans(obj ^ 4, na.rm = T)#-rowMeans(obj^3,na.rm=T)
  vec2 <- rowMeans(obj ^ 2, na.rm = T) ^ 2
  vec3 <- rowMeans(obj ^ 3, na.rm = T)
  #vec4<-rowMeans(obj2^2,na.rm=T)^2
  
  n.v <- length(vec1)
  vec1 <- (vec1 / (vec2) - 3)#/(1+vec3/vec2^1.5)
  ks <- median(rowSums(!is.na(obj)))
  vec1 <- as.vector(vec1)
  #vec1<-vec1#/var(as.vector(vec1),na.rm=T)^2
  #qqnorm(vec1)
  #qqnorm(vec1^bc1$lambda)
  n.v <- length(vec1)
  rs.expected1 <-
    sort(sapply(
      1:n.v,
      FUN = function(z, k0 = ks) {
        d1 <-
          rnorm(k0)
        (mean(d1 ^ 4) / mean(d1 ^ 2) ^ 2 - 3)
      }
    ))#/(1+mean(d1^3)/mean(d1^2)^1.5)} )
  rs.expected2 <-
    sort(sapply(
      1:n.v,
      FUN = function(z, k0 = ks) {
        d1 <-
          rnorm(k0)
        (mean(d1 ^ 4) / mean(d1 ^ 2) ^ 2 - 3)
      }
    ))#/(1+mean(d1^3)/mean(d1^2)^1.5)} )
  rs.expected3 <-
    sort(sapply(
      1:n.v,
      FUN = function(z, k0 = ks) {
        d1 <-
          rnorm(k0)
        (mean(d1 ^ 4) / mean(d1 ^ 2) ^ 2 - 3)
      }
    ))#/(1+mean(d1^3)/mean(d1^2)^1.5)} )
  rs.expected <- (rs.expected1 + rs.expected2 + rs.expected3) / 3
  #expected<-sort(qchisq((1:n.v)/(n.v+1),df=1))^2
  expected <- sort(rs.expected)
  #bc2<-EnvStats::boxcox(expected,optimize=T)
  #bc1<-EnvStats::boxcox(vec1,optimize=T)
  vec1 <- vec1#^(bc1$lambda)
  expected <- expected#^(bc2$lambda)
  ind <- sort(vec1, dec = F, ind = T)
  obs <- ind$x
  
  #obs<-obs-mean(obs)
  #expected<-expected-mean(expected)
  
  diff <- obs - expected
  trim.num <- floor(trim * n.v)
  diff <- diff[-c(1:trim.num)]
  diff <- rev(rev(diff)[-c(1:trim.num)])
  #diff<-diff-mean(diff)
  
  output <-
    list(
      "observed" = obs,
      "expected" = expected,
      "diff" = diff,
      "index" = ind
    )
}

make.rank <- function(z) {
  (0:(length(z) - 1)) / (length(z) - 1) * 100
}

check.cor <- function(X, thresh = 0, nruns = 2) {
  # set thresh unless more than 1 value is provided
  # settings is to set 0.001 for Xs with p > 100, stricter limit if p is smaller
  
  X <- as.matrix(X)
  X<-cleanNAs(X)
  keeps.out <- rep(FALSE, ncol(X))
  dropints <- apply(X, 2, sd_cpp2) == 0
  X <- X[, !dropints]
  
  if (length(thresh) == 0)
    thresh <- ifelse(ncol(X) > 100, 0.001, 0.0001)
  
  # if of only column-1, nothing to check
  if (NCOL(X) == 1)
    return(list("keeps" = TRUE))
  
  # track across columns
  drops.run <- rep(0, ncol(X))
  
  # for a given number of runs (defaults to 2)
  for (j in 1:nruns) {
    # container for each of p-columns
    drops <- rep(0, ncol(X))
    
    ## generate one n-vector of stdnorm error then see how it correlates with
    ## each column of X.
    cor1 <- as.vector(abs(cor(X, rnorm(nrow(
      X
    )))))
    
    
    for (i in 1:(ncol(X) - 1)) {
      ind_subseq_cols <- (i + 1):ncol(X)
      
      # the noise-correlation mags for subsequent columns
      cor.temp <- cor1[ind_subseq_cols]
      
      ## mark as (-1) any column whose noise correlation is indistinguishable,
      ## as defined from thresh, with the column i.
      drops[ind_subseq_cols][abs(cor.temp - cor1[i]) < thresh] <-
        (-1)
    }
    
    # sum all -1s together j times
    drops.run <- drops + drops.run
  }
  
  ## if column was marked as too similar with some other colum for all nruns-runs,
  ## drop it. else keep it.
  keeps <- drops.run != (-nruns)
  keeps.out[!dropints][keeps] <- TRUE
  out1 <- list("keeps" = keeps.out)
  return(out1)
}


### GCV for sparsereg

sparsereg_GCV <- function(y0,X0,id0=NULL, usecpp=TRUE){
  n<- length(y0)
  p <- ncol(X0)
  if(ncol(as.matrix(X0))<3 ){
    model.out <- sparsereg(y=y0,X=X0, id=id0 )
    return(model.out)
   }
  if(usecpp==FALSE){
  alpha.max <- max(n*log(ncol(X0)), ncol(X0)*1.25)
  alpha.min <- min(ncol(X0)*1.25, n*log(ncol(X0))/2)
  alpha.schedule <- seq((alpha.max), alpha.min,length=10)
  GCV.out<-sapply(alpha.schedule, FUN=function(z) 
    tryCatch(sparsereg(y=y0,X=X0, id=id0, edf=T, alpha.prior = "custom", alpha.use = z )$GCV,
             error = function(h) 1e7, warning = function(h) 1e7
    )
  )# Close out sapply
  # a.min <-alpha.schedule[GCV.out==min(GCV.out)]
  lm1<-lm(GCV.out~alpha.schedule+I(alpha.schedule^2))
  alpha.lm <- -(lm1$coef[2])/(2*lm1$coef[3])
  alpha.min <- ifelse(
    lm1$coef[3] > 0 ,
    -(lm1$coef[2])/(2*lm1$coef[3]),
    alpha.schedule[which.min(GCV.out)]
  )
  alpha.min <- rcppClamp(alpha.min,alpha.schedule[1],rev(alpha.schedule)[1])
  model.out <- sparsereg(y=y0,X=X0, id=id0, edf=T, alpha.prior = "custom", alpha.use = alpha.min )
  
  } else{
  
  ## New using Cpp bayeslasso to find alpha!
  alpha.schedule <-
    seq(log(8*n*log(p)),log(p),length=8)
  alpha.schedule <-
    seq(log(8*n*log(p)),log(p),length=5)
  #alpha.max <- max(n*log(ncol(X0)), ncol(X0)*1.25)
  #alpha.min <- min(ncol(X0)*1.25, n*log(ncol(X0))/2)
  #alpha.schedule <- log(seq((alpha.max), alpha.min,length=10))
  X0.temp <- apply(X0,2,scale2)
  gcv.out <- log(sapply(alpha.schedule, FUN=function(z) bayesLasso(y0,cbind(1,X0.temp),exp(z))$GCV ))
  lm1<-lm(gcv.out~alpha.schedule+I(alpha.schedule^2))
  alpha.lm <- -(lm1$coef[2])/(2*lm1$coef[3])
  alpha.min <- ifelse(
    lm1$coef[3] > 0,
    -(lm1$coef[2])/(2*lm1$coef[3]),
    alpha.schedule[which.min(gcv.out)]
  )
  alpha.min <- rcppClamp(alpha.min,alpha.schedule[1],rev(alpha.schedule)[1])
  #print(alpha.min)
  #model.out <- sparsereg(y=y0,X=X0, id=id0, edf=T, alpha.prior = "custom", alpha.use =exp( alpha.min ))
  model.out <- sparsereg(y=y0,X=X0, id=id0, edf=T, alpha.prior = "custom", alpha.use =exp( alpha.min ))
  
  } # Closes out if
  model.out
}




make.soe <- function(basesy.all.inner,id,sy,y,fits.y,fits.REs.y,replaceme){
  
  if(length(sy$chol.inv)<=1) return(NULL)
  basesy.all.temp <-cbind(1,basesy.all.inner)
  basesy.all.temp <- cbind(basesy.all.temp[,sy$beta!=0])
  
  mat.all<-basesy.all.temp%*%(sy$chol.inv)
  
  mat.control<-cbind(1,fits.REs.y,fits.y)
  weights.est <- 1*(replaceme <3)
  y2<-fastLm(y~mat.control,weights=weights.est)$res
  mat.all2<-apply(mat.all,2,FUN=function(z) scale2(fastLm(z~mat.control,weights=weights.est)$res))
  
  mat.orthog<-as.matrix(orthog.me(y2,mat.all2,weights.lm=weights.est))
  if(ncol(mat.orthog)==1) return(mat.orthog[weights.est==0,])
  colnames(mat.orthog)<-paste("X",1:ncol(mat.orthog),sep="_")
  
  svd.temp <- svd(mat.orthog[weights.est==1,])
  mat.orthog2 <- mat.orthog%*%t(svd.temp$v)
  #print(svd.temp$d/svd.temp$d[1])
  keeps <- (cumsum(svd.temp$d^2)/sum(svd.temp$d^2)) < 1
  #print(c(ncol(basesy.all.inner),sum(keeps)))
  if(sum(keeps)>.75*sum(replaceme<3)) {
    maxcol <- ceiling(.75*sum(replaceme<3))
    keeps[-c(1:maxcol)] <- FALSE
  }
  mat.orthog2[weights.est==0,keeps]
}      

#fit.interference = TRUE, 
#fit.treatment.heteroskedasticity = TRUE

allbases <- function(y,y2.b,treat,treat.b,treat.y,X,id, replaceme, 
                     fit.interference, 
                     fit.treatment.heteroskedasticity,
                     inter.schedule.obj = NULL
                     ){
  
  n<-length(y)
  replaceme <- make.replaceme(y, id)

  #generate.bases <- function(y2.b, y, basesy0, X, id=NULL, replaceme) {
  ## Make bases for outcome, treatment----
  #cat("#####  Step 1: Constructing conditional mean bases\n")

  
  basest0 <- generate.bases(treat.b, treat, X, X, id=NULL,replaceme)
  basest0 <- apply(basest0,2,scale2)
  treat.b.2 <-
    treat.b - basest0%*%sparsereg_GCV(treat.b[replaceme>4], basest0[replaceme>4,])$coef#, EM = T, verbose = F, id=NULL)$fitted
  #basest0 <- cbind(basest0, generate.bases(treat.b.2, treat, X, X, id=NULL,replaceme))
  basest0 <- cleanNAs(basest0)
  
  basesy0 <- generate.bases(y2.b, y2.b, cbind(X), cbind(X),id=NULL, replaceme)
  y2.b.2 <- y2.b - sparsereg(y2.b, basesy0, EM = T, verbose = F, id=NULL)$fitted
  #basesy0.2 <- generate.bases(y2.b.2, y2.b.2, cbind(X), cbind(X),id=NULL, replaceme) ## Don't need REs--already partialed out!
  #basesy0<-cbind(basesy0,basesy0.2)
  basesy0 <- cleanNAs(basesy0)
  
  ## Strip out linear dependencies--it just helps
  keeps.t <- which(!is.na(lm(treat ~ basest0)$coef[-1]))
  keeps.y <- which(!is.na(lm(y ~ treat + basesy0)$coef[-c(1:2)]))
  
  basest0 <- as.matrix(basest0[, keeps.t])
  basesy0 <- as.matrix(basesy0[, keeps.y])
  # 
  
  baset0 <- basest0[,check.cor(basest0)$k]
  basey0 <- basesy0[,check.cor(basesy0)$k]
  
  basest0 <- orthog.me(treat, basest0, weights.lm = 1*(replaceme>4))
  basesy0 <- orthog.me(y, basesy0, weights.lm = 1*(replaceme>4))
  
  
  ## Initial outcome and propensity model ----
  sy <- sparsereg(y,
                  cbind(basesy0, basest0),
                  EM = T,
                  verbose = F, 
                  id=NULL)
  res1.y <- as.vector(scale2(lm(y ~ sy$fit)$res))
  
  st <- sparsereg(treat, basest0, EM = T, verbose = F, id=NULL)
  res1 <- as.vector(scale2(treat - st$fit))
  
  ## Make residual nonparametric terms ----
  res2.b <- (res1)
 # generate.bases <- function(y2.b, y, basesy0, X, id=NULL, replaceme) {
  basest2.0 <- NULL    
  res1.2 <- scale2(res1)
  
  if(fit.treatment.heteroskedasticity){
  basest2.0 <- generate.bases(res2.b, res2.b, X, X, id=NULL, replaceme,alwaysinter = res2.b)
  basest2.0 <- cleanNAs(basest2.0)
   
   st.2 <-
     sparsereg(res1, apply(
       basest2.0,
       2,
       FUN = function(z)
         scale2(z) * scale2(res1)
     ), id=NULL)
   res1.2 <- scale2(res1 - st.2$fit)
   
   res2.b.2 <- (res1.2)
   
    basest2.0 <- cbind(basest2.0, generate.bases(res2.b.2, res2.b.2, X,X, id=NULL,replaceme,alwaysinter = res2.b.2))
    basest2.0 <- cleanNAs(basest2.0)
   basest2.0 <- orthog.me(res1.2, basest2.0, weights.lm = 1*(replaceme>4))
}
  
  ## Split sample ----
  replaceme <- replaceme0 <- rep(c(1:6), times = floor(n / 2))[1:(n + 4)]
  replaceme <- make.replaceme(y,id)#sample(replaceme[1:n])
  
  ## Make interference bases ----
  X.interfy <- X.interft <- NULL
  if(fit.interference){
  X.interfy <- generate.Xinterf(res1.y, X, treat, replaceme)
  X.interft <- generate.Xinterf(res1.2, X, NULL, replaceme)
  }
  
  ## Generate splits ----
  # fits.replaceme <- lm(lm(res1~st.2$fit)$res~X.interft$Xbasis)$fit#fastLm(treat ~ st$fit)$res
  # fits.replaceme.y <- fastLm(y ~ sy$fit)$res
  
  if(length(basesy0)==0) basesy0<-cbind(rnorm(n),rnorm(n))
  if(length(basest0)==0) basest0<-cbind(rnorm(n),rnorm(n))
  if(length(basest2.0)==0) basest2.0<-cbind(rnorm(n),rnorm(n))
  

  #if(abs(cor(fits.replaceme.y,y))>abs(cor(fits.replaceme,treat))) fits.replaceme<-fits.replaceme.y
  output<-list("basesy0"=basesy0,
               "basest0"=basest0,
               "basest2.0"=basest2.0,
               "X.interfy"=X.interfy,
               "X.interft"=X.interft,
               "replaceme" = replaceme
  )
}

sd_cpp2 <- function(z) sd(z)
