#' Estimating a Marginal Effect in the Partially Linear Causal Effect Model.
#'
#' @param y Vector of the dependent, or outcome variable.
#' @param treat Treatment variable, for which the marginal effect is desired.
#' May be continous, binary, or count.
#' @param X Matrix of covariates to be adjusted for in estimating the marginal effect. 
#' A machine learning method is used to adjust not just for the covariates but also
#' nonlinear, interactive functions of the covariates.
#' @param id A vector with elements denoting groupings to fit via random effects.
#' @param num.fit The number of cross-fitting iterations to run.  Default is 25.
#' @param var.type The type of variance estimate, passed to \code{vcovHC} in sandwich. Options are
#' \code{HC0}, \code{HC1}, \code{HC2}, \code{HC3}, with \code{HC3} the default.
#' @param sens Whether to fit the sensitivity analysis using \code{sensemakr}.
#' @param printevery How frequently to print output during estimation.  Default is 10.
#' @param fit.interference Whether to adjust for interference.  Default is \code{TRUE}.
#' @param fit.treatment.heteroskedasticity Whether to adjust for bias induced through
#' heteroskedasticity in the treatment variable.  Default is \code{TRUE}.
#'
#' @export
#'
#'
#'
#' @return \describe{
#' \item{point}{The estimated marginal effect.}
#' \item{se}{Standard error of the estimate of the marginal effect.}
#' \item{sens}{Results from the sensitivity analysis.}
#' \item{treat.res}{Residuals from modeling the treatment model. Used in the
#' positivity diagnostic.}
#'}
#'
#' @references Ratkovic, Marc.  2021. "Relaxing Assumptions, Improving Inference:
#' Utilizing Machine Learning for Valid Causal Inference."  Working Paper.
#'
#' @rdname plce
#'
#' @importFrom splines bs
#' @importFrom sandwich vcovHC
#' @importFrom sensemakr sensemakr
#' @importFrom splines2 bSpline
#' @importFrom RcppEigen fastLm
#' @importFrom ccaPP corSpearman
#' @importFrom ccaPP corPearson
#' @importFrom lme4 lmer
#' @examples
#'  \dontrun{
#' ## This example takes you through an implementation
#' # of the PLCE model. We first create the sample size and
#' # number of covariates.
#' set.seed(1234)
#' n <- 1000
#' p <- 5
#' 
#' ## Generate covariate matrix
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' X<- apply(X,2,scale)
#' 
#' ## Generate random effects
#' ids.map <- sample(letters[1:10], n, TRUE)
#' res.map <- rnorm(10)
#' names(res.map) <- letters[1:10]
#' res.true <-  (res.map[ids.map])
#' 
#' ## Generate the treatment and outcome
#' treat <- (X[, 1]) +  res.true + rnorm(n)
#' Y <- 1 + treat * (X[, 1] ^ 2) + res.true + rnorm(n)
#' 
#' ## Fit the PLCE model
#' plce1 <-
#' plce(
#'  y=Y,
#'  treat=treat,
#'  X=X, id = ids.map,
#'  printevery = 1,
#'  num.fit = 5
#'   )
#' 
#' ## Results: Point estimate, standard error,
#' ##   sensitivity analysis
#' plce1$point
#' plce1$se
#' plce1$sens
#'  }



plce <-
  function(y,
           treat,
           X,
           id = NULL,
           num.fit = 25,
           var.type = "HC3",
           sens = TRUE,
           printevery = 10,
           fit.interference = TRUE,
           fit.treatment.heteroskedasticity = TRUE) {
    n <- length(y)
    k <- ncol(X)
    y0 <- y
    treat0 <- treat
    X0 <- X
    
    REs <- F
    if (length(id) > 0) {
      #X.REs <- model.matrix(~id-1)
      REs <- T
    }
    ## Format inputs  -----
    treat <- as.vector(scale2(treat0))
    y <- as.vector(scale2(lm(y0 ~ treat)$res))
    if (length(id) > 0) {
      treat <- suppressMessages(scale2(treat - predict(lmer(treat ~ (
        1 | id
      )))))
      y <- suppressMessages(scale2(y - predict(lmer(y ~ (
        1 | id
      )))))
      
    }
    
    X <- cleanX(X)
    k <- ncol(X)
    
    if (REs) {
      y <- suppressMessages(scale2(y - predict(lmer(y ~ (
        1 | id
      )))))
      treat <-
        suppressMessages(scale2(treat - predict(lmer(treat ~ (
          1 | id
        )))))
      
    }
    

    
    y2.b <- y
    treat.y <- treat.b <- treat
    
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
    treat.res.run <- matrix(NA, nrow = n, ncol = (num.fit * 2))
    
    beta.t.init <- beta.y.init <- NULL
    # if(length(id)==0){
    #   beta.t.init <- solve_cpp(crossprod(cbind(1,basest0))+sd_cpp2(treat)^2*diag(ncol(basest0)+1),crossprod(cbind(1,basest0),treat) )
    #   beta.y.init <- solve_cpp(crossprod(cbind(1,basesy0))+sd_cpp2(y)^2*diag(ncol(basesy0)+1),crossprod(cbind(1,basesy0),y0) )
    # }
    cat("#####  Beginning cross-fitting for the marginal effect\n")
    
    for (i.hoe in 1:num.fit) {
      ## Create matrices for fitting ----
      replaceme <- make.replaceme(y, id)
     allbases.obj <-
        allbases(
          y,
          y2.b,
          treat,
          treat.b,
          treat.y,
          X,
          id,
          replaceme,
          fit.interference,
          fit.treatment.heteroskedasticity
        )
      
      basest <- basest0 <- allbases.obj$basest0
      basesy <- basesy0 <- allbases.obj$basesy0
      basest2 <- basest2.0 <- allbases.obj$basest2.0
      X.interfy <- allbases.obj$X.interfy
      X.interft <- allbases.obj$X.interft
      replaceme <- allbases.obj$replaceme
      
      
      
      id.temp <- NULL
      #if(length(id)>0) id.temp <- id[replaceme>4]
      ## Propensity model for interference, heteroskedasticity estimated off splits 5-6 ----
      colnames(basest) <- paste("X", 1:ncol(basest), sep = "_")
      st <-
        sparsereg_GCV(treat[replaceme > 4], X0 = basest[replaceme > 4,], id0 = id.temp)
      fits.all <-
        cbind(basest[, intersect(colnames(basest), names(st$coef))]) %*% st$coef[intersect(colnames(basest), names(st$coef))]
      if (length(id.temp) > 0) {
        fits.all <-
          fits.all + as.vector(st$REs[id])#X.REs[,names(st$REs)]%*%st$REs
      }
      fits <- fits.all[replaceme < 5]
      # m1<-mget(ls())
      # save(m1,file="diagnose.Rda")
      res <- as.vector(scale2(treat[replaceme < 5] - fits))
      
      ## Outcome model for interference estimated off splits 5-6 ----
      colnames(basesy) <- paste("X", 1:ncol(basesy), sep = "_")
      sy <-
        sparsereg_GCV(y[replaceme > 4], X0 = basesy[replaceme > 4,], id0 = id.temp)
      names.inter.y <- intersect(colnames(basesy), names(sy$coef))
      
      fits.y <-
        scale2(basesy[, names.inter.y] %*% sy$coef[names.inter.y])
      if (length(id.temp) > 0) {
        fits.y <- fits.y + as.vector(sy$REs[id])#X.REs%*%sy$REs
      }
      
      resy <- scale2(y - fits.y)
      
      
      
      ## Construct interference bases off this split----
      ## First, the residuals we need
      res1 <- scale2(treat - fits.all)
      colnames(basest2.0) <-
        paste("X", 1:ncol(basest2.0), sep = "_")
      
      ## No need for REs--already taken out.
      st.2 <-
        sparsereg_GCV(res1[replaceme > 4], apply(
          basest2.0[replaceme > 4,],
          2,
          FUN = function(z)
            scale2(z * scale2(res1[replaceme > 4]))
        ))
      fits.st.2 <- basest2.0 %*% st.2$coef
      res1.2 <- scale2(res1 - fits.st.2)
      
      ## Make interference bases ----
      if (fit.interference) {
        X.resy <- cleanNAs(makebases.interf(X.interfy, resy, replaceme))
        X.rest <-
          cleanNAs(makebases.interf(X.interft, res1.2, replaceme))
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
      id.temp <- NULL
      if (length(id) > 0)
        id.temp <- id[replaceme < 5]
      
      # m1<-mget(ls())
      # save(m1,file="diagnose.Rda")
      
      ## Bases to be passed to selection and estimation set ----
      if (fit.treatment.heteroskedasticity) {
        bases.t.temp <-
          cbind(basest[replaceme < 5,], res * basest2[replaceme < 5,])
      } else{
        bases.t.temp <- basest[replaceme < 5,]
      }
      bases.y.temp <- basesy0[replaceme < 5,]
      bases.t.temp <- cleanNAs(bases.t.temp)
      bases.y.temp <- cleanNAs(bases.y.temp)
      
      
      if (length(X.rest) > 0)
        X.rest.temp <-
        as.matrix(X.rest)[replaceme < 5,]
      else
        X.rest.temp <- NULL
      if (length(X.resy) > 0)
        X.resy.temp <-
        as.matrix(X.resy)[replaceme < 5,]
      else
        X.resy.temp <- NULL
      
      
      h1 <-
        hoe.inner(
          y = y.temp,
          treat = treat.temp,
          y0 = y0.temp,
          treat0 = treat0.temp,
          basest = bases.t.temp,
          basesy = bases.y.temp,
          replaceme = replaceme.temp,
          var.type,
          lin.adj = NULL,
          lin.adj.t = NULL,
          sens = sens,
          X.rest.temp,
          X.resy.temp,
          beta.t.init,
          beta.y.init,
          id = id.temp,
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
    
    sens.out <- rowMeans(sens.run)
    sens.out[1] <- point.est
    sens.out[2] <- (var.est) ^ .5
    sens.out[3] <- sens.out[1] / sens.out[2]
    out1 <-
      list(
        "point" = point.est,
        "se" = (var.est) ^ .5,
        "treatpoint" = points.run,
        "sens" = sens.out,
        "treat.res" = treat.res.run
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
           replaceme,
           var.type,
           lin.adj = NULL,
           lin.adj.t = NULL,
           sens = NULL,
           X.rest,
           X.resy,
           beta.t.init,
           beta.y.init,
           id,
           fit.interference,
           fit.treatment.heteroskedasticity) {
    ## Save original inputs ----
    n <- length(y)
    
    ##Declare holders for loop
    sens.run <- point.out <- se.out <- NULL
    treat.res.inner <- matrix(NA, nrow = n, ncol = 2)
    
    if (fit.interference) {
      basest.all <- cbind(basest, X.rest)
      basesy.all <- cbind(basesy, X.resy)
    } else{
      basest.all <- cbind(basest)
      basesy.all <- cbind(basesy)
    }
    colnames(basest.all) <-
      paste("X", 1:ncol(basest.all), sep = "_")
    colnames(basesy.all) <-
      paste("X", 1:ncol(basesy.all), sep = "_")
    
    basest.all <- basest.all[, check.cor(basest.all, thresh = 0.001)$k]
    basesy.all <- basesy.all[, check.cor(basesy.all, thresh = 0.001)$k]
    
    basest.all <- apply(basest.all, 2, scale2)
    basesy.all <- apply(basesy.all, 2, scale2)
    
    y.orthog <- y
    treat.orthog <- treat
    if (length(id) > 0) {
      X.dum <- cbind(sample(y), sample(y))
      colnames(X.dum) <- c("X1", "X2")
      X.dum <- apply(
        X.dum,
        2,
        FUN = function(z)
          fastLm(z ~ cbind(1, y))$res
      )
      y.orthog <- y - sparsereg_GCV(y, X.dum, id0 = NULL)$fit
      
      colnames(X.dum) <- c("X1", "X2")
      X.dum <-
        apply(
          X.dum,
          2,
          FUN = function(z)
            fastLm(z ~ cbind(1, treat))$res
        )
      treat.orthog <- treat - sparsereg_GCV(treat, X.dum, id0 = NULL)$fit
      
    }
    
    
    for (i.outer in 1:2) {
      ## Select bases off splits 1,2 ----
      basest.all.inner <-
        orthog.me(treat.orthog, basest.all, weights.lm = (replaceme < 3))
      basesy.all.inner <-
        orthog.me(y.orthog, basesy.all, weights.lm = (replaceme < 3))
      
      #basest.all.inner <-svd(apply(basest.all,2,scale2))$u
      basest.all.inner <-
        basest.all[, check.cor(basest.all[replaceme < 3, ])$k] #orthog.me(treat.orthog, basest.all, weights.lm = (replaceme<3))
      basesy.all.inner <-
        basesy.all[, check.cor(basesy.all[replaceme < 3, ])$k] #orthog.me(y.orthog, basesy.all, weights.lm = (replaceme<3))
      
      id.temp <- NULL
      if (length(id) > 0)
        id.temp <- id[replaceme < 3]
      
      st <-
        sparsereg(treat[replaceme < 3],
                  basest.all.inner[replaceme < 3,],
                  id.temp,
                  edf = T,
                  alpha.prior = "parametric")
      select.lin.t <- (abs(st$coef) > 1e-3)
      
      sy <-
        sparsereg(y[replaceme < 3], basesy.all.inner[replaceme < 3,], id.temp, edf =
                    T)
      select.lin.y <- (abs(sy$coef) > 1e-3)
      
      soe.t <- soe.y <- NULL
      
      threshind <-
        floor(min(n / 2 * .8, length(st$coef) + length(sy$coef) - 2))
      
      thresh1 <-
        sort(c(abs(st$coef), abs(sy$coef)), decreasing = T)[threshind]
      select.lin.y[abs(sy$coef) < thresh1] <- F
      select.lin.t[abs(st$coef) < thresh1] <- F
      
      fits.t <- cbind(basest.all.inner) %*% st$coef
      fits.y <- cbind(basesy.all.inner) %*% sy$coef
      
      fits.REs.t <- fits.REs.y <- NULL
      
      
      soe.t <-
        make.soe(apply(basest.all.inner, 2, scale2),
                 id,
                 st,
                 treat,
                 fits.t,
                 fits.REs.t,
                 replaceme)
      #y.temp <-lm(y~soe.t,weights=1*(replaceme < 3) )$res
      soe.y <-
        make.soe(apply(basesy.all.inner, 2, scale2),
                 id,
                 sy,
                 y,
                 fits.y,
                 fits.REs.y,
                 replaceme)
      
      if (TRUE) {
        Xmat.adjust <- as.matrix(cbind(# basest.all.inner[replaceme > 2,],
          # basesy.all.inner[replaceme > 2,],
          fits.t[replaceme > 2],
          fits.y[replaceme > 2],
          soe.t,
          soe.y,
          NULL))
        
        
        if (length(Xmat.adjust) == 0) {
          Xmat.adjust <- matrix(rnorm(sum(replaceme > 2)))
          Xmat.adjust <-
            matrix(lm(Xmat.adjust ~ y[replaceme > 2] + treat[replaceme > 2])$res)
        }
        
      } else {
        Xmat.adjust <- as.matrix(cbind(lin.adj[replaceme > 2,]))
        
      }
      
      if (F) {
        cat("#########################\n")
        cat("Dimension Check: Candidate Bases\n")
        cat(ncol(basesy.all.inner) + ncol(basest.all.inner), "\n")
        cat("Dimension Check:: Selected Bases\n")
        cat(ncol(Xmat.adjust), "\n")
      }
      
      if (length(id) > 0) {
        id.2 <- id[replaceme > 2]
        y2 <- lm(y0[replaceme > 2] ~ Xmat.adjust)$res
        treat2 <- lm(treat0[replaceme > 2] ~ Xmat.adjust)$res
        res.y.temp <- suppressMessages(predict(lmer(y2 ~ (1 | id.2))))
        res.t.temp <-
          suppressMessages(predict(lmer(treat2 ~ (1 | id.2))))
        
        
        Xmat.adjust <- cbind(Xmat.adjust, res.y.temp, res.t.temp)
        lm1 <-
          lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      } else {
        lm1 <-
          lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
      }
      lm1$coef[2]
      keeps <- which(!is.na(lm1$coef[-c(1:2)]))
      Xmat.adjust <- cleanNAs(Xmat.adjust[, keeps])
      lm1 <-
        lm(y0[replaceme > 2] ~ treat0[replaceme > 2] + Xmat.adjust)
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
        sens1 <- as.numeric(unlist(sens1)[-1])
        sens.run <- cbind(sens.run, sens1)
        rownames(sens.run) <-
          names(sensemakr(lm.s, name.use)$sensitivity_stats)[-1]
      }
      
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
        "sens" = (sens.run),
        "treat.res" = treat.res.inner,
        "loo" = NULL
      )
    return(output)
  }
