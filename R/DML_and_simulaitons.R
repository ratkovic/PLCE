


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
#'
#' @importFrom randomForest randomForest
#' @importFrom AER ivreg
#' @noRd

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
  function(outcome,
           predictor,
           method = "RF",
           train_idx,
           ...) {
    if (method == "RF") {
      # double machine learning
      args <- list(...)
      # train and predict the outcome with covariates X
      rfout <- randomForest(x = predictor[!train_idx, ],
                            y = outcome[!train_idx],
                            ntrees = args[['n_trees']])
      predict_diffed <-
        outcome[train_idx] - predict(rfout, newdata = predictor[train_idx, ])
      
    } else if (method == "LASSO") {
      #lsout <- cv.glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1, nfolds = 5)
      # lsout <- rlasso(outcome[!train_idx] ~ predictor[!train_idx, ])
      # lsout <- glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1)
      # predict_diffed <- outcome[train_idx] - predict(lsout, newdata = predictor[train_idx,])
    } else if (method == "SL") {
      #lsout <- cv.glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1, nfolds = 5
      # p1<-predictor[!train_idx, ]
      # p2<-predictor[train_idx, ]
      # lsout <- SuperLearner(outcome[!train_idx] ,data.frame(p1),newX=data.frame(p2),
      #                       SL.library =c("SL.glm", "SL.glmnet","SL.randomForest", "SL.loess",  "SL.polymars", "SL.mean"))
      # # lsout <- glmnet(x = predictor[!train_idx, ], y = outcome[!train_idx], alpha = 1)
      # predict_diffed <- outcome[train_idx] - lsout$SL.predict
    } else if (method == "KRLS") {
      # p1<-predictor[!train_idx, ]
      # p2<-predictor[train_idx, ]
      # lsout <- krls(y=outcome[!train_idx] ,X=(p1))
      # predict_diffed <- outcome[train_idx] - predict(lsout,newdata=p2)$fit
    } else {
      stop("Proper method required!\n")
    }
    
    return(predict_diffed)
  }




#' Generate simulation data
#' Generate data used in simulations.
#'
#'
#' @importFrom MASS mvrnorm
#' @noRd
#'
#'
#'
make.simdata <-
  function(n,
           k,
           meanhet = F,
           errhet = F,
           REs = T,
           inter = FALSE,
           rho = 0.5) {
    # n <- 500; k <- 5; meahet=F; errhet=F; REs=T; inter=FALSE; rho=0.5
    
    ##Format data
    
    ##Format data ----
    var.mat <- diag(k)
    var.mat[var.mat==0] <- rho
    # covariates
    X <- MASS::mvrnorm(n, rep(0, k), Sig = var.mat)
    X<-apply(X,2,scale)
    X[,1]<-X[,1]/(mean(X[,1]^2)^.5)
    X<-apply(X,2,FUN=function(x) x/(mean(x^2))^.5)
    
    
    ## REs ----
    ids.map<-sample(as.factor(letters[1:20]),n,T)
    ids.map <- as.factor(ceiling(sort(X[,1],ind=T)$ix/20))
    ids.map <- as.factor(ceiling(sort(rnorm(n),ind=T)$ix/20))
    
    res.map <- rnorm(length(unique(ids.map)))
    names(res.map)<-sort(unique(ids.map))
    res.true <-  (res.map[ids.map])*REs
    X.REs<-model.matrix(~ids.map-1)
    
    if(errhet)	sd.use<-((1+X[,1]+(length(unique(X[,1]))==2) )^2/2)^.5 else sd.use<-1
    errst<-rnorm(n,sd=sd.use)#rt(n,8)*sd.use#
    treat<-(X[,1])+errst+res.true
    errsy<-rnorm(n,sd=1)#rt(n,8)*sd.use#
    if(meanhet) Y<-1+treat*(X[,1]^2)+errsy else Y<-1+treat+(X[,1]^2)+errsy
    Y <- Y - res.true
    ##Shift X
    if(inter){
      m1<-exp(-5*(outer(X[,1],X[,1],"-")^2))
      diag(m1)<-0
      oY<-m1%*%X[,1]^2/rowSums(m1)
      oT<-m1%*%X[,1]^2/rowSums(m1)
      
      oY<-m1%*%(treat)/rowSums(m1)
      oT<-m1%*%(X[,1]^2)/rowSums(m1)
      
      
      Y<-Y+.5*oY
      treat<-treat+.5*oT
    }
    lm(Y~treat+X)
    
    X2.model<-X
    X2.model[,1]<-X[,1]-.5*X[,2]
    X2.model[,2]<-X[,2]-.5*X[,1]
    ## Rotates X's
    X<-X2.model
    colnames(X)<-paste("X",1:ncol(X),sep="_")
    
    output <- list(
      "Y" = Y,
      "treat" = treat,
      "X" = X,
      "X.all" = cbind(X, X.REs),
      "ids.map" = ids.map
    )
    return(output)
  }