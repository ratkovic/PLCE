##  Function for making basis functions
makebases <-
  function(treat,
           X,
           het,
           SIS.use = NULL,
           replaceme,
           alwaysinter = NULL) {
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
      inter.schedule[inter.schedule[, 1] > inter.schedule[, 2], ]
    sample.cor <- 1:length(treat)
    sample.cor <- which(replaceme > 4)
    if (length(sample.cor) > 1000) {
      sample.cor <- sample(sample.cor, 1000, FALSE)
    }
    
    ## Rcpp implementation
    cors <-
      corbases(treat = treat[sample.cor],
               X = Xt[sample.cor, ],
               inter.schedule,
               1:length(sample.cor))$cors

    ## Sure screen ----
    n <- length(treat)
    if (length(SIS.use) == 0) {
      SIS.use <- floor(25 * (1 + n ^ (1 / 5)))
      SIS.use <- min(SIS.use, 400)
      SIS.use <- max(SIS.use, 50)
      #SIS.use<-min(SIS.use,n/4)
      #SIS.use<-floor(SIS.use/2)
    }
    
    keeps <- sort(cors, decreasing = T, ind = T)$ix[1:SIS.use]
    t1 <-
      apply(
        inter.schedule[keeps, ],
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


make.replaceme <- function(x, id) {
  n <- length(x)
  out <- NULL
  if (length(id) == 0) {
    for (i in 1:(ceiling(n / 6) + 2)) {
      out[6 * (i - 1) + 1:6] <- sample(1:6)
    }
    out <- out[1:n]
  } else{
    unique.id <- unique(id)
    for (i in unique.id) {
      out.temp <- which(id == i)
      length.temp <- length(out.temp)
      out[out.temp] <-
        sample(rep(sample(1:6), ceiling(length.temp / 6)))[1:length.temp]
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
  
  ## Try not orthogonalizing?
  
  # return(X)
  
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
    which.curr <- sort(covs.curr, decreasing = T, ind = T)$ix[1]
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
  median_cpp(c(diff1, x))
}

#' @importFrom glmnet cv.glmnet
trimbases.boot <-
  function(y,
           basesy0,
           trimboot.num = 5,
           wts.use = NULL,
           id = NULL,
           replaceme) {
    
    glmnet1 <- cv.glmnet(as.matrix(basesy0)[replaceme>4,],y[replaceme>4])
    keeps <- which(coef(glmnet1)[-1]!=0)
    if(length(keeps)==0) keeps <- 1:3
    return(basesy0[,keeps])
    
    n <- length(y)
    if (length(wts.use) == 0)
      wts.use <- rep(1, n)
    coefsy.samp <- NULL
    basesy0 <- apply(basesy0, 2, scale2)
    y <- scale2(y)
    basesy0int <- cbind(1, basesy0)
    beta.init <- NULL
    
  
    
    indsamp <- which(replaceme > 4)
    numsamp <- length(indsamp)
    numsamp <- min(numsamp, 1000)
    
    
    lastline <- 214000
    # m1 <- mget(ls())
    # save(m1, file="diagnose.Rda" ); rm(m1)
    sy0 <-
     sparsereg_GCV(y,
                    basesy0)

    # trimboot.num <- rcppClamp(max(ceiling(sy0$dof), 2) * 5, 10, 50)
    coefsy.samp <- matrix(NA, nrow=trimboot.num+1, ncol=length(sy0$coef))
    coefsy.samp[1, ] <-abs(sy0$coef)
    for (i.samp in 1:trimboot.num) {
      #samp.curr <- which(sample(1:3, n, T) == 1)
      
      samp.curr <- sample(indsamp, numsamp, TRUE)
      basesy0.boot <- basesy0[samp.curr, ]
      sds.boot <- apply(basesy0.boot, 2, sd_cpp2)
      basesy0.boot[, sds.boot < 0.01] <-
        rnorm(length(samp.curr) * sum(sds.boot < 0.01))
      
      lastline <- 233000
      # m1 <- mget(ls())
      # save(m1, file="diagnose.Rda" ); rm(m1)
      sy0 <-
        sparsereg_GCV(y[samp.curr],
                      basesy0.boot)

        coefsy.samp[i.samp+1, ] <- abs(sy0$coef)
      }
    
    cutoff.sis <-
      mean(sapply(
        1:500,
        FUN = function(x)
          max(abs(rnorm(trimboot.num+1))) / n ^ .5
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
  col.sched <- col.sched[col.sched[, 1] != col.sched[, 2], ]
  X.svd <- NULL
  for (i in 1:nrow(col.sched)) {
    inter.temp <- as.vector(apply(X0[, col.sched], 1, prod))
    X.svd <- cbind(X.svd, svd(X0[, col.sched])$u[, 1:2])
  }
  #svd.X <- svd(X)
  #spec <- cumsum(svd(X)$d) / sum((svd(X)$d))
  #keeps <- rep(T, ncol(X))
  #keeps[keeps > 0.9] <- F
  #X.svd <- cbind(X, apply(svd.X$u, 2, scale)[, keeps])
  X.svd
}

## Construct bases ----
generate.bases <-
  function(y2.b,
           y,
           basesy0,
           X,
           id = NULL,
           replaceme,
           alwaysinter = NULL) {
    # tic("starting generate.bases")
    n <- length(y2.b)
    X <- basesy0
    if (length(id) > 0) {
      X.dum <- cbind(sample(y), sample(y))
      colnames(X.dum) <- c("X1", "X2")
      X.dum <- apply(
        X.dum,
        2,
        FUN = function(z)
          fastLm(z ~ cbind(1, y))$res
      )
      y <- y2.b <- y2.b - sparsereg(y2.b, X.dum, id = NULL)$fit
    }
    # toc()
    #makebases <- function(treat, X, het, SIS.use = NULL, replaceme)
    ## Stage 1 SIS
    # tic("make bases stage 1")
    basesy0 <-
      suppressMessages(makebases(
        y2.b,
        basesy0,
        het = F,
        SIS.use = NULL,
        replaceme,
        alwaysinter
      )$bases)
    keeps.try <- check.cor(basesy0, thresh = 0.0000001)$k
    keeps.try[1:ncol(X)] <- T
    #basesy0<-basesy0[,keeps.try]
    
    basesy <- cbind(X, basesy0[, keeps.try])
    keeps.try <- check.cor(basesy, thresh = 0.000001)$k
    keeps.try[1:ncol(X)] <- T
    basesy <- as.matrix(basesy[, keeps.try])
    basesy <- orthog.me(y, basesy)
    basesy <- apply(basesy, 2, scale2)
    # toc()
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
    
    # tic("trimbases.boot")
    basesy <- trimbases.boot(y, basesy, replaceme = replaceme)
    # toc()
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
  if (n > 200) {
    sample.use <- sample(1:n, 200, FALSE)
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
                     sample.use.0 = sample.use) {
        z1 <- scale2(X.interf.mat.0[sample.use.0, x[1]])
        z2 <- scale2(X.interf.outer.0[sample.use.0, x[2]])
        create.inter.basis(theta = rot.est, resy.0[sample.use.0], z1, z2, replaceme.0[sample.use.0])$logcor
      }
    )
  ind2 <- sort(cors.screen, decreasing = T, index = T)$ix[1:10]
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
                 0)
  cor.out <- log(1+abs(cor1))
  # cor.out <- log(max(1e-5, cor1))
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
    
    replaceme.temp <- replaceme
    if (sum(replaceme > 4) > 1000) {
      which.gt4 <- which(replaceme.temp > 4)
      replaceme.temp[replaceme.temp > 4] <- 0
      replaceme.temp[sample(which.gt4, 1000, FALSE)] <- 5
    }
    
    if (length(replaceme) > 3000) {
      which.zero <- sample(1:length(replaceme), 3000, FALSE)
      replaceme.temp[which.zero] <- 0
    }
    
    
    rot.est <-
      log(1.06 * (sum(replaceme.temp > 4)) ^ -.2 * (1 / 3) ^ -.2)
    
    for (i.interf in 1:ncol(X.interfy$Xvars)) {
      opt.y <-
        optimize(
          cor1,
          interval = c(-2, 2) + rot.est,
          resy = resy[replaceme.temp != 0],
          res1.2 = X.interfy$Xvars[replaceme.temp != 0, i.interf],
          x = X.interfy$Xmatvars[replaceme.temp != 0 , i.interf],
          replaceme.temp[replaceme.temp != 0],
          maximum = T,
          tol = 0.001
        )
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
    out <- (z - mean(z)) / sdz
    #out <- as.vector(scale(z))
  }
  out
}

bs.me_mdei <- function(x, varname = "X") {
  x <- x - mean(x)
  if (length(unique(x)) <= 2){
    m1 <- as.matrix(x)
    m1 <- apply(m1, 2, scale)
    colnames(m1)<-paste(varname,"linear",sep="_")
    return(m1)
  }
  if (length(unique(x)) <= 3){
    m1<-as.matrix(cbind(x, x ^ 2))
    m1 <- apply(m1, 2, scale)
    colnames(m1)<-paste(varname,c("linear","quadratic"),sep="_")
    return(m1)
  }
  if (length(unique(x)) <= 4){
    m1<-cbind(x, x ^ 2, x ^ 3-3*x^2)
    m1 <- apply(m1, 2, scale)
    colnames(m1)<-paste(varname,c("linear","quadratic","cubic"),sep="_")
    return(m1)
  }
  
  b1 <- bSpline2(x, df = 3)
  colnames(b1) <- paste(varname,"spline",1:ncol(b1),sep="_")
  x2 <- scale(x)
  #m1 <-
  #  cbind(x, b1, ibs(x, df=3)[,-3], ibs(x, df=5)[,-5])
  
  m1 <-
    cbind(x,  bSpline2(x, df = 3),  ibs(x, df=3)[,-3])#, ibs(x, df=5)[,-5])    
  colnames(m1)[1]<-paste(varname,"linear",sep="_")
  x2 <- scale(x)
  sd.x <- sd(x)
  x3 <- x2-1
  x4 <- x2+1
  m2 <- cbind(#x,b1,
    #m1,
    x,
    x2^2-1, x2^3-3*x2, x2^4-6*x2^2+3,
    x3^2-1, x3^3-3*x3, x3^4-6*x3^2+3,
    x4^2-1, x4^3-3*x4, x4^4-6*x4^2+3
  )
  m1 <- cbind(m2)
  m1 <- apply(m1, 2, scale)
  #m1 <- svd(m1)$u
  #m1 <- apply(m1, 2, scale)
  colnames(m1) <- paste(varname,"spline",1:ncol(m1),sep="_")
  return(m1)
}


bs.me <- function(x, degree = 5) {
  x<-scale(x)
  basis.out <-
    cbind(x, bSpline2(x, degree=3, knots=0))#, bSpline2(x, degree=5, knots=0))
  basis.out <- basis.out[, check.cor(basis.out, 0.0001)$k]
  apply(as.matrix(basis.out),2,FUN=function(x) x/sd(x))
  
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




cleanNAs <- function(X) {
  return(X)
  if (length(X) == 0)
    return(X)
  X <- as.matrix(X)
  if (ncol(X) < 2)
    return(as.matrix(X))
  return(as.matrix(X[, colSums(is.na(X)) == 0]))
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
  obj <- obj / sd(as.vector(obj), na.rm=T)
  
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
  drops.vec <- is.na(vec1)
  vec1 <- vec1[!drops.vec]
  expected <- expected[!drops.vec]#^(bc2$lambda)
  ind <- sort(vec1, decreasing = F, ind = T)
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
  X <- cleanNAs(X)
  keeps.out <- rep(FALSE, ncol(X))
  dropints <- apply(X, 2, sd_cpp2) == 0
  X <- X[,!dropints]
  
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

sparsereg_GCV <- function(y0, X0, id0 = NULL, usecpp = TRUE) {
  norms <- which(apply(X0, 2, sd)==0)
  for(i.fill in norms){
    x.temp <- rnorm(length(y0))
    x.temp <- lm(x.temp~y0)$res
    X0[,i.fill] <- x.temp
  }
  model.out <- sparsereg(y = y0, X = X0, id = id0)
  return(model.out)
  
  n <- length(y0)
  p <- ncol(X0)
  if (ncol(as.matrix(X0)) < 3) {
    model.out <- sparsereg(y = y0, X = X0, id = id0)
    return(model.out)
  }
  if (usecpp == FALSE) {
    alpha.max <- max(n * log(ncol(X0)), ncol(X0) * 1.25)
    alpha.min <- min(ncol(X0) * 1.25, n * log(ncol(X0)) / 2)
    alpha.schedule <- seq((alpha.max), alpha.min, length = 10)
    GCV.out <- sapply(
      alpha.schedule,
      FUN = function(z)
        tryCatch(
          sparsereg(
            y = y0,
            X = X0,
            id = id0,
            edf = T,
            alpha.prior = "custom",
            alpha.use = z
          )$GCV,
          error = function(h)
            1e7,
          warning = function(h)
            1e7
        )
    )# Close out sapply
    # a.min <-alpha.schedule[GCV.out==min(GCV.out)]
    lm1 <- lm(GCV.out ~ alpha.schedule + I(alpha.schedule ^ 2))
    alpha.lm <- -(lm1$coef[2]) / (2 * lm1$coef[3])
    alpha.min <- ifelse(lm1$coef[3] > 0 ,-(lm1$coef[2]) / (2 * lm1$coef[3]),
                        alpha.schedule[which.min(GCV.out)])
    alpha.min <-
      rcppClamp(alpha.min, alpha.schedule[1], rev(alpha.schedule)[1])
    model.out <-
      sparsereg(
        y = y0,
        X = X0,
        id = id0,
        edf = T,
        alpha.prior = "custom",
        alpha.use = alpha.min
      )
    
  } else{
    ## New using Cpp bayeslasso to find alpha!
    alpha.schedule <-
      seq(log(8 * n * log(p)), log(p), length = 8)
    alpha.schedule <-
      seq(log(8 * n * log(p)), log(p), length = 5)
    #alpha.max <- max(n*log(ncol(X0)), ncol(X0)*1.25)
    #alpha.min <- min(ncol(X0)*1.25, n*log(ncol(X0))/2)
    #alpha.schedule <- log(seq((alpha.max), alpha.min,length=10))
    X0.temp <- apply(X0, 2, scale2)
    gcv.out <-
      log(sapply(
        alpha.schedule,
        FUN = function(z)
          bayesLasso(y0, cbind(1, X0.temp), exp(z))$GCV
      ))
    lm1 <- lm(gcv.out ~ alpha.schedule + I(alpha.schedule ^ 2))
    alpha.lm <- -(lm1$coef[2]) / (2 * lm1$coef[3])
    alpha.min <- ifelse(lm1$coef[3] > 0,-(lm1$coef[2]) / (2 * lm1$coef[3]),
                        alpha.schedule[which.min(gcv.out)])
    alpha.min <-
      rcppClamp(alpha.min, alpha.schedule[1], rev(alpha.schedule)[1])
    #print(alpha.min)
    #model.out <- sparsereg(y=y0,X=X0, id=id0, edf=T, alpha.prior = "custom", alpha.use =exp( alpha.min ))
    model.out <-
      sparsereg(
        y = y0,
        X = X0,
        id = id0,
        edf = T,
        alpha.prior = "custom",
        alpha.use = exp(alpha.min)
      )
    
  } # Closes out if
  model.out
}




make.soe <-
  function(basesy.all.inner,
           id,
           sy,
           y,
           fits.y,
           fits.REs.y,
           replaceme) {
    if (length(sy$chol.inv) <= 1)
      return(NULL)
    basesy.all.temp <- cbind(1, basesy.all.inner)
    basesy.all.temp <- cbind(basesy.all.temp[, sy$beta != 0])
    
    mat.all <- basesy.all.temp %*% (sy$chol.inv)
    
    mat.control <- cbind(1, fits.REs.y, fits.y)
    weights.est <- 1 * (replaceme < 3)
    y2 <- fastLm(y ~ mat.control, weights = weights.est)$res
    mat.all2 <-
      apply(
        mat.all,
        2,
        FUN = function(z)
          scale2(fastLm(z ~ mat.control, weights = weights.est)$res)
      )
    
    mat.orthog <-
      as.matrix(orthog.me(y2, mat.all2, weights.lm = weights.est))
    if (ncol(mat.orthog) == 1)
      return(mat.orthog[weights.est == 0, ])
    colnames(mat.orthog) <- paste("X", 1:ncol(mat.orthog), sep = "_")
    
    svd.temp <- svd(mat.orthog[weights.est == 1, ])
    mat.orthog2 <- mat.orthog %*% t(svd.temp$v)
    #print(svd.temp$d/svd.temp$d[1])
    keeps <- (cumsum(svd.temp$d ^ 2) / sum(svd.temp$d ^ 2)) < 1
    #print(c(ncol(basesy.all.inner),sum(keeps)))
    if (sum(keeps) > .75 * sum(replaceme < 3)) {
      maxcol <- ceiling(.75 * sum(replaceme < 3))
      keeps[-c(1:maxcol)] <- FALSE
    }
    mat.orthog2[weights.est == 0, keeps]
  }


allbases <- function(y,
                     y2.b,
                     treat,
                     treat.b,
                     treat.y,
                     X,
                     id,
                     replaceme,
                     fit.interference,
                     fit.treatment.heteroskedasticity,
                     inter.schedule.obj = NULL) {
  n <- length(y)
  replaceme <- make.replaceme(y, id)
  
  #generate.bases <- function(y2.b, y, basesy0, X, id=NULL, replaceme) {
  ## Make bases for outcome, treatment----
  #cat("#####  Step 1: Constructing conditional mean bases\n")
  # tic("generate.bases")
  basest0 <- generate.bases(treat.b, treat, X, X, id = NULL, replaceme)
  basest0 <- apply(basest0, 2, scale2)
  sds <- apply(basest0, 2, FUN=function(z) min(tapply(z,replaceme, sd)))
  basest0 <- basest0[,sds>0]
  lastline <- 949000
  # m1 <- mget(ls())
  # save(m1, file="diagnose.Rda" ); rm(m1)
  treat.b.2 <-
    treat.b - basest0 %*% sparsereg_GCV(treat.b[replaceme > 4], basest0[replaceme >
                                                                          4, ])$coef#,   id=NULL)$fitted
  #basest0 <- cbind(basest0, generate.bases(treat.b.2, treat, X, X, id=NULL,replaceme))
  basest0 <- cleanNAs(basest0)
  
  basesy0 <-
    generate.bases(y2.b, y2.b, cbind(X), cbind(X), id = NULL, replaceme)
  sds <- apply(basesy0, 2, FUN=function(z) min(tapply(z,replaceme, sd)))
  basesy0 <- basesy0[,sds>0]
  y2.b.2 <- y2.b - sparsereg(y2.b, basesy0,   id = NULL)$fitted
  #basesy0.2 <- generate.bases(y2.b.2, y2.b.2, cbind(X), cbind(X),id=NULL, replaceme) ## Don't need REs--already partialed out!
  #basesy0<-cbind(basesy0,basesy0.2)
  basesy0 <- cleanNAs(basesy0)
  
  ## Strip out linear dependencies--it just helps
  keeps.t <- which(!is.na(lm(treat ~ basest0)$coef[-1]))
  keeps.y <- which(!is.na(lm(y ~ treat + basesy0)$coef[-c(1:2)]))
  
  basest0 <- as.matrix(basest0[, keeps.t])
  basesy0 <- as.matrix(basesy0[, keeps.y])
  #
  
  baset0 <- basest0[, check.cor(basest0)$k]
  basey0 <- basesy0[, check.cor(basesy0)$k]
  # toc()
  
  # tic("orthog.me generate bases")
  basest0 <- orthog.me(treat, basest0, weights.lm = 1 * (replaceme > 4))
  basesy0 <- orthog.me(y, basesy0, weights.lm = 1 * (replaceme > 4))
  # toc()
  
  ## Initial outcome and propensity model ----
  sy <- sparsereg(y,
                  cbind(basesy0, basest0),
                  
                  
                  id = NULL)
  res1.y <- as.vector(scale2(lm(y ~ sy$fit)$res))
  
  st <- sparsereg(treat, basest0,   id = NULL)
  res1 <- as.vector(scale2(treat - st$fit))
  
  ## Make residual nonparametric terms ----
  res2.b <- (res1)
  # generate.bases <- function(y2.b, y, basesy0, X, id=NULL, replaceme) {
  basest2.0 <- NULL
  res1.2 <- scale2(res1)
  
  if (fit.treatment.heteroskedasticity) {
    # tic("treatment hetero")
    basest2.0 <-
      generate.bases(res2.b,
                     res2.b,
                     X,
                     X,
                     id = NULL,
                     replaceme,
                     alwaysinter = res2.b)
    basest2.0 <- cleanNAs(basest2.0)
    
    st.2 <-
      sparsereg(res1,
                apply(
                  basest2.0,
                  2,
                  FUN = function(z)
                    scale2(z) * scale2(res1)
                ),
                id = NULL)
    res1.2 <- scale2(res1 - st.2$fit)
    
    res2.b.2 <- (res1.2)
    
    basest2.0 <-
      cbind(
        basest2.0,
        generate.bases(
          res2.b.2,
          res2.b.2,
          X,
          X,
          id = NULL,
          replaceme,
          alwaysinter = res2.b.2
        )
      )
    basest2.0 <- cleanNAs(basest2.0)
    basest2.0 <-
      orthog.me(res1.2, basest2.0, weights.lm = 1 * (replaceme > 4))
    # toc()
  }
  
  ## Split sample ----
  replaceme <-
    replaceme0 <- rep(c(1:6), times = floor(n / 2))[1:(n + 4)]
  replaceme <- make.replaceme(y, id)#sample(replaceme[1:n])
  
  ## Make interference bases ----
  X.interfy <- X.interft <- NULL
  if (fit.interference) {
    # tic("interference bases")
    X.interfy <- generate.Xinterf(res1.y, X, treat, replaceme)
    X.interft <- generate.Xinterf(res1.2, X, NULL, replaceme)
    # toc()
  }
  
  ## Generate splits ----
  # fits.replaceme <- lm(lm(res1~st.2$fit)$res~X.interft$Xbasis)$fit#fastLm(treat ~ st$fit)$res
  # fits.replaceme.y <- fastLm(y ~ sy$fit)$res
  
  if (length(basesy0) == 0)
    basesy0 <- cbind(rnorm(n), rnorm(n))
  if (length(basest0) == 0)
    basest0 <- cbind(rnorm(n), rnorm(n))
  if (length(basest2.0) == 0)
    basest2.0 <- cbind(rnorm(n), rnorm(n))
  
  replaceme3<-ceiling(replaceme/2)
  for(i.r in 1:3){
    basesy0[replaceme3==i.r,]<-apply(basesy0[replaceme3==i.r,],2,scale2)
    basest0[replaceme3==i.r,]<-apply(basest0[replaceme3==i.r,],2,scale2)
    basest2.0[replaceme3==i.r,]<-apply(basest2.0[replaceme3==i.r,],2,scale2)
    if(fit.interference){
      if(is.matrix(X.interfy)){
        print(dim(X.interfy))
       X.interfy = as.matrix(X.interfy)
       X.interfy[replaceme3==i.r,]<-apply(X.interfy[replaceme3==i.r,],2,scale2)
      }
      if(is.matrix(X.interft)){
        X.interft = as.matrix(X.interft)
        X.interft[replaceme3==i.r,]<-apply(X.interft[replaceme3==i.r,],2,scale2)
      }
     }
    
  }
  #if(abs(cor(fits.replaceme.y,y))>abs(cor(fits.replaceme,treat))) fits.replaceme<-fits.replaceme.y
  output <- list(
    "basesy0" = basesy0,
    "basest0" = basest0,
    "basest2.0" = basest2.0,
    "X.interfy" = X.interfy,
    "X.interft" = X.interft,
    "replaceme" = replaceme
  )
}

sd_cpp2 <- function(z)
  sd(z)
