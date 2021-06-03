#' Q Function
#'
#'

## Update beta ----
#' update beta
#'
#'
#' @noRd
#'
updatebeta <-
  function(y,
           Dtau,
           sigma.sq,
           EM,
           XprimeX,
           Xprimey,
           type,
           blockupdate,
           beta,
           trim,
           i.iter,
           burnin,
           tauthresh,
           sd.y) {
    if (length(XprimeX) == 1) {
      beta.mean <- beta <- mean(y)
      if (!EM)
        beta <- beta + sigma.sq ^ .5 * rnorm(1) / XprimeX[1] ^ .5
      return(list("beta" = beta))
    }
    
    
    
    which.update <- (abs(beta) > 1e-5 * sd.y)
    which.update[1] <- TRUE
    if (sum(which.update) == 1)
      which.update[1:2] <- TRUE
    if (!EM)
      which.update <- rep(TRUE, length(beta))
    if (trim == FALSE)
      which.update <- rep(TRUE, length(beta))
    
    
    Dtau.inv <- 1 / Dtau[which.update]
    Dtau.inv[Dtau.inv > tauthresh] <- tauthresh
    Dtau.inv[is.na(Dtau.inv)] <- 1e5
    Dtau.inv[!is.finite(Dtau.inv)] <- 1e5
    Dtau.inv <- diag(Dtau.inv)
    
    # beta[which.update] <-
    #   solve_cpp(XprimeX[which.update, which.update] + Dtau.inv, Xprimey[which.update, ])$beta
    beta[which.update] <-
      #solve(XprimeX[which.update, which.update] + Dtau.inv, Xprimey[which.update, ])
      qr.coef(qr(XprimeX[which.update, which.update] + Dtau.inv), Xprimey[which.update,])
    beta[is.na(beta)] <- 1e-6 * sd(y)
    beta.mean <- NA
    
    
    output <- list("beta" = beta, "beta.mean" = beta.mean)
  }


#' update weights via EM
#'
#' @param x a vector with m in the first element and gamma0 in the second
#'
#' @details A lookup table with pre-computed integrals, `precomp_wts`, is used internally.
#' The table is generated in the data-raw directory.
#'
#' @noRd
#'
updatewts.EM <- function(x, seq.obj0) {
  pow <- 1
  m <- x[1]
  gamma0 <- x[2]
  m <- max(m, 1e-6)
  seq.m <- seq.obj0[[1]]
  seq.gamma <- seq.obj0[[2]]
  #seq.m <- seq(-6, 10, by = 0.025)
  #seq.gamma <- seq(-6, 12, by = 0.025)
  
  logm <- log(m)
  #logm <- rcppClamp(logm,-5.99,9.95)#pmax(-5.95, pmin(logm, 9.9))
  loggamma <- log(gamma0)
  #loggamma <-rcppClamp(loggamma,-5.99,11.95)# pmax(-5.95, pmin(loggamma, 11.95))
  
  
  #ind.m <- which(abs(logm - .0125 - seq.m) == min(abs(logm - .0125 - seq.m)))
  #ind.gamma <- which(abs(loggamma - .0125 - seq.gamma) == min(abs(loggamma - .0125 - seq.gamma)))
  ind.m <-
    floor((6 + logm) * 40) + 1 #which_maxCpp(-abs(logm - .0125 - seq.m))+1
  ind.gamma <-
    floor((6 + loggamma) * 40) + 1#which_maxCpp(-abs(loggamma - .0125 - seq.gamma))+1
  
  
  m.diff <-
    c(abs(logm - seq.m[ind.m]), abs(logm - seq.m[ind.m + 1]))
  m.diff <- m.diff / sum(m.diff)
  w1 <-
    sum(c(log(precomp_wts[[ind.m]][ind.gamma]), log(precomp_wts[[ind.m + 1]][ind.gamma])) * m.diff)
  
  gamma.diff <-
    c(abs(loggamma - seq.gamma[ind.gamma]), abs(loggamma - seq.gamma[ind.gamma + 1]))
  gamma.diff <- gamma.diff / sum(gamma.diff)
  w2 <-
    sum(c(log(precomp_wts[[ind.m]][ind.gamma]), log(precomp_wts[[ind.m]][ind.gamma + 1])) * gamma.diff)
  return(exp(mean(c(w1, w2))))
}
