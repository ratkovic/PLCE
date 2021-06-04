library(PLCE)

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
