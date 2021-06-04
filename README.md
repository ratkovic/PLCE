<!-- badges: start -->
[![R build status](https://github.com/ratkovic/PLCE/workflows/R-CMD-check/badge.svg)](https://github.com/ratkovic/PLCE/actions)
<!-- badges: end -->

# PLCE
Software for the Partially Linear Causal Effect Model, as given in Ratkovic (2021).

### Install latest version

You can install the latest version by running:
```R
devtools::install_github('ratkovic/PLCE')
```

### Troubleshooting installation

This version uses `Rcpp` extensively for speed reasons. These means you need to have the right compilers on your machine.

#### Windows
If you are on Windows, you will need to install [RTools](https://cran.r-project.org/bin/windows/Rtools/) if you haven't already. If you still are having difficulty with installing and it says that the compilation failed, try installing it without support for multiple architectures:
```R
devtools::install_github('lukesonnet/KRLS', args=c('--no-multiarch'))
```

#### Mac OSX

In order to compile the `C++` in this package, `RcppArmadillo` will require you to have compilers installed on your machine. You may already have these, but you can install them by running:

```bash
xcode-select --install
```

If you are having problems with this install on Mac OSX, specifically if you are getting errors with either `lgfortran` or `lquadmath`, then try open your Terminal and try the following:

```bash
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Also see section 2.16 [here](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf)

#### An Example


```R
library(PLCE)

# Generate data.
n <- 1000
p <- 10
X <- matrix(rnorm(n * p), n, p)

```
