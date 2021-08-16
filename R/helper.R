if(F){
Rcpp::compileAttributes('~/Dropbox/Github/PLCE')           # this updates the Rcpp layer from C++ to R
roxygen2::roxygenize('~/Dropbox/Github/PLCE', roclets="namespace", clean=TRUE)  # this updates the documentation based on roxygen comments
roxygen2::roxygenize('~/Dropbox/Github/PLCE')  # this updates the documentation based on roxygen comments
}