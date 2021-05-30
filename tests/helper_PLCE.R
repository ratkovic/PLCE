Rcpp::compileAttributes('~/Dropbox/InfluenceFunctions/APSRsubmission/02_Resubmission/04a_Replication/Code/PLCE')           # this updates the Rcpp layer from C++ to R
roxygen2::roxygenize('~/Dropbox/InfluenceFunctions/APSRsubmission/02_Resubmission/04a_Replication/Code/PLCE', roclets="namespace", clean=TRUE)  # this updates the documentation based on roxygen comments
