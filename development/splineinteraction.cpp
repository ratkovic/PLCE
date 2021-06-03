// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace Rcpp;
using namespace RcppArmadillo;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::List corbases(arma::vec treat, 
                       arma::mat X,
                       arma::mat inter_schedule,
                       arma::vec one_to_n
                       ) {
  
  // Declare inputs 
  int nr = inter_schedule.n_rows;
  arma::vec output = arma::zeros(nr);
  arma::vec tempbasis = treat;
  arma::vec rankvec = treat;
  
  arma::vec treat_ranked =  one_to_n(sort_index(sort_index(treat)));
  int col1 = 0;
  int col2 = 0;
  
  //Start Loop loop
  for(int i_nr = 0; i_nr<nr; i_nr++){
  col1 = inter_schedule(i_nr,0);
  col2 = inter_schedule(i_nr,1);
  tempbasis = X.col(col1-1)%X.col(col2-1);
  rankvec = one_to_n(sort_index(sort_index(tempbasis)));
    output.row(i_nr) = arma::abs(arma::cor(treat_ranked,rankvec ));
  }
  return Rcpp::List::create(
    Rcpp::Named("cors") = output
  );
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
 corbases((treat),X,inter.schedule,1:n)$cors[1:10]
  cors(inter.schedule)[1:10]
  z1<-rnorm(20);z2<-rnorm(20)
  corbases((z1),cbind(1,z2),matrix(c(1,2),nr=1),1:20)
  corSpearman(z1,z2)
  microbenchmark("cpp"=corbases((treat),X,inter.schedule,1:n),"R"=cors(inter.schedule),times=10 )
corSpearman(X[,1],X[,2])
  cor(rank(X[,1]),rank(X[,2]))
  
  microbenchmark("cpp"= corbases((z1),cbind(1,z2),matrix(c(1,2),nr=1),1:20),"cpp_notme"= corSpearman(z1,z2),times=1000 )
  
  */










