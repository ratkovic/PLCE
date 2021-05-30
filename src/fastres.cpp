// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// fastLm.cpp: Rcpp/Armadillo glue example of a simple lm() alternative
//
// Copyright (C)  2010 - 2017  Dirk Eddelbuettel, Romain Francois and Douglas Bates
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;


//' Fast partialing out using weighted least squares
//' 
//' @param X A single integer.
//' @param y A single integer.
//' @param w A single integer.
//' @export
// [[Rcpp::export]]

List fastres_cpp(const arma::mat& X, const arma::mat& y, const arma::colvec& w) {

    int nc = X.n_cols;
    int nr = X.n_rows;
    int ncy = y.n_cols;
    
    arma::colvec sqw=sqrt(w);
    arma::mat Xw= X;//%sqrt(w);
    arma::mat yw= y;
    
    
    for(int i_nr = 0; i_nr<nr; i_nr++){
         for(int i_nc = 0; i_nc<nc; i_nc++){
            Xw(i_nr,i_nc) = X(i_nr,i_nc)*sqw(i_nr,0);
            }
         for(int i_ncy = 0; i_ncy<ncy; i_ncy++){
             yw(i_nr,i_ncy)=y(i_nr,i_ncy)*sqw(i_nr,0);
         }
    }
    arma::mat coef = arma::pinv(Xw.t()*Xw)*(Xw.t()*yw);    // fit model y ~ X
    arma::mat res  = y - X*coef;           // residuals
    
    
    return List::create(Named("res") = res
                            );
}


// [[Rcpp::export]]
NumericVector rcppClamp(NumericVector x, double mi, double ma) {
    return clamp(mi, x, ma);
}



// [[Rcpp::export]]
int which_maxCpp(NumericVector v) {
    int z = which_max(v);
    return z;
}


// [[Rcpp::export]]
arma::vec arma_var(const arma::vec& v1) {
    return arma::cov(v1);
};

// [[Rcpp::export]]
arma::mat arma_matvar(const arma::mat& v1) {
    return arma::cov(v1);
};

// [[Rcpp::export]]
arma::mat arma_ginv(const arma::mat& v1) {
    return arma::pinv(v1);
};

// [[Rcpp::export]]
arma::mat arma_symm(const arma::mat& v1) {
    return .5*(v1.t()+v1);
};
