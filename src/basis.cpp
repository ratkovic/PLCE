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

//' Create a basis vector
//' 
//' @param matvec A single integer.
//' @param resvec A single integer.
//' @param onesvec A single integer.
//' @param theta A single integer.
//' @export
// [[Rcpp::export]]

Rcpp::List basis_rcpp(const arma::colvec & matvec,
                     const arma::colvec & resvec,
                     const arma::colvec & onesvec,
                     const arma::colvec & theta
                ) {
    
    
    //Gather lengths
    int nr = matvec.n_elem;
    //int nc = vec2.n_elem;
    double exptheta = exp(theta[0]);
    
    arma::mat kern = arma::zeros(nr,nr);
    
    for(int i_nr = 0; i_nr<nr; i_nr++){
        for(int i_nc = 0; i_nc<nr; i_nc++){
            if(i_nr != i_nc) kern(i_nr,i_nc) = std::exp(-exptheta*std::pow(matvec[i_nr]-matvec[i_nc],2));
        }
    }
    

    arma::vec numer = kern*resvec;
    arma::vec denom = kern*onesvec;
    
    arma::vec basis = numer/denom;

    return Rcpp::List::create(
        Rcpp::Named("basis") = basis
    );
}


//' Create a basis vector but without making a kernel matrix
//' 
//' @param matvec A single integer.
//' @param resvec A single integer.
//' @param onesvec A single integer.
//' @param theta A single integer.
//' @export
// [[Rcpp::export]]


Rcpp::List basisvec_rcpp(const arma::colvec & matvec,
                       const arma::colvec & resvec,
                       const arma::colvec & onesvec,
                       const arma::colvec & theta
 ) {


     //Gather lengths
     int nr = matvec.n_elem;
     //int nc = vec2.n_elem;
     double exptheta = exp(theta[0]);
     arma::vec numer = arma::zeros<arma::vec>(nr);
     arma::vec denom = arma::zeros<arma::vec>(nr);
     arma::vec basis = arma::zeros<arma::vec>(nr);

     double powcurr = 0;

     for(int i_nr = 0; i_nr<nr; i_nr++){
         for(int i_nc = 0; i_nc<nr; i_nc++){
             powcurr = -exptheta*pow(matvec[i_nr]-matvec[i_nc],2);
             if( (i_nr!=i_nc)  & (onesvec(i_nc,0)>0) ){
                 numer(i_nr,0) += exp(powcurr)*resvec(i_nc,0)*onesvec(i_nc,0);
                 denom(i_nr,0) += exp(powcurr)*onesvec(i_nc,0);
             }
         }
         basis(i_nr,0) = numer(i_nr,0)/denom(i_nr,0);
     }


     return Rcpp::List::create(
         Rcpp::Named("basis") = basis
     );
 }

//' Check Spearman correlations between interactions in X and treatment
//' 
//' @param treat A vector of outcomes.
//' @param X A matrix of spline bases.
//' @param inter.schedule An interaction schedule.
//' @param one_to_n A vector of length 1:length(treat).
//' @export

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