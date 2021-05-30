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


//' Check Spearman correlations between interactions in X and treatment
//' 
//' @param y A vector of outcomes.
//' @param X A matrix of spline bases.
//' @param alpha.schedule The prior on lambda
//' @export

// [[Rcpp::export]]
Rcpp::List bayesLasso(arma::vec y, 
                      arma::mat X,
                      arma::vec alpha
) {
    
    // Declare inputs 
    int n = X.n_rows;
    int p = X.n_cols;
    arma::vec Etausqinv = arma::ones(p);
    arma::vec Ewtsqtausq = arma::ones(p);
    arma::mat XpX = X.t()*X;
    arma::mat Xpy = X.t()*y;
    arma::mat XpXsolve = XpX;
    arma::vec fits = arma::zeros(n);
    arma::vec lambda = arma::ones(1);
    lambda(0) = 1;
    
    double sdy = arma::stddev(y);
    double prec = 1/pow(sdy,2);
    double sigma_sq = 1/prec;
    double conv = 1;
    double lambda_temp = 1;
    double edf = 0;
    double GCV = 0;
    
    arma::vec beta = arma::zeros(p)+sdy*10;
    arma::vec beta_last = beta;
    for(int i_outer = 0; i_outer < 5; i_outer++){
        if(conv > 0.0000001*sdy){
            //*  Update XpX
            for(int i_p = 1; i_p<p; i_p++){
                XpXsolve(i_p,i_p) = XpX(i_p,i_p)+Etausqinv(i_p);
            }
            
            //* Update beta along nonzero indices
            beta_last = beta;
            arma::uvec update_ind = arma::find(arma::abs(beta) > 0.0001*sdy );
            beta(update_ind) = arma::solve(XpXsolve.submat(update_ind,update_ind),Xpy.rows(update_ind));
            fits = X*beta;
            
            //* Update Etausq and Etausqinv
            for(int i_p = 1; i_p<p; i_p++){
                Ewtsqtausq(i_p)  = abs(beta(i_p))/lambda(0)*sqrt(prec)+pow(lambda(0),-2);
                Etausqinv(i_p)  = lambda(0)/abs(beta(i_p))*sqrt(sigma_sq);  
            }
            // Ewtsqtausq  = arma::abs(beta)/lambda(0)*sqrt(prec)+pow(lambda(0),-2);
            // Etausqinv  = lambda(0)/arma::abs(beta)*sqrt(sigma_sq);
            Etausqinv(0) = 0;
            Ewtsqtausq(0) = 0;
            // Ewtsqtausq <- pmax(abs(beta), 1e-5) / lambda * sqrt(prec) * wts + 1 / lambda.sq
            //  Etausqinv <- lambda / pmax(abs(beta), 1e-5) * sqrt(sigma.sq) * wts
            
            
            //* Update lambda
            lambda_temp = sqrt(
                (alpha(0) - 1) / (sum(Ewtsqtausq) / 2 + 1)
            );
            lambda(0) = lambda_temp;
            
            //* Update sigma_sq
            prec = (n/2+(p-1)/2+1)/(
                sum((y-fits)%(y-fits))+sum(beta%beta% Etausqinv/2)
            );
            sigma_sq = 1/prec;
            // prec <- (
            //    ( n/2+ p/2 +1)/
            //       ( sum((y-fits.curr-fits.res)^2) +sum(beta[2:(p+1)]^2 / (Dtau[2:(p+1)])) / 2)
            // )
            
            
            conv = max(arma::abs(beta-beta_last));
        }//close out while loop
        arma::uvec update_ind = arma::find(arma::abs(beta) > 0.0001*sdy );
        edf = arma::trace(XpX.submat(update_ind,update_ind)*arma::pinv(XpXsolve.submat(update_ind,update_ind)));
        GCV = sum((y-fits)%(y-fits))/((n-log(n)/2*edf)*(n-log(n)/2*edf));
    }// Close out i_outer
    
    return Rcpp::List::create(
        Rcpp::Named("intercept") = beta(0),
        Rcpp::Named("coefficients") = beta.rows(1,p-1),
        Rcpp::Named("fitted.values") = fits,
        Rcpp::Named("GCV") = GCV,
        Rcpp::Named("Etausqinv") = Etausqinv
        
        
    );
}
