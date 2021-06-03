# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace RcppArmadillo ;


//' Update beta
//' 
//' @noRd
// [[Rcpp::export()]]
Rcpp::List solve_cpp(arma::mat XpX,
                          arma::colvec Xpy
) {
  
 arma::mat output = solve(XpX,Xpy);
 
  return Rcpp::List::create(Rcpp::Named("beta") = output);
}


//' Armadillo Standard Deviation
//' 
//' @param v1 vector to take sd
//' @noRd
// [[Rcpp::export()]]
double sd_cpp(arma::colvec v1
) {
  
  arma::mat output = arma::pow(cov(v1),.5);
  return output[0];
}



//' Armadillo Median
//' 
//' @noRd
// [[Rcpp::export()]]
double median_cpp(arma::colvec v1
) {
  
  return(arma::median(v1));
}