// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
double fr_R_cpp( arma::vec x,  arma::mat y, double R){
  
  int y1 = y.n_rows;
  int y2 = y.n_cols;
  
  arma::vec val(y2);
  double tot=0;
  
  for (int i = 0; i < y1; i++) {
    for (int j = 0; j < y2; j++){
      val(j)=sin(R*(x(j)-y(i,j)))/(arma::datum::pi*(x(j)-y(i,j)));
    }
    tot+=prod(val);
  }
  
  double out= tot/y1;
  return out;
  }

// [[Rcpp::export]]
List fr_R_cpp_vec(arma::vec x, arma::mat y, double R) {

    int y1 = y.n_rows;
    int y2 = y.n_cols;

    arma::vec val(y2);
    NumericVector out2;

    double tot = 0;

    for (int i = 0; i < y1; i++) {
        for (int j = 0; j < y2; j++) {
            val(j) = sin(R * (x(j) - y(i, j))) / (arma::datum::pi * (x(j) - y(i, j)));
        }
        double prod_val = prod(val);
        tot += prod_val;
        out2.push_back(prod_val);
    }

    double out = tot / y1;
    return(List::create(
        _["val"] = out,
        _["vec"] = out2
        ));
}

// [[Rcpp::export]]
arma::vec fr_Rm_cpp( arma::vec x,  arma::mat y, arma::vec R){
  
  int y1 = y.n_rows;
  int y2 = y.n_cols;
  int m= R.size();
  
  arma::vec val(y2);
  arma::vec  tot(m,  arma::fill::zeros);
  
  for (int i = 0; i < y1; i++) {
    for(int k=0;k<m;k++){
    for (int j = 0; j < y2; j++){
      val(j)=sin(R(k)*(x(j)-y(i,j)))/(arma::datum::pi*(x(j)-y(i,j)));
      }
    tot(k)+=prod(val);
    }

  }
  
  arma::vec out= tot/y1;
  return out;
}