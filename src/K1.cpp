// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double K1(double x, arma::vec egvalues, double q)
{
    double res = 0.0;
    const int en = egvalues.size();
    
    for(int i = 0; i < en; i++)
    {
        res = res + egvalues(i)/(1-2*egvalues(i)*x);
    }
    
    res = res - q;
    
    return res;
}

