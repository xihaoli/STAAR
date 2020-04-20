// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double K2(double x, arma::vec egvalues)
{
    double res = 0.0;
    const int en = egvalues.size();

    for(int i = 0; i < en; i++)
    {
        res = res + pow(egvalues(i),2)/pow(1-2*egvalues(i)*x,2.0);
    }

    res = res*2;

    return res;
}

