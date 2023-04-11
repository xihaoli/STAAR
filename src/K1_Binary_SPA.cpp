// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double K1_Binary_SPA(double x, arma::vec muhat, arma::vec G, double q)
{
    double res = 0.0;
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
        res = res - muhat(i) * G(i);
		res = res + muhat(i) * G(i)/(muhat(i) + (1 - muhat(i))*exp(-x * G(i)));
    }
	
    res = res - q;

    return res;
}

