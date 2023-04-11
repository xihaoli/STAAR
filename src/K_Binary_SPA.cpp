// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double K_Binary_SPA(double x, arma::vec muhat, arma::vec G)
{
    double res = 0.0;
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
		res = res - x * muhat(i) * G(i);
        res = res + log(1 - muhat(i) + muhat(i) * exp(x * G(i)));
    }

    return res;
}

