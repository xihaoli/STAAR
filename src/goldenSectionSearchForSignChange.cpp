// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// declare K1_Binary_SPA (first derivative)
double K1_Binary_SPA(double x, arma::vec muhat, arma::vec G, double q);
// declare K1_Binary_SPA_alt (first derivative, alternative way if not converge)
double K1_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G, double q);
// check na
bool check_is_na(double x); 
// check same sign
bool haveSameSign(double a, double b);

// [[Rcpp::export]]
arma::vec goldenSectionSearchForSignChange(double a, double b, arma::vec muhat, arma::vec G, double q, double tol, int max_iter)
{
	int iter = 0;

    // results init
	arma::vec x;
	x.ones(2);
	x(0) = a;
	x(1) = b;

	// K1 init
	double K1x1 = 0.0;
	double K1x0 = 0.0;
	
	double phi = (1 + sqrt(5)) / 2;
	
	K1x0 = K1_Binary_SPA(x(0), muhat, G, q);
	if((R_finite(K1x0)==0)||(check_is_na(K1x0)))
	{
		K1x0 = K1_Binary_SPA_alt(x(0), muhat, G, q);
	}
	
	K1x1 = K1_Binary_SPA(x(1), muhat, G, q);
	if((R_finite(K1x1)==0)||(check_is_na(K1x1)))
	{
		K1x1 = K1_Binary_SPA_alt(x(1), muhat, G, q);
	}
	
	while ((iter < max_iter) && (fabs(x(1) - x(0)) > tol) && haveSameSign(K1x0,K1x1)) 
	{
		x(1) = b - (b - a) / phi;
        x(0) = a + (b - a) / phi;
		
		K1x0 = K1_Binary_SPA(x(0), muhat, G, q);
		if((R_finite(K1x0)==0)||(check_is_na(K1x0)))
		{
			K1x0 = K1_Binary_SPA_alt(x(0), muhat, G, q);
		}
	
		K1x1 = K1_Binary_SPA(x(1), muhat, G, q);
		if((R_finite(K1x1)==0)||(check_is_na(K1x1)))
		{
			K1x1 = K1_Binary_SPA_alt(x(1), muhat, G, q);
		}
	
        if (K1x1 < K1x0) 
		{
            b = x(0);
        }else 
		{
            a = x(1);
        }
        iter++;
    }
	
	return x;
	
}


