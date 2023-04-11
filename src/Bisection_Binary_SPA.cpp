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
double Bisection_Binary_SPA(arma::vec muhat, arma::vec G, double q, double xmin, double xmax, double tol)
{

    // the range of x to search
    double xupper = xmax;
    double xlower = xmin;
	
	double K1left = 0.0;
	double K1right = 0.0;
	
    double x0 = 0.0;
    double K1x0 = 1.0;
	
	K1left = K1_Binary_SPA(xlower, muhat, G, q);
	if((R_finite(K1left)==0)||(check_is_na(K1left)))
	{
		K1left = K1_Binary_SPA_alt(xlower, muhat, G, q);
	}
	
	K1right = K1_Binary_SPA(xupper, muhat, G, q);
	if((R_finite(K1right)==0)||(check_is_na(K1right)))
	{
		K1right = K1_Binary_SPA_alt(xupper, muhat, G, q);
	}
	
	if((!((R_finite(K1right)==0)||(check_is_na(K1right))))&&(!((R_finite(K1left)==0)||(check_is_na(K1left)))))
	{	
		if(!haveSameSign(K1left, K1right))
		{
			while ((fabs(xupper-xlower) > tol)&&(fabs(K1x0) > tol))
			{
				x0 = (xupper + xlower)/2.0;
				K1x0 = K1_Binary_SPA(x0, muhat, G, q);
				
				if((R_finite(K1x0)==0)||(check_is_na(K1x0)))
				{
					K1x0 = K1_Binary_SPA_alt(x0, muhat, G, q);
				}
			
				if(K1x0 == 0)
				{
					break;
				}else
				{
					if(haveSameSign(K1left, K1x0))
					{
						xlower = x0;
					}else
					{
						xupper = x0;
					}
				
					K1left = K1_Binary_SPA(xlower, muhat, G, q);
					if((R_finite(K1left)==0)||(check_is_na(K1left)))
					{
						K1left = K1_Binary_SPA_alt(xlower, muhat, G, q);
					}
	
					K1right = K1_Binary_SPA(xupper, muhat, G, q);
					if((R_finite(K1right)==0)||(check_is_na(K1right)))
					{
						K1right = K1_Binary_SPA_alt(xupper, muhat, G, q);
					}
				}
			}
		}		
	}
    return x0;
}
