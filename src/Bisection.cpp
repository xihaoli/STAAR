// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// declare K1 (first derivative)
double K1(double x, arma::vec egvalues, double q);

// [[Rcpp::export]]
double Bisection(arma::vec egvalues, double q, double xmin, double xmax)
{
   
    // the range of x to search
    double xupper = xmax;
    double xlower = xmin;
    
    double x0 = 0.0;
    double K1x0 = 1.0;
	
   
    while (fabs(xupper-xlower) > 1e-08)
    {
		x0 = (xupper + xlower)/2.0;
		K1x0 = K1(x0,egvalues,q);
		
		if(K1x0 == 0){
			break;
        }else if (K1x0 > 0){
            xupper = x0;
        }else{
            xlower = x0;
        }
    }
    
    return x0;
}

