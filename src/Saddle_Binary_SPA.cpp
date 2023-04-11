// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// declare K_Binary_SPA
double K_Binary_SPA(double x, arma::vec muhat, arma::vec G);
// declare K_Binary_SPA_alt (alternative way if not converge)
double K_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G);
// declare K1_Binary_SPA (first derivative)
double K1_Binary_SPA(double x, arma::vec muhat, arma::vec G, double q);
// declare K1_Binary_SPA_alt (first derivative, alternative way if not converge)
double K1_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G, double q);
// declare K2_Binary_SPA (second derivative)
double K2_Binary_SPA(double x, arma::vec muhat, arma::vec G);
// declare K2_Binary_SPA_alt (second derivative, alternative way if not converge)
double K2_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G);
// declare NR
double NR_Binary_SPA(arma::vec muhat, arma::vec G, double q, double init, double tol, int max_iter);
// check na
bool check_is_na(double x); 
// check same sign
bool haveSameSign(double a, double b);

// [[Rcpp::export]]
double Saddle_Binary_SPA(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, bool lower)
{
	bool logp = false;
	
	// init for NR
	double init = 0.0;

	double ki = 0.0;
	double res = 0.0;

	// Saddle point approximation
    double xhat = NR_Binary_SPA(muhat,G,q,init,tol,max_iter);
    double w = sqrt(2*(xhat*q - K_Binary_SPA(xhat,muhat,G)));
	
	if((R_finite(w)==0)||(check_is_na(w)))
	{
		w = sqrt(2*(xhat*q - K_Binary_SPA_alt(xhat,muhat,G)));
	}
	
    if (xhat < 0){
        w = -w;
    }

	ki = xhat*sqrt(K2_Binary_SPA(xhat,muhat,G));
	if((R_finite(ki)==0)||(check_is_na(ki)))
	{
		ki = xhat*sqrt(K2_Binary_SPA_alt(xhat,muhat,G));
	}
	
	if(fabs(xhat)<1e-04)
	{
		res = 1.0;
	}else
	{
		res =  R::pnorm(w+log(ki/w)/w,0.0,1.0,lower,logp);
	}
    
	return res;
}

