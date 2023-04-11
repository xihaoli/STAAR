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
// declare golden section search method for selecting initial value of bisection algorothm
arma::vec goldenSectionSearchForSignChange(double a, double b, arma::vec muhat, arma::vec G, double q, double tol, int max_iter);
// declare Bisection method
double Bisection_Binary_SPA(arma::vec muhat, arma::vec G, double q, double xmin, double xmax, double tol);
// check na
bool check_is_na(double x); 
// check same sign
bool haveSameSign(double a, double b);

// [[Rcpp::export]]
double Saddle_Binary_SPA_Bisection(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, double xmin, double xmax, bool lower)
{
	bool logp = false;
	
	double ki = 0.0;
	double res = 1.0;

	// Saddle point approximation
	// golden section search for selecting init value of bisection algorithm
	arma::vec xlimit;
	xlimit.ones(2);
	
	xlimit = goldenSectionSearchForSignChange(xmin, xmax, muhat, G, q, tol, max_iter);
	
	if(xlimit(0) < xlimit(1))
	{
		xmin = xlimit(0);
		xmax = xlimit(1);
	}else
	{
		xmin = xlimit(1);
		xmax = xlimit(0);
	}
	
    double xhat = Bisection_Binary_SPA(muhat, G, q, xmin, xmax, tol);
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



