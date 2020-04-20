// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// declare K
double K(double x, arma::vec egvalues);
// declare K2 (second derivative)
double K2(double x, arma::vec egvalues);
// declare bisection
double Bisection(arma::vec egvalues, double q, double xmin, double xmax);

// [[Rcpp::export]]
double Saddle(double q, arma::vec egvalues)
{
	bool lower = false;
	bool logp = false;
	double xmin = 0.0;
	double xmax = 0.0;
	double v = 0.0;
	double res = 0.0;

	int i;
	const int en = egvalues.size();
	double lambdamax = max(egvalues);
	q = q/lambdamax;

	for(i = 0; i < en; i++)
	{
		egvalues[i] = egvalues[i]/lambdamax;
	}
	lambdamax = 1.0;

    if (q > arma::sum(egvalues))
	{
        xmin = -0.01;
	}else
	{
        xmin = -en/(2 * q);
    }

    xmax = 1/(2*lambdamax) * 0.99999;

    double xhat = Bisection(egvalues,q,xmin,xmax);
    double w = sqrt(2*(xhat*q-K(xhat,egvalues)));
    if (xhat < 0){
        w = -w;
    }

	v = xhat*sqrt(K2(xhat,egvalues));
	if(fabs(xhat)<1e-04)
	{
		res = 2;
	}else
	{
		res =  R::pnorm(w+log(v/w)/w,0.0,1.0,lower,logp);
	}
    return res;


}

