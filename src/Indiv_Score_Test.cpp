// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Indiv_Score_Test(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, arma::vec residuals)
{
	int i;

	// number of markers
	int p = G.n_cols;
	int c = X.n_cols;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;

	arma::vec pvalue;
	pvalue.zeros(p);

	arma::vec Uscore_se;
	Uscore_se.zeros(p);

	double test_stat = 0;

	arma::mat tX_G;
	tX_G.zeros(c,p);


	arma::mat Cov;
	Cov.zeros(p,p);

	if(fam == 0)
	{
		tX_G = trans(X)*G;
		Cov = trans(G)*G - trans(tX_G)*inv(trans(X)*X)*tX_G;

	}else
	{
		tX_G = trans(X)*(arma::diagmat(working))*G;
		Cov = trans(arma::diagmat(working)*G)*G - trans(tX_G)*inv(trans(X)*arma::diagmat(working)*X)*tX_G;
	}

	for(i = 0; i < p; i++)
	{
		Uscore_se(i) = sqrt(Cov(i , i));

		if (Cov(i , i) == 0)
		{
			pvalue(i) = 1;
		}
		else
		{
			test_stat = pow(Uscore(i),2)/Cov(i,i);
			pvalue(i) = R::pchisq(test_stat,1,false,false);
		}

	}

	return List::create(Named("Uscore") = trans(Uscore), Named("Uscore_se") = Uscore_se, Named("pvalue") = pvalue);
}

