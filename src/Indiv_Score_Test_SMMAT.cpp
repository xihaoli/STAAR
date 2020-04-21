// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Indiv_Score_Test_SMMAT(arma::sp_mat G, arma::mat P, arma::vec residuals)
{
	int i;

	// number of markers
	int p = G.n_cols;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;

	arma::vec pvalue;
	pvalue.zeros(p);

	arma::vec Uscore_se;
	Uscore_se.zeros(p);

	double test_stat = 0;

	arma::mat Cov;
	Cov.zeros(p,p);

	Cov = trans(P*G)*G;

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

