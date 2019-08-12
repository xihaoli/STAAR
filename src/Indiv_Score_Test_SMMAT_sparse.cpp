// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List Indiv_Score_Test_SMMAT_sparse(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals)
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
	
	int q = Sigma_iX.n_cols;
	
	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(q,p);

	arma::mat Cov;
	Cov.zeros(p,p);
	
	tSigma_iX_G = trans(Sigma_iX)*G;
	Cov = trans(Sigma_i*G)*G - trans(tSigma_iX_G)*cov*tSigma_iX_G;

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















